! This file is part of Nagual software.
!
!    Nagual is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, version 3 of the License.
!
!    Nagual is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Nagual.  If not, see <http://www.gnu.org/licenses/>.
!
! Cite as follows:
! 
! R. Flores-Moreno, H. N. Gonzalez-Ramirez, J. F. H. Lew-Yee, J. M. del Campo,
! J. J. Villalobos-Castro, J. A. Flores-Ramos, J. A. Guerrero-Cruz,
! B. Zuniga-Gutierrez, Nagual 1, Guadalajara Jal., Mexico (2020)
!
!###################################################
!   Nagual: Multicomponent many body calculations.
!   Copyright (C) 2006-2020 Nagual developers.
!
!   Contact: r.flores@academicos.udg.mx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nmscf
! Self-Consistent Field: RHF, UHF
  
  use nmparameter
  use nmtypes
  use nmmatrix
  use nmintegrals
  use nmfile
  use nmstring
  use nmshell
  use nmatom
  use nmmolecule
  use nmtime
  use nmsystem
  use nmbasis
  use nmaux
  use nmunits
  use nmadft
  use nmxf
  use nmcf
  use nmecp

  implicit none

    public :: scf
    public :: scf_build_density_matrix
    public :: scf_build_pulay_matrix

    private

      integer :: tape

contains

  logical function scf(task,sys)
  ! Self Consistent Field Iterations
  ! Roberto Flores Moreno,  2008, 2009, 2014, 2018
  !
  ! C. C. J. Roothaan, 
  ! Rev. Mod. Phys. 23, 69 (1951).
  !
  ! G. G. Hall,
  ! Proc. Roy. Soc. Ser. A 205, 541 (1951).
  !
  ! J. A. Pople and R. K. Nesbet,
  ! J. Chem. Phys. 22, 571 (1954).
  !
  ! A. Szabo and N. S. Ostlund,
  ! Modern Quantum Chemistry
  ! (Dover Pubs. Inc., Mineola, 1996).
  !
  ! R. Flores-Moreno, 
  ! PhD Thesis (Cinvestav-IPN, Mexico, 2006)
  implicit none
    type(nscftask) :: task
    type(nsystem) :: sys

    integer :: i,is,k,js,nit
    logical :: converged,dft,gotz
    real(8) :: alpha,beta,coree,coret,error,energy,exc,penergy,&
               nucree,scfmix,tf
    real(8) :: eig(nmaxnbas)

    integer :: allocs,dabas,dbas,dbas2,dorb
    real(8),allocatable :: s(:,:),w(:,:)
    real(8),allocatable :: f(:,:),p(:,:),p2(:,:)
    real(8),allocatable :: g(:,:),x(:),z(:)

    call file_open_ascii("Nagual.scf",tape)
    call file_copyright(tape)
    call file_header(tape,"Self-Consistent Field")

    call system_print(tape,sys)
    call file_string(tape,1,'SCF',1)

    ! Nuclear repulsion energy
    call getnre(sys%mol,nucree)

    do is=1,sys%ns
      dbas = basis_nfun(sys%basis(is))
      if (sys%orbital_type.ne.'cartesian') then
        dorb = basis_nfun_spherical(sys%basis(is))
      else
        dorb = dbas
      end if

      allocate(s(dbas,dbas),w(dbas,dbas),stat=allocs)
      if (allocs.gt.0) call file_error('SCF: Allocation failure')

      ! Evaluate overlap matrix 
      call build_overlap_matrix(sys%mol,sys%basis(is),s,dbas)
      write(tape,'(t2,"Overlap matrix for subsystem ",I2)') is 
      call matrix_print(tape,s,dbas,dbas,1,1,dbas,dbas)
      call file_matrix_bas(s,dbas,is,'overlap','write')

      w = s
      call matrix_svd_power(w,eig,dbas,0.5,1.0e-14)
      call file_matrix_bas(w,dbas,is,'overlaph','write')

      ! Symmetric orthogonalization matrix 
      w = s
      if (sys%orbital_type.ne.'cartesian') then
        call cartesian_to_spherical(sys%basis(is),w,s,dbas)
      end if
      call matrix_svd_power(s(1:dorb,1:dorb),eig,dorb,-0.5,1.0e-14)
      !rfm call build_orthogonalization_matrix(s,dbas,dorb,-0.5)
      write(tape,'(t2,"X matrix for subsystem ",I2)') is
      call matrix_print(tape,s,dbas,dbas,1,1,dbas,dbas)
      call file_matrix_bas(s,dbas,is,'ortho','write')

      ! Nuclear atraction integrals
      call build_core_matrix(sys,sys%basis(is),w,dbas,is)
      call file_matrix_bas(w,dbas,is,'core','write')
      write(tape,'(t2,"Core matrix for subsystem ",I2)') is
      call matrix_print(tape,w,dbas,dbas,1,1,dbas,dbas)

      ! Initialize Fock matrix
      call file_matrix_bas(w,dbas,is,'fock','write')

      ! Initialize density matrices
      w(:,:) = 0.0
      call file_matrix_bas(w,dbas,is,'density','write')

      call doguess(sys,w,s,dbas,is,task%guess)

      deallocate(s,w,stat=allocs)
      if (allocs.gt.0) call file_error('SCF: Deallocation failure')
    end do

    ! Prepare G
    if (task%eris.eq.3) call aux_setup(tape,sys)

    dft = .false.
    if (xf_jacob(sys)+cf_jacob(sys).gt.0) dft = .true.

    call file_header(tape,"SCF loop")

    exc = 0.0
    energy = 0.0
    error = 10.0
    nit = 0
    converged = .false.
    do while (.not.converged)
      nit = nit + 1

      do is=1,sys%ns
        dbas = basis_nfun(sys%basis(is))
        if (sys%orbital_type.ne.'cartesian') then
          dorb = basis_nfun_spherical(sys%basis(is))
        else
          dorb = dbas
        end if

        allocate(f(dbas,dbas),p(dbas,dbas),stat=allocs)
        if (allocs.gt.0) call file_error('SCF: Allocation failure')

        if (nit.gt.1) then
          call file_matrix_bas(f,dbas,is,'orbitals','read')
          call scf_build_density_matrix(sys,f,p,dbas,is)
          call file_matrix_bas(f,dbas,is,'density','read')
          p = (1.0-task%damping)*f + task%damping*p
          call file_matrix_bas(p,dbas,is,'density','write')
        end if

        deallocate(f,p,stat=allocs)
        if (allocs.gt.0) call file_error('SCF: Deallocation failure')
      end do

      ! Build Fock matrix
      gotz = .false.
      if ((task%guess.ne.'core').or.(nit.gt.1)) then
        if (task%eris.eq.3) call vfcp(tape,sys)
        call adft_vxc(tape,sys,exc,dft)
        if (dft) gotz = .true.
        do is=1,sys%ns
          dbas = basis_nfun(sys%basis(is))
          allocate(f(dbas,dbas),stat=allocs)
          if (allocs.gt.0) call file_error('SCF: Allocation failure')
          call build_fock_matrix(tape,sys,is,f,dbas,nit,task%eris,dft)
          deallocate(f,stat=allocs)
          if (allocs.gt.0) call file_error('SCF: Deallocation failure')
        end do
      end if

      ! energy
      penergy = energy

      coree = 0.0
      energy = exc
      do is=1,sys%ns
        dbas = basis_nfun(sys%basis(is))
        allocate(f(dbas,dbas),p(dbas,dbas),stat=allocs)
        if (allocs.gt.0) call file_error('SCF: Allocation failure')
        call file_matrix_bas(p,dbas,is,'density','read')
        call file_matrix_bas(f,dbas,is,'fock','read')
        energy = energy + 0.5*sum(p*f)
        call file_matrix_bas(f,dbas,is,'core','read')
        coret = sum(p*f)
        coree = coree + coret
        energy = energy + 0.5*coret
        deallocate(f,p,stat=allocs)
        if (allocs.gt.0) call file_error('SCF: Deallocation failure')
        if (gotz) then
          dabas = basis_nfun(sys%auxis(is))
          allocate(g(dabas,dabas),x(dabas),z(dabas),stat=allocs)
          if (allocs.gt.0) call file_error('SCF: Allocation failure')
          call file_matrix_aux(g,dabas**2,is,is,'g','read')
          call file_vector_aux(z,dabas,is,0,'x','read')
          x = matmul(g,z)
          call file_vector_aux(z,dabas,is,0,'z','read')
          energy = energy - 0.5*sum(x(1:dabas)*z(1:dabas))
          deallocate(g,x,z,stat=allocs)
          if (allocs.gt.0) call file_error('SCF: Deallocation failure')
        end if
      end do
      energy = energy + nucree
      sys%energy = energy

      ! Check convergence on energy
      error = penergy - energy
      write(tape,'(t2,i4," Energy: ",f24.8," Error: ",f24.8)')&
        nit,energy,error

      if (nit.gt.1) then
        if (abs(error).lt.task%tol) converged = .true.
      end if
      call file_flush(tape)
      if (nit.gt.task%nitmax) go to 2000

      ! MOs
      do is=1,sys%ns
        dbas = basis_nfun(sys%basis(is))
        if (sys%orbital_type.ne.'cartesian') then
          dorb = basis_nfun_spherical(sys%basis(is))
        else
          dorb = dbas
        end if
        allocate(f(dbas,dbas),p(dbas,dbas),w(dbas,dbas),stat=allocs)
        if (allocs.gt.0) call file_error('SCF: Allocation failure')
        call scf_diis(f,p,w,dbas,nit,task%diis,is)
        call build_orbital_matrix(sys,f,p,w,dbas,dorb,is) 
        deallocate(f,p,w,stat=allocs)
        if (allocs.gt.0) call file_error('SCF: Deallocation failure')
      end do

    end do
 2000 continue

    write(tape,'(/,t2," Core energy (a.u.): ",f24.8)') coree
    write(tape,'(t2," Core energy (eV): ",f24.8)') units_au_to_ev(coree)

    if (converged) then
      call file_string(tape,1,' SCF CONVERGED ',0)
    else
      call file_string(tape,1,' SCF NOT CONVERGED ',0)
    end if
    call file_flush(tape)

    call print_orbitals(tape,sys)
    call file_close(tape)

    scf = converged
  end function

  subroutine build_fock_matrix(tape,sys,is,f,dbas,nit,eris,dft)
  ! Build AO Fock matrix (F)
  ! Roberto Flores Moreno, 2008, 2009, 2018
  implicit none
    integer :: dbas,eris,is,nit,tape
    logical :: dft
    real(8) :: f(dbas,dbas)
    type(nsystem) :: sys

    logical :: found
    integer :: js

    ! Search for built copy (RHF)
    found = .false.
    do js=1,is-1
      if (sys%species_name(js).eq.sys%species_name(is)) then
        if (sys%basis(js)%filename.eq.sys%basis(is)%filename) then
          if (sys%nmoo(js).eq.sys%nmoo(is)) then
            call file_matrix_bas(f,dbas,js,'fock','read')
            found = .true.
          end if
        end if
      end if
    end do

    if (.not.found) then
      call file_matrix_bas(f,dbas,is,'core','read')
      if (eris.eq.4) then
        call build_fock_matrix_dir4(tape,sys,is,f,dbas)
      else
        call build_fock_matrix_dir3(tape,sys,is,f,dbas,dft)
      end if
    end if

    call file_matrix_bas(f,dbas,is,'fock','write')
  end subroutine

  subroutine build_fock_matrix_dir4(tape,sys,is,f,dbas)
  ! Build AO Fock matrix (F)
  ! Roberto Flores Moreno, 2008, 2009, 2018
  implicit none
    type(nsystem) :: sys
    integer :: dbas,is,tape
    real(8) :: f(dbas,dbas)

    integer :: dshell,i,iset,ishell,js,jset,jshell,ks,kset,kshell,lset,lshell
    integer :: ba,bb,bc,bd,na,nb,nc,nd,ibas,jbas,kbas,lbas,ii,jj,kk,ll
    logical :: gocoul,goex,skip(2)
    logical :: sused(nmaxns)
    real(8) :: wex
    real(8) :: ra(3),rb(3),rc(3),rd(3)
    real(8) :: blk(nmaxnco,nmaxnco,nmaxnco,nmaxnco)
    type(nshell) :: sa,sb,sc,sd

    integer :: allocs,dbas2
    real(8),allocatable :: p(:,:),w(:,:)

    sused(1:sys%ns) = .false.
    do js=1,sys%ns
      dbas2 = basis_nfun(sys%basis(js))
      allocate(p(dbas2,dbas2),w(dbas2,dbas2),stat=allocs)
      if (allocs.gt.0) call file_error('fock: Allocation failure')
      gocoul = .false.
      if (.not.sused(js)) then
        gocoul = .true.
        sused(js) = .true.
        call file_matrix_bas(p,dbas2,js,'density','read')
        do ks=js+1,sys%ns
          if (sys%species_name(js).eq.sys%species_name(ks)) then
            if (sys%basis(js)%filename.eq.sys%basis(ks)%filename) then
              allocate(w(dbas2,dbas2),stat=allocs)
              call file_matrix_bas(w,dbas2,ks,'density','read')
              p(1:dbas2,1:dbas2) = p(1:dbas2,1:dbas2) +  w(1:dbas2,1:dbas2) 
              sused(ks) = .true.
            end if
          end if
        end do
      end if
      if (js.eq.is) then
        wex = sys%interaction(is,is)%exchange(1)
        if (wex.eq.0.0) then
          goex = .false.
        else
          goex = .true.
          call file_matrix_bas(w,dbas,is,'density','read')
          if (wex.ne.1.0) w(:,:) = wex*w(:,:)
        end if
      else
          goex = .false.
      end if
              
    if (gocoul.or.goex) then

    call integrals_set_ppi_type(sys%interaction(is,js)%simple)
    ba = 0
    ii = 0
    do iset=1,sys%basis(is)%nsets
      ra(1:3) = sys%mol%atom(sys%basis(is)%set(iset)%atom)%pos(1:3)
      do ishell=1,sys%basis(is)%set(iset)%nshell
        ii = ii + 1
        sa = sys%basis(is)%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        bb = 0
        jj = 0
        do jset=1,sys%basis(is)%nsets
          rb(1:3) = sys%mol%atom(sys%basis(is)%set(jset)%atom)%pos(1:3)
          do jshell=1,sys%basis(is)%set(jset)%nshell
            jj = jj + 1
            sb = sys%basis(is)%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            if ((jset.gt.iset).or.                                    &
     &          ((jset.eq.iset).and.(jshell.ge.ishell))) then
              bc = 0
              kk = 0
              do kset=1,sys%basis(js)%nsets
                rc(1:3) = sys%mol%atom(sys%basis(js)%set(kset)%atom)%pos(1:3)
                do kshell=1,sys%basis(js)%set(kset)%nshell
                  kk = kk + 1
                  sc = sys%basis(js)%set(kset)%shell(kshell)
                  nc = (sc%l+1)*(sc%l+2)/2
                  bd = 0
                  ll = 0
                  do lset=1,sys%basis(js)%nsets
                    rd(1:3) = sys%mol%atom(sys%basis(js)%set(lset)%atom)%pos(1:3)
                    do lshell=1,sys%basis(js)%set(lset)%nshell
                      ll = ll + 1
                      sd = sys%basis(js)%set(lset)%shell(lshell)
                      nd = (sd%l+1)*(sd%l+2)/2
                      if ((lset.gt.kset).or.  &
     &                    ((lset.eq.kset).and.(lshell.ge.kshell))) then
                        call integrals_ppi(ra,rb,rc,rd,&
                                           sa,sb,sc,sd,blk,ntolint,skip,0,0,0)
                        if (skip(1)) go to 1000
                        if (.not.skip(2)) then
                          if (gocoul) call ppi_sum(ba,bb,bc,bd,na,nb,&
                                      nc,nd,f,p,blk,dbas,dbas2)
                          if (goex) call exchange_sum(ba,bb,bc,bd,na,&
                                                  nb,nc,nd,f,w,blk,dbas)
                        end if
                      end if
                      bd = bd + nd
                    end do
                  end do
                  bc = bc + nc
                end do
              end do
            end if
 1000       continue
            bb = bb + nb
          end do
        end do
        ba = ba + na
      end do
    end do

    end if


      deallocate(p,w,stat=allocs)
      if (allocs.gt.0) call file_error('fock: Deallocation failure')
    end do

    ! Symmetrize
    call matrix_symmetrize(f,dbas,dbas,'uplow')
    
  end subroutine 

  subroutine ppi_sum(f1,f2,f3,f4,n1,n2,n3,n4,f,p,b,dbas,dbas2)
  ! Add up four center integrals to Fock matrix (Coulomb)
  ! Roberto Flores-Moreno, 2008, 2012, 2018
  implicit none
    integer :: f1,f2,f3,f4,n1,n2,n3,n4,dbas,dbas2
    real(8) :: f(dbas,dbas),p(dbas2,dbas2)
    real(8) :: b(nmaxnco,nmaxnco,nmaxnco,nmaxnco)

    integer :: i,is,j,js,k,l,b1,b2,b3,b4
    real(8) :: s
    real(8) :: pblk(nmaxnco,nmaxnco)

    do k=f3+1,f3+n3
      b3 = k - f3
      do l=max(f4+1,k),f4+n4
        b4 = l - f4
        pblk(b3,b4) = p(k,l)
        if (l.ne.k) pblk(b3,b4) = 2.0*pblk(b3,b4)
      end do
    end do

    do i=f1+1,f1+n1
      b1 = i - f1
      do j=max(f2+1,i),f2+n2
        b2 = j - f2
        s = 0.0
        do k=f3+1,f3+n3
          b3 = k - f3
          do l=max(f4+1,k),f4+n4
            b4 = l - f4
            s = s + pblk(b3,b4)*b(b1,b2,b3,b4)
          end do
        end do
        f(i,j) = f(i,j) + s
      end do
    end do
  end subroutine 

  subroutine exchange_sum(f1,f2,f3,f4,n1,n2,n3,n4,f,p,b,dbas)
  ! Add up four center integrals to Fock matrix (Exchange)
  ! Roberto Flores-Moreno, 2008, 2012, 2018
  implicit none
    integer :: f1,f2,f3,f4,n1,n2,n3,n4,dbas
    real(8) :: f(dbas,dbas),p(dbas,dbas)
    real(8) :: b(nmaxnco,nmaxnco,nmaxnco,nmaxnco)

    integer :: i,j,k,l,b1,b2,b3,b4
    real(8) :: s,w

    do i=f1+1,f1+n1
      b1 = i - f1
      do l=max(f4+1,i),f4+n4
        b4 = l - f4
        s = 0.0
        do j=max(f2+1,i),f2+n2
          b2 = j - f2
          do k=f3+1,min(l,f3+n3)
            b3 = k - f3
            s = s + p(k,j)*b(b1,b2,b3,b4)
          end do
        end do
        f(i,l) = f(i,l) - s
      end do
    end do

    do i=f2+1,f2+n2
      b2 = i - f2
      do l=max(f4+1,i),f4+n4
        b4 = l - f4
        s = 0.0
        do j=f1+1,min(f1+n1,i-1)
          b1 = j - f1
          do k=f3+1,min(l,f3+n3)
            b3 = k - f3
            s = s + p(k,j)*b(b1,b2,b3,b4)
          end do
        end do
        f(i,l) = f(i,l) - s
      end do
    end do

    do i=f1+1,f1+n1
      b1 = i - f1
      do l=max(f3+1,i),f3+n3
        b3 = l - f3
        s = 0.0
        do j=max(f2+1,i),f2+n2
          b2 = j - f2
          do k=max(l+1,f4+1),f4+n4
            b4 = k - f4
            s = s + p(k,j)*b(b1,b2,b3,b4)
          end do
        end do
        f(i,l) = f(i,l) - s
      end do
    end do

    do i=f2+1,f2+n2
      b2 = i - f2
      do l=max(f3+1,i),f3+n3
        b3 = l - f3
        s = 0.0
        do j=f1+1,min(i-1,f1+n1)
          b1 = j - f1
          do k=max(l+1,f4+1),f4+n4
            b4 = k - f4
            s = s + p(k,j)*b(b1,b2,b3,b4)
          end do
        end do
        f(i,l) = f(i,l) - s
      end do
    end do
  end subroutine 

  subroutine scf_build_density_matrix(sys,c,p,dbas,is)
  ! Build AO density matrix (P)
  ! Roberto Flores Moreno 2008-2012, 2018
  implicit none
    type(nsystem) :: sys
    integer :: dbas,is
    real(8) :: c(dbas,dbas),p(dbas,dbas)

    integer :: i,j,mu,nu

    do mu=1,dbas
      do nu=mu,dbas
        p(mu,nu) = 0.0
        do i=1,sys%nmoo(is)
          p(mu,nu) = p(mu,nu) + sys%moocc(i,is)*c(mu,i)*c(nu,i)
        end do
        p(nu,mu) = p(mu,nu) 
      end do
    end do
  end subroutine 

  subroutine scf_build_pulay_matrix(sys,c,w,dbas,is)
  ! Build AO energy weighted density matrix (W)
  ! Roberto Flores Moreno, 2008-2012, 2018
  implicit none
    type(nsystem) :: sys
    integer :: dbas,is
    real(8) :: c(dbas,dbas),w(dbas,dbas)

    integer :: i,j,mu,nu

    do mu=1,dbas
      do nu=mu,dbas
        w(mu,nu) = 0.0
        do i=1,sys%nmoo(is)
          w(mu,nu) = w(mu,nu) + sys%moocc(i,is)*sys%moe(i,is)*c(mu,i)*c(nu,i)
        end do
        w(nu,mu) = w(mu,nu)
      end do
    end do
  end subroutine

  subroutine build_orthogonalization_matrix(x,dx,dorb)
  ! Build canonical orthogonalization (X) matrix
  ! Roberto Flores Moreno, 2008
  implicit none
    integer :: dorb,dx
    real(8) :: x(dx,dx)

    integer :: i,j
    real(8) :: eigv(nmaxnbas)

    call matrix_diagonalize(x,eigv,dx,dorb)

    do i=1,dorb
      if (eigv(i).le.ntolnum) then
        eigv(i) = 0.0
        do j=1,dorb
          x(j,i) = 0.0
        end do
      else
        do j=1,dorb
          x(j,i) = x(j,i)/sqrt(eigv(i))
        end do
      end if
    end do
  end subroutine 

  subroutine build_overlap_matrix(m,b,s,ds)
  ! Build overlap (S) matrix
  ! Roberto Flores Moreno, 2008, 2018
  implicit none
    type(nmolecule) :: m
    type(nbasis) :: b
    integer :: ds
    real(8) :: s(ds,ds)

    integer :: iset,ibas,ishell,jset,jbas,jshell,na,nb,fa,fb
    real(8) :: ra(3),rb(3)
    real(8) :: blk(nmaxnco,nmaxnco)
    type(nshell) :: sa,sb

    fa = 0
    do iset=1,b%nsets
      ra(1:3) = m%atom(b%set(iset)%atom)%pos(1:3)
      do ishell=1,b%set(iset)%nshell
        sa = b%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        fb = 0
        do jset=1,b%nsets
          rb(1:3) = m%atom(b%set(jset)%atom)%pos(1:3)
          do jshell=1,b%set(jset)%nshell
            sb = b%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            if ((jset.gt.iset).or.                                    &
     &          ((jset.eq.iset).and.(jshell.ge.ishell))) then
              call integrals_overlap(ra,rb,sa,sb,blk,nmaxnco,0)
              do ibas=1,na
                do jbas=1,nb
                  s(fa+ibas,fb+jbas) = blk(ibas,jbas)
                end do
              end do
            end if
            fb = fb + nb
          end do
        end do
        fa = fa + na
      end do
    end do

    ! Symmetrize
    call matrix_symmetrize(s,ds,ds,'uplow')
  end subroutine

  subroutine build_core_matrix(sys,b,h,dbas,is)
  ! Build mono-electronic hamiltonian matrix including:
  !  - kinetic energy matrix
  !  - nuclear atraction matrix
  !  - ECPs
  !  - MCPs
  ! Roberto Flores Moreno, 2008, 2018
  implicit none
    type(nsystem) :: sys
    type(nbasis) :: b
    integer :: dbas,is
    real(8) :: h(dbas,dbas)

    integer :: i,iset,ibas,ishell,jset,jbas,jshell,katom,na,nb,fa,fb
    real(8) :: ra(3),rb(3),rc(3)
    real(8) :: blk(nmaxnco,nmaxnco),sblk(nmaxnco,nmaxnco)
    real(8) :: dblk(nmaxnco,nmaxnco,3)
    type(nshell) :: sa,sb

    fa = 0
    do iset=1,b%nsets
      ra(1:3) = sys%mol%atom(b%set(iset)%atom)%pos(1:3)
      do ishell=1,b%set(iset)%nshell
        sa = b%set(iset)%shell(ishell)
        na = (sa%l+1)*(sa%l+2)/2
        fb = 0
        do jset=1,b%nsets
          rb(1:3) = sys%mol%atom(b%set(jset)%atom)%pos(1:3)
          do jshell=1,b%set(jset)%nshell
            sb = b%set(jset)%shell(jshell)
            nb = (sb%l+1)*(sb%l+2)/2
            if ((jset.gt.iset).or.                                    &
     &          ((jset.eq.iset).and.(jshell.ge.ishell))) then
              ! Nuclear atraction integrals
              sblk(:,:) = 0.0
              do katom=1,sys%mol%natom
                if (sys%mol%atom(katom)%zeff.ne.0.0) then
                  rc(1:3) = sys%mol%atom(katom)%pos(1:3)
                  call integrals_core(ra,rb,rc,sa,sb,blk,0,0)
                  sblk(1:na,1:nb) = sblk(1:na,1:nb) + blk(1:na,1:nb)*     &
     &            sys%mol%atom(katom)%zeff*sys%species_charge(is)
                end if
              end do
              ! Save block
              do ibas=1,na
                do jbas=1,nb
                   h(fa+ibas,fb+jbas) = sblk(ibas,jbas)
                end do
              end do
              ! Kinetic energy integrals
              sblk(:,:) = 0.0
              call integrals_kinetic(ra,rb,sa,sb,sblk,0)
              do ibas=1,na
                do jbas=1,nb
                  h(fa+ibas,fb+jbas) = h(fa+ibas,fb+jbas) + &
                  sblk(ibas,jbas)/sys%species_mass(is)
                end do
              end do

            end if
            fb = fb + nb
          end do
        end do
        fa = fa + na
      end do
    end do

    ! Effective core potentials (only for electrons)
    !rfm if (string_to_lowercase(sys%species_name(is)).eq.'electron') &
    !rfm  call ecp(sys%mol,sys%ecps,sys%basis(is),h,dbas)

    ! Symmetrize
    call matrix_symmetrize(h,dbas,dbas,'uplow')
  end subroutine 

  subroutine getnre(m,nre)
  ! Get nuclear repulsion energy
  ! Roberto Flores-Moreno, 2008
  implicit none
    type(nmolecule) :: m
    real(8) :: nre

    integer :: iatom,jatom
    real(8) :: za,zb,ra(3),rb(3)

    nre = 0.0
    do iatom=1,m%natom
      ra = m%atom(iatom)%pos
      za = m%atom(iatom)%zeff
      if (za.ne.0.0) then
        do jatom=iatom+1,m%natom
          rb = m%atom(jatom)%pos
          zb = m%atom(jatom)%zeff
          if (zb.ne.0.0) then
            nre = nre + za*zb/sqrt(sum((ra-rb)**2))
          end if
        end do
      end if
    end do
  end subroutine getnre

  subroutine build_orbital_matrix(sys,f,c,w,dbas,dorb,is)
  ! Build MO coefficient matrix
  ! Roberto Flores-Moreno, 2008, 2018
  implicit none
    type(nsystem) :: sys
    integer :: dbas,dorb,is
    real(8) :: f(dbas,dbas),c(dbas,dbas),w(dbas,dbas)

    integer :: ba,i,iatom,iao,imo,ii,ishell,j,jmax,na
    real(8) :: delta,smax

    call file_matrix_bas(f,dbas,is,'fock','read')

    ! Cartesian to spherical
    if (sys%orbital_type.ne.'cartesian') then
      c = f
      call cartesian_to_spherical(sys%basis(is),c,f,dbas)
    end if

    ! Transform Fock matrix
    call file_matrix_bas(w,dbas,is,'ortho','read')
    call matrix_transform(f(1:dorb,1:dorb),w(1:dorb,1:dorb),&
          c(1:dorb,1:dorb),dorb)

    ! Diagonalize
    call matrix_diagonalize(f,sys%moe(1,is),dbas,dorb)

    call file_matrix_bas(f,dbas,is,'oorbitals','write')

    ! Orthogonal -> non-orthogonal
    call file_matrix_bas(w,dbas,is,'ortho','read')
    c(1:dbas,1:dorb) = matmul(w(1:dbas,1:dbas),f(1:dbas,1:dorb))

    ! Transform back MO coefficients to cartesian representation
    if (sys%orbital_type.ne.'cartesian') then
      w = c
      call spherical_to_cartesian(sys%basis(is),w,c,dbas,dorb)
    end if

    call file_matrix_bas(c,dbas,is,'orbitals','write')
  end subroutine 
!
  subroutine cartesian_to_spherical(b,car,sph,dbas)
  ! Transform Fock matrix from the cartesian to the spherical basis
  ! Roberto Flores Moreno, Apr 2009
  !
  ! H. B. Schlegel and M. J. Frisch,
  ! Int. J. Quantum Chem. 54, 83 (1995).
  !
  implicit none
    integer :: dbas
    real(8) :: car(dbas,dbas),sph(dbas,dbas)
    type(nbasis) :: b

    integer :: aica,aicb,aisa,aisb,cfa,cfb,cna,cnb,iset,ica,icb,isa,isb
    integer :: ishell,jset,jshell,sfa,sfb,sna,snb
    type(nshell) :: sa,sb

    ! Initialize
    sph(:,:) = 0.0

    cfa = 0
    sfa = 0
    do iset=1,b%nsets
      do ishell=1,b%set(iset)%nshell
        sa = b%set(iset)%shell(ishell)
        cna = (sa%l+1)*(sa%l+2)/2
        sna = 2*sa%l+1
        cfb = 0
        sfb = 0
        do jset=1,b%nsets
          do jshell=1,b%set(jset)%nshell
            sb = b%set(jset)%shell(jshell)
            cnb = (sb%l+1)*(sb%l+2)/2
            snb = 2*sb%l+1
            do isa=1,sna
              aisa = isa + sfa
              do isb=1,snb
                aisb = isb + sfb
                do ica=1,cna
                  aica = ica + cfa
                  do icb=1,cnb
                    aicb = icb + cfb
                    sph(aisa,aisb) = sph(aisa,aisb) + car(aica,aicb)*   &
     &              ctostm(isa,ica,sa%l)*ctostm(isb,icb,sb%l)
                  end do
                end do
              end do
            end do
            cfb = cfb + cnb
            sfb = sfb + snb
          end do
        end do
        cfa = cfa + cna
        sfa = sfa + sna
      end do
    end do
  end subroutine

  subroutine spherical_to_cartesian(b,sph,car,dbas,dorb)
  ! Transform MO coefficients from the spherical AO basis to the cartesian
  ! Roberto Flores Moreno, Apr 2009
  !
  ! H. B. Schlegel and M. J. Frisch,
  ! Int. J. Quantum Chem. 54, 83 (1995).
  !
  implicit none
    integer :: dbas,dorb
    real(8) :: sph(dbas,dbas),car(dbas,dbas)
    type(nbasis) :: b

    integer :: aica,aisa,cfa,cna,iset,ica,imo,isa,ishell,sfa,sna
    type(nshell) :: sa

    ! Initialize 
    car(:,:) = 0.0

    do imo=1,dorb
      cfa = 0
      sfa = 0
      do iset=1,b%nsets
        do ishell=1,b%set(iset)%nshell
          sa = b%set(iset)%shell(ishell)
          cna = (sa%l+1)*(sa%l+2)/2
          sna = 2*sa%l+1
          do isa=1,sna
            aisa = isa + sfa
            do ica=1,cna
              aica = ica + cfa
              car(aica,imo) = car(aica,imo) + ctostm(isa,ica,sa%l)*     &
     &                        sph(aisa,imo)
            end do
          end do
          cfa = cfa + cna
          sfa = sfa + sna
        end do
      end do
    end do
  end subroutine

  subroutine doguess(sys,c,p,dbas,is,guess)
  ! Do special guess
  ! Roberto Flores Moreno, Nov 2008, 2018
  implicit none
    character*(*) :: guess
    integer :: dbas,is
    real(8) :: c(dbas,dbas),p(dbas,dbas)
    type(nsystem) :: sys

    !rfm logical :: loaded(nmaxatom)

    ! Restart
    if (guess.eq.'restart') then
      ! Make sure to copy old MOS file
      call file_matrix_bas(c,dbas,is,'orbitals','read')
      call scf_build_density_matrix(sys,c,p,dbas,is)
    ! Core
    else
      c = 0.0
      p = 0.0
    end if

    call file_matrix_bas(c,dbas,is,'orbitals','write')
    call file_matrix_bas(p,dbas,is,'density','write')
  end subroutine

  subroutine print_orbitals(tape,sys)
  ! Print MOs to tape
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: tape
    type(nsystem) :: sys

    integer :: i,is

    integer :: allocs,dbas,dorb
    real(8),allocatable :: c(:,:)

    write(tape,'(//,t2,"Orbital energies and coefficients",/)') 
    do is=1,sys%ns
      dbas = basis_nfun(sys%basis(is))
      if (sys%orbital_type.ne.'cartesian') then
        dorb = basis_nfun_spherical(sys%basis(is))
      else
        dorb = dbas
      end if
      write(tape,'(/,t2,"Set:",i2,/)') is
      write(tape,'(/,t2,"Number of MOs: ",i5,/)') dorb
      write(tape,'(/,t2,"Energies for set ",i3,/)') is
      do i=1,dorb
        write(tape,'(t2,i8,5x,f5.2,5x,f25.13,x,a2)') i,sys%moocc(i,is),&
        sys%moe(i,is)
      end do
      allocate(c(dbas,dbas),stat=allocs)
      if (allocs.gt.0) &
        call file_error('scf/print_orbitals: Allocation failure')
      write(tape,'(/,t2,"Coefficients for set ",i3,/)') is
      call file_matrix_bas(c,dbas,is,'orbitals','read')
      call matrix_print(tape,c,dbas,dbas,1,1,dbas,dorb)
      deallocate(c,stat=allocs)
      if (allocs.gt.0) &
        call file_error('scf/print_orbitals: Deallocation failure')
    end do
  end subroutine

  subroutine scf_diis(f,w1,w2,dbas,nit,nmat,is)
  ! Perform DIIS for SCF
  ! Felipe Lew-Yee, Jorge Martin del Campo-Ramirez, R. Flores-Moreno, 2019
  implicit none
    integer :: dbas,is,nit,nmat
    real(8) :: f(dbas,dbas),w1(dbas,dbas),w2(dbas,dbas)

    integer :: i,info,j,nf,nitm

    integer :: allocs,db
    integer,allocatable :: ipiv(:)
    real(8),allocatable :: b(:,:),c(:)

    if ((nmat.lt.0).or.(nit.le.1)) return
    nitm = nit - 1

    call file_matrix_bas(f,dbas,is,'fock','read')
    call file_matrix_diis(f,dbas,is,'f','write',nitm,nmat)

    call file_matrix_bas(w1,dbas,is,'density','read')
    f = matmul(f,w1)

    call file_matrix_bas(w1,dbas,is,'overlap','read')
    f = matmul(f,w1)

    call file_matrix_bas(w1,dbas,is,'ortho','read')
    f = matmul(transpose(w1),matmul(f,w1))

    f = f - transpose(f)

    call file_matrix_diis(f,dbas,is,'r','write',nitm,nmat)

    if (nit.lt.3) return

    !FLY: If not enough saved matrices, or if use all matrices option
    if ((nmat.ge.nitm).or.(nmat.eq.0)) then
      nf = nitm
    !FLY: Use the number of fock matrices requested in input
    else if (nmat.lt.nitm) then
      nf = nmat
    end if
    db = nf + 1
    allocate(b(db,db),c(db),ipiv(db),stat=allocs)
    if (allocs.gt.0) call file_error('scf_diis: allocation')

    b(db,1:nf) = -1.0
    b(1:nf,db) = -1.0
    b(db,db) = 0.0
    do i=1,nf
      call file_matrix_diis(w1,dbas,is,'r','read',i,nmat)
      do j=1,nf
        call file_matrix_diis(w2,dbas,is,'r','read',j,nmat)
        b(i,j) = sum(w1*w2)
      end do
    end do

    c(1:nf) = 0.0
    c(db) = -1.0
    call dgesv(db,1,b,db,ipiv,c,db,info)
    if (info.ne.0) call file_error('scf_diis: dgesv ')

    f(:,:) = 0.0
    do i=1,nf
      call file_matrix_diis(w1,dbas,is,'f','read',i,nmat)
      f(:,:) = f(:,:) + c(i)*w1(:,:)
    end do

    call file_matrix_bas(f,dbas,is,'fock','write')

    deallocate(b,c,ipiv,stat=allocs)
    if (allocs.gt.0) call file_error('scf_diis: deallocation')
  end subroutine

  subroutine vfcp(tape,sys)
  ! Evaluate coulomb fitting coefficients
  ! Roberto Flores-Moreno, 2018
  implicit none
    integer :: tape
    type(nsystem) :: sys

    integer :: is,js
    logical :: found

    integer :: allocs,dabas,dbas
    real(8),allocatable :: p(:,:),g(:,:),j(:),x(:)

    do is=1,sys%ns
      dabas = basis_nfun(sys%auxis(is))
      allocate(j(dabas),x(dabas),stat=allocs)
      if (allocs.gt.0) call file_error('vfcp: Allocation failure')

      found = .false.
      do js=1,is-1
        if (sys%species_name(js).eq.sys%species_name(is)) then
          if (sys%basis(js)%filename.eq.sys%basis(is)%filename) then
            if (sys%nmoo(js).eq.sys%nmoo(is)) then
              call file_vector_aux(x,dabas,js,0,'x','read')
              found = .true.
            end if
          end if
        end if
      end do

      if (.not.found) then
        dbas = basis_nfun(sys%basis(is))
        allocate(p(dbas,dbas),g(dabas,dabas),stat=allocs)
        if (allocs.gt.0) call file_error('vfcp: Allocation failure')

        call file_matrix_bas(p,dbas,is,'density','read')
        call integrals_set_ppi_type(sys%interaction(is,is)%simple)
        call aux_j_vector(p,j,sys%mol,sys%basis(is),sys%auxis(is),dbas,dabas)
        call file_matrix_aux(g,dabas**2,is,is,'ig','read')
        x = matmul(g,j)

        deallocate(p,g,stat=allocs)
        if (allocs.gt.0) call file_error('vfcp: Deallocation failure')
      end if
      call file_vector_aux(x,dabas,is,0,'x','write')

      deallocate(j,x,stat=allocs)
      if (allocs.gt.0) call file_error('vfcp: Deallocation failure')
    end do
  end subroutine

  subroutine  build_fock_matrix_dir3(tape,sys,is,f,dbas,dft)
  ! Build Fock matrix using VFCP for J and K
  ! D. Mejia-Rodriguez et al. J. Chem. Phys. 141, 124114 (2014)
  ! Roberto Flores-Moreno, 2018
  implicit none
    type(nsystem) :: sys
    integer :: dbas,is,tape
    logical :: dft
    real(8) :: f(dbas,dbas)

    integer :: ibatch,imo,iset,jmo,js,ks,ll,n,nbatch,nmob,ul
    logical :: gocoul,goex
    logical :: sused(nmaxns)
    real(8) :: wex

    integer :: allocs,dabas,dabas2,dbas2,dmo,dset
    logical,allocatable :: ldfaux(:),ldfbas(:),auxptr(:,:),basptr(:,:)
    real(8),allocatable :: g(:,:),j(:),x(:)
    real(8),allocatable :: p(:,:),w(:,:),eri(:,:,:)

    gocoul = .true.
    wex = sys%interaction(is,is)%exchange(1)
    if (wex.eq.0.0) then
      goex = .false.
    else
      goex = .true.
    end if

    ! Exchange interaction
    if (goex) then
      call file_error(' 3-ERIs in exchange ')
    end if

    ! Classical interaction
    if (gocoul) then
      sused(1:sys%ns) = .false.
      do js=1,sys%ns
        if (.not.sused(js)) then
          sused(js) = .true.
          dabas = basis_nfun(sys%auxis(js))
          allocate(x(dabas),j(dabas),stat=allocs)
          if (allocs.gt.0) call file_error('vfcp_fock: Allocation failure')
          call file_vector_aux(x,dabas,js,0,'x','read')
          if (dft.and.(is.eq.js)) then
            call file_vector_aux(j,dabas,is,0,'z','read')
            x(1:dabas) = x(1:dabas) + j(1:dabas)
          end if
          do ks=js+1,sys%ns
            if (sys%species_name(js).eq.sys%species_name(ks)) then
              if (sys%basis(js)%filename.eq.sys%basis(ks)%filename) then
                call file_vector_aux(j,dabas,ks,0,'x','read')
                x(1:dabas) = x(1:dabas) +  j(1:dabas)
                if (dft.and.(is.eq.ks)) then
                  call file_vector_aux(j,dabas,is,0,'z','read')
                  x(1:dabas) = x(1:dabas) + j(1:dabas)
                end if
                sused(ks) = .true.
              end if
            end if
          end do
          call integrals_set_ppi_type(sys%interaction(is,js)%simple)
          call aux_build_fock_matrix(tape,sys%mol,sys%basis(is),sys%auxis(js),&
                                     f,x,dbas,dabas)
          deallocate(x,j,stat=allocs)
          if (allocs.gt.0) call file_error('vfcp_fock: Deallocation failure')
        end if
      end do
    end if
  end subroutine

end module
