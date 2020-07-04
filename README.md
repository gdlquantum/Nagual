<p align="center">
  <img src="logoNGL.jpg" />
</p>

## License

 1. If you have got an integer number version
   you are bound to GNU/GPL version 3 license,
   a copy of which should be included in the
   package, in the file named GNU_GPLv3.txt .

 2. If you have got a fractional number version then
   you need to ask the authors (r.flores@academicos.udg.mx)
   for authorization
   on any use you are planning for Nagual/GNagual software.

## Installation

### Compiling Nagual 1
 
 Put yourself in the “nagual/src” directory and once there type the following:
 ```
 >>> ./configure
 >>> make
 ```
 If everything goes fine the executable, named Nagual.1.exe, will be created and
 located on the “bin” directory. If NOT, look inside the configure file for 
 adjustments of configuration. In particular check that libraries like LAPACK have 
 the corresponding directory properly addressed.

### Compiling GNagual 1

 Put yourself in the “gnagual/src” directory and once there type the following:
 ```
 >>> qmake
 >>> make
 ```
 If everything goes fine the executable, named GNagual.1.exe, will be created and
 located on the “bin” directory. If NOT, look inside the gnagual.pro file for 
 adjustments of QT configuration. Additionally you will find useful to look into QT
 documentation.

## Running calculations

 For a direct calculation

  1) Prepare an input named "Nagual.inp". Read the manual for 
     preparation instructions
  2) Locate your sesion at the same directory where "Nagual.inp" is located
  3) Execute nagual.x.y.exe

## Any further questions? Please contact

 Roberto Flores-Moreno  
 Chemistry Department  
 University of Guadalajara  
 Blvd. Marcelino Garcia Barragan 1421  
 Col. Olimpica, Guadalajara Jalisco  
 C.P. 44430, Mexico  

 Phone: +52 33 14957095  
 E-mail: r.flores@academicos.udg.mx  
 Twiter: @gdlquantum
