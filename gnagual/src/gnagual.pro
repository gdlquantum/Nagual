TEMPLATE = app
TARGET = 
DEPENDPATH += .
INCLUDEPATH += .
LIBS += -lGLU

# Input
QT += opengl 
HEADERS += \
           DwElement.h \
           DwGeo.h \
           DwImage.h \
           DwPanel.h \
           DwPlot.h \
           Element.h \
           Evaluator.h \
           GNagual.h \
           Grid.h \
           Math.h \
           Matrix.h \
           Model.h \
           Molecule.h \
           Parameter.h \
           Plotter.h \
           Set.h \
           Shell.h \
           Spin.h \
           Surface.h \
           System.h \
           Units.h \
           Vector.h \ 
           Viewer.h 
           
SOURCES += \
           Atom.cpp \
           DwElement.cpp \
           DwGeo.cpp \
           DwImage.cpp \
           DwPanel.cpp \
           DwPlot.cpp \
           Element.cpp \
           Evaluator.cpp \
           GNagual.cpp \
           Grid.cpp \
           main.cpp \
           Math.cpp \
           Matrix.cpp \
           Model.cpp \
           Molecule.cpp \
           Plotter.cpp \
           Set.cpp \
           Shell.cpp \
           Spin.cpp \
           Surface.cpp \
           System.cpp \
           Units.cpp \
           Vector.cpp \ 
           Viewer.cpp 
