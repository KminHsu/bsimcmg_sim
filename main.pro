######################################################################
# Automatically generated by qmake (2.01a) Sat Apr 18 07:06:14 2015
######################################################################

TEMPLATE = app
TARGET = bsimcmg_module
DEPENDPATH += .
INCLUDEPATH += . \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./ParallelDistPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./DeviceModelPKG/NeuronModels/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./DeviceModelPKG/EXTSC/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./DeviceModelPKG/TCADModels/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./DeviceModelPKG/ADMS/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./DeviceModelPKG/OpenModels/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./DeviceModelPKG/Core/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./TopoManagerPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./UtilityPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./DakotaLinkPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./NonlinearSolverPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./LinearAlgebraServicesPKG/ksparse/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./LinearAlgebraServicesPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./AnalysisPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./IOInterfacePKG/Output/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./IOInterfacePKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./ErrorHandlingPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./CircuitPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./TimeIntegrationPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./MultiTimePDEPKG/include \
               /home/kminhsu/Desktop/reg/Xyce-6.2-build/src/./LoaderServicesPKG/include \
               /usr/local/XyceLibs/Serial/include

QMAKE_CXXFLAGS += -g 
QMAKE_CXX = /usr/bin/g++
QMAKE_LINK = /usr/bin/g++

CONFIG += core

# Input
#HEADERS += 
SOURCES += main.cpp

LIBS += -L/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7 -L/usr/local/XyceLibs/Serial/lib -L/usr/local/lib -lxyce -lloca -lnox -lbelosepetra -lbelos -laztecoo -lifpack -ltrilinoscouplings -lamesos -lumfpack -lcolamd -lepetraext -lepetra -ltriutils -lteuchosparameterlist -lteuchoscomm -lteuchosnumerics -lteuchosremainder -lteuchoscore -lamd -llapack -lblas -L/usr/local/lib -L/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7 -L/usr/local/XyceLibs/Serial/lib -L/lib -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lsuitesparseconfig -lcholmod -ldl -lrt -lutil -lgfortranbegin -lgfortran -lm -lkokkoscore -lfftw3


