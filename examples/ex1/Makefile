CONFIG_MK = $(MFEM_DIR)/../../../share/mfem/config.mk
MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

CPPFLAGS= $(MFEM_CPPFLAGS) -g -Wall -W -pedantic                    

MFEM_ISYSTEM := $(subst -I,-isystem ,$(MFEM_INCFLAGS))

INCFLAGS = $(MFEM_ISYSTEM) \
          -isystem $(MFEMMGIS_DIR)/../../../include      \
          -isystem $(MFrontGenericInterface_DIR)/../../../include

CXXFLAGS =  -O2 -std=c++17 # $(MFEM_CXXFLAGS)

LDFLAGS= -L$(MFEMMGIS_DIR)/../../../lib  \
	 -L$(MFrontGenericInterface_DIR)/../../../lib  \
	 -lMFEMMGIS -lMFrontGenericInterface $(MFEM_LIBS) 

all: UniaxialTensileTest src/libBehaviour.so 

UniaxialTensileTest: UniaxialTensileTest.cxx
	$(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) UniaxialTensileTest.cxx -o UniaxialTensileTest $(LDFLAGS)

src/libBehaviour.so: Plasticity.mfront
	mfront --obuild --interface=generic Plasticity.mfront

clean:
	rm -fr include/ src/ UniaxialTensileTest
