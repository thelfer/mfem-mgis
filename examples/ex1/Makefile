CPPFLAGS= -g -Wall -W -pedantic                             \
          -I $(MFEMMGIS_DIR)/../../../include      \
          -I $(MFEM_DIR)/../../../include            \
          -I $(MFrontGenericInterface_DIR)/../../../include
          

CXXFLAGS= -O2 -std=c++17

LDFLAGS= -L$(MFEMMGIS_DIR)/../../../lib  \
	 -L$(MFEM_DIR)/../../../lib  \
	 -L$(MFrontGenericInterface_DIR)/../../../lib  \
	 -lMFEMMGIS -lMFrontGenericInterface -lmfem 

ifeq ($(strip $(MFEM_LIB_MAKE)),)
else
	LDFLAGS+=$(MFEM_LIB_MAKE)
endif	

ifeq ($(strip $(MFEM_INC_MAKE)),)
else
	CPPFLAGS+=$(MFEM_INC_MAKE)
endif	


all: UniaxialTensileTest src/libBehaviour.so 

UniaxialTensileTest: UniaxialTensileTest.cxx
	$(MFEMMGIS_CXX) $(CPPFLAGS) $(CXXFLAGS) UniaxialTensileTest.cxx -o UniaxialTensileTest $(LDFLAGS)

src/libBehaviour.so: Plasticity.mfront
	mfront --obuild --interface=generic Plasticity.mfront

clean:
	rm -fr include/ src/ UniaxialTensileTest
