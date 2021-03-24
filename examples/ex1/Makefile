CPPFLAGS= -g -Wall -W -pedantic                             \
          -isystem $(MFEMMGIS_DIR)/../../../include      \
          -isystem $(MFEM_DIR)/../../../include            \
          -isystem $(HYPRE_DIR)/include            \
          -isystem $(METIS_DIR)/include            \
          -isystem $(MFrontGenericInterface_DIR)/../../../include
          

CXXFLAGS= -O2 -std=c++17

LDFLAGS= -L$(MFEMMGIS_DIR)/../../../lib  \
	 -L$(MFEM_DIR)/../../../lib  \
	 -L$(MFrontGenericInterface_DIR)/../../../lib  \
	 -L$(METIS_DIR)/lib  \
	 -L$(HYPRE_DIR)/lib  \
	 -lMFEMMGIS -lMFrontGenericInterface -lmfem -lHYPRE -lmetis 

all: UniaxialTensileTest src/libBehaviour.so 

UniaxialTensileTest: UniaxialTensileTest.cxx
	$(MFEMMGIS_CXX) $(CPPFLAGS) $(CXXFLAGS) UniaxialTensileTest.cxx -o UniaxialTensileTest $(LDFLAGS)

src/libBehaviour.so: Plasticity.mfront
	mfront --obuild --interface=generic Plasticity.mfront

clean:
	rm -fr include/ src/ UniaxialTensileTest
