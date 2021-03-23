CPPFLAGS= -g -Wall -W -pedantic -DMFEM_USE_MPI=ON                              \
	  -isystem /home/hc265945/spack_codes/mfem-mgis/install/include      \
	  -isystem /home/hc265945/spack_codes/mfem/install/include            \
          -isystem /usr/include/suitesparse                                    \
	  -isystem $(MGIS_DIR)/include
	  

CXXFLAGS= -O2 -std=c++17

LDFLAGS= -L/home/hc265945/spack_codes/mfem-mgis/install/lib  \
	 -L/home/hc265945/spack_codes/mfem/install/lib  \
	 -lMFEMMGIS -lMFrontGenericInterface -lmfem -lHYPRE -lmetis -L${HYPRE_DIR}/lib -L${METIS_DIR}/lib

all: UniaxialTensileTest MPIUniaxialTensileTest src/libBehaviour.so 

UniaxialTensileTest: UniaxialTensileTest.cxx
	mpicxx $(CPPFLAGS) $(CXXFLAGS) UniaxialTensileTest.cxx -o UniaxialTensileTest $(LDFLAGS)

MPIUniaxialTensileTest: MPIUniaxialTensileTest.cxx
	mpicxx $(CPPFLAGS) $(CXXFLAGS) MPIUniaxialTensileTest.cxx -o MPIUniaxialTensileTest $(LDFLAGS)

src/libBehaviour.so: Plasticity.mfront
	mfront --obuild --interface=generic Plasticity.mfront

clean:
	rm -fr include/ src/ UniaxialTensileTest
