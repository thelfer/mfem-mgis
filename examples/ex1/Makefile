DEBUG?=0

CONFIG_MK = $(MFEM_DIR)/../../../share/mfem/config.mk
MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

CPPFLAGS= $(MFEM_CPPFLAGS) -Wall -W -pedantic
MFRONT_FLAGS=

ifeq ($(DEBUG),1)
	CPPFLAGS += -g
        MFRONT_FLAGS += -g
endif

# Change "-I" into "-isystem" to avoid too much warnings.
# Avoid this change for /usr/include directory because it fails.
#$(info $$MFEM_INCFLAGS is [${MFEM_INCFLAGS}])
MFEM_ISYSTEM := $(shell echo $(MFEM_INCFLAGS) | sed "s%-I%-isystem %g;s%-isystem\ \ */usr/include%-I/usr/include%g") 

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
	mfront --obuild --interface=generic $(MFRONT_FLAGS) Plasticity.mfront

clean:
	rm -fr include/ src/ UniaxialTensileTest
