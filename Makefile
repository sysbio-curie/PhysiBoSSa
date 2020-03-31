VERSION := $(shell grep . VERSION.txt | cut -f1 -d:)
PROGRAM_NAME := ECM_try

CC := g++
# CC := g++-mp-7 # typical macports compiler name
# CC := g++-7 # typical homebrew compiler name 

# Check for environment definitions of compiler 
# e.g., on CC = g++-7 on OSX
ifdef PHYSICELL_CPP 
	CC := $(PHYSICELL_CPP)
endif

### MaBoSS configuration 
# MaBoSS max nodes
ifndef MABOSS_MAX_NODES
MABOSS_MAX_NODES = 64
endif

# MaBoSS directory
MABOSS_DIR = addons/PhysiBoSSa/MaBoSS-env-2.0/engine
CUR_DIR = $(shell pwd)
LIB := -L$(CUR_DIR)/$(MABOSS_DIR)/lib -lMaBoSS -ldl
INC := -I$(CUR_DIR)/$(MABOSS_DIR)/include -DMAXNODES=$(MABOSS_MAX_NODES)

# If max nodes > 64, change lib path 
ifeq ($(shell expr $(MABOSS_MAX_NODES) '>' 64), 1)
LIB := -L$(CUR_DIR)/$(MABOSS_DIR)/lib -lMaBoSS_$(MABOSS_MAX_NODES)n -ldl
endif

ARCH := native # best auto-tuning
# ARCH := core2 # a reasonably safe default for most CPUs since 2007
# ARCH := corei7
# ARCH := corei7-avx # earlier i7 
# ARCH := core-avx-i # i7 ivy bridge or newer 
# ARCH := core-avx2 # i7 with Haswell or newer
# ARCH := nehalem
# ARCH := westmere
# ARCH := sandybridge # circa 2011
# ARCH := ivybridge   # circa 2012
# ARCH := haswell     # circa 2013
# ARCH := broadwell   # circa 2014
# ARCH := skylake     # circa 2015
# ARCH := bonnell
# ARCH := silvermont
# ARCH := skylake-avx512
# ARCH := nocona #64-bit pentium 4 or later 

# CFLAGS := -march=$(ARCH) -Ofast -s -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11
CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11

COMPILE_COMMAND := $(CC) $(CFLAGS) 

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o \
BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o 

PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o \
PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o

PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o \
PhysiCell_pugixml.o PhysiCell_settings.o

# put your custom objects here (they should be in the custom_modules directory)

PhysiBoSS_module_OBJECTS := readXML.o vector3d.o writeSVG.o 

tinyxml_OBJECTS := tinyxml2.o

PhysiBoSSa_OBJECTS := utils.o maboss_network.o boolean_network.o

PhysiCell_custom_module_OBJECTS := custom.o custom_main.o

pugixml_OBJECTS := pugixml.o

PhysiCell_OBJECTS := $(BioFVM_OBJECTS)  $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS)
ALL_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECTS) $(PhysiBoSSa_OBJECTS) $(PhysiBoSS_module_OBJECTS) $(tinyxml_OBJECTS)

# compile the project 

all: MaBoSS main.cpp $(ALL_OBJECTS)
	$(COMPILE_COMMAND) $(INC) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp $(LIB)

# PhysiCell core components	

PhysiCell_phenotype.o: ./core/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_phenotype.cpp
	
PhysiCell_digital_cell_line.o: ./core/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_digital_cell_line.cpp

PhysiCell_cell.o: ./core/PhysiCell_cell.cpp
	$(COMPILE_COMMAND) $(INC) -c ./core/PhysiCell_cell.cpp 

PhysiCell_cell_container.o: ./core/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) $(INC) -c ./core/PhysiCell_cell_container.cpp 
	
PhysiCell_standard_models.o: ./core/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) $(INC) -c ./core/PhysiCell_standard_models.cpp 
	
PhysiCell_utilities.o: ./core/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) $(INC) -c ./core/PhysiCell_utilities.cpp 
	
PhysiCell_custom.o: ./core/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_custom.cpp 

PhysiCell_constants.o: ./core/PhysiCell_constants.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_constants.cpp
	
# BioFVM core components (needed by PhysiCell)
	
BioFVM_vector.o: ./BioFVM/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_vector.cpp 

BioFVM_agent_container.o: ./BioFVM/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_agent_container.cpp 
	
BioFVM_mesh.o: ./BioFVM/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_mesh.cpp 

BioFVM_microenvironment.o: ./BioFVM/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) $(INC) -c ./BioFVM/BioFVM_microenvironment.cpp 

BioFVM_solvers.o: ./BioFVM/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_solvers.cpp 

BioFVM_utilities.o: ./BioFVM/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_utilities.cpp 
	
BioFVM_basic_agent.o: ./BioFVM/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_basic_agent.cpp 
	
BioFVM_matlab.o: ./BioFVM/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: ./BioFVM/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_MultiCellDS.cpp
	
pugixml.o: ./BioFVM/pugixml.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/pugixml.cpp
	
# standard PhysiCell modules

PhysiCell_SVG.o: ./modules/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_SVG.cpp

PhysiCell_pathology.o: ./modules/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) $(INC) -c ./modules/PhysiCell_pathology.cpp

PhysiCell_MultiCellDS.o: ./modules/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) $(INC) -c ./modules/PhysiCell_MultiCellDS.cpp

PhysiCell_various_outputs.o: ./modules/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) $(INC) -c ./modules/PhysiCell_various_outputs.cpp

	
PhysiCell_pugixml.o: ./modules/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pugixml.cpp
	
PhysiCell_settings.o: ./modules/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_settings.cpp	
	
# user-defined PhysiCell modules

MaBoSS:
	cd $(MABOSS_DIR)/src && make MAXNODES=$(MABOSS_MAX_NODES) install_alib

maboss_network.o: ./addons/PhysiBoSSa/src/maboss_network.cpp
	$(COMPILE_COMMAND) $(INC) -c ./addons/PhysiBoSSa/src/maboss_network.cpp

utils.o: ./addons/PhysiBoSSa/src/utils.cpp
	$(COMPILE_COMMAND) -c ./addons/PhysiBoSSa/src/utils.cpp

boolean_network.o: ./addons/PhysiBoSSa/src/boolean_network.cpp
	$(COMPILE_COMMAND) $(INC) -c ./addons/PhysiBoSSa/src/boolean_network.cpp

custom.o: ./custom_modules/custom.cpp 
	$(COMPILE_COMMAND) $(INC) -c ./custom_modules/custom.cpp

custom_main.o: ./custom_modules/custom_main.cpp 
	$(COMPILE_COMMAND) $(INC) -c ./custom_modules/custom_main.cpp


# PhysiBoSS modules from Gaelle's code:

readXML.o: ./core/base/readXML.cpp
	$(COMPILE_COMMAND) -c ./core/base/readXML.cpp

vector3d.o: ./core/base/vector3d.cpp
	$(COMPILE_COMMAND) -c ./core/base/vector3d.cpp

writeSVG.o: ./core/base/writeSVG.cpp
	$(COMPILE_COMMAND) -c ./core/base/writeSVG.cpp

tinyxml2.o: ./core/tinyxml/tinyxml2.cpp
	$(COMPILE_COMMAND) -c ./core/tinyxml/tinyxml2.cpp


# cleanup

reset:
	rm -f *.cpp 
	cp ./sample_projects/Makefile-default Makefile 
	rm -f ./custom_modules/*
	touch ./custom_modules/empty.txt 
	touch ALL_CITATIONS.txt 
	rm ALL_CITATIONS.txt 
	cp ./config/PhysiCell_settings-backup.xml ./config/PhysiCell_settings.xml 
	rm -rf ./config/boolean_network/
	rm -rf ./scripts

MaBoSS-clean:
	cd $(MABOSS_DIR)/src && make clean
	
clean: MaBoSS-clean
	rm -f *.o
	rm -f $(PROGRAM_NAME)*
	
data-cleanup:
	rm -f *.mat
	rm -f *.xml
	rm -f *.svg
	rm -f ./output/*
	touch ./output/empty.txt
	
# archival 

checkpoint: 
	zip -r $$(date +%b_%d_%Y_%H%M).zip Makefile *.cpp *.h config/*.xml custom_modules/* 

zip:
	zip -r latest.zip Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/* 
	cp latest.zip $$(date +%b_%d_%Y_%H%M).zip
	cp latest.zip VERSION_$(VERSION).zip 
	mv *.zip archives/
	
tar:
	tar --ignore-failed-read -czf latest.tar Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/* 
	cp latest.tar $$(date +%b_%d_%Y_%H%M).tar
	cp latest.tar VERSION_$(VERSION).tar
	mv *.tar archives/

unzip: 
	cp ./archives/latest.zip . 
	unzip latest.zip 
	
untar: 
	cp ./archives/latest.tar .
	tar -xzf latest.tar
