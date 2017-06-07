CXX := $(shell ${ROOTSYS}/bin/root-config --cxx)
MARCH := $(shell ${ROOTSYS}/bin/root-config --arch)
LD:=$(CXX)

SRC=./src
LIB=./lib
INC=./include
OBJ=./obj
EXE=./exe

VERSION      := $(shell $(ROOTSYS)/bin/root-config --version | cut -b1-4)
ifeq ($(VERSION),5.27)
VERSIONP=
else
VERSIONP = $(VERSION)
endif

ifndef  SLC6system
SLC6system=
ifneq "$(wildcard /etc/redhat-release)" ""
# SLC6system := $(shell cat /etc/redhat-release | grep "Scientific Linux CERN SLC release 6")
  SLC6system := $(shell cat /etc/redhat-release | grep "Scientific Linux release 6")
endif
endif
ifndef  SLC64system
SLC64system=
ifneq "$(wildcard /etc/redhat-release)" ""
 SLC64system := $(shell cat /etc/redhat-release | grep "Scientific Linux CERN SLC release 6.4")
endif
ifndef SLC64system
ifneq "$(wildcard /etc/redhat-release)" ""
 SLC64system := $(shell cat /etc/redhat-release | grep "Scientific Linux CERN SLC release 6.5")
endif
endif
endif


#CFLAGS:= -g -Wno-deprecated
CFLAGS:= -O3 -Wno-deprecated

ifdef PGTRACK
	AUXFLAGS += -D_PGTRACK_
endif

ifdef ONMACOSX
        AUXFLAGS += -D_ONMACOSX_
endif

SOFLAGS=
ifeq "$(system)" "Darwin"
      SOFLAGS = --dynamic --shared --dynamiclib -undefined dynamic_lookup
else
      SOFLAGS = -shared -fPIC
endif

CFLAGS += $(shell $(ROOTSYS)/bin/root-config --auxcflags) -fPIC
CFLAGS += $(AUXFLAGS)

ifdef ONMACOSX
	CFLAGSFORCINT=$(subst -stdlib=libc++ -std=c++11,,$(CFLAGS))
else
	CFLAGSFORCINT=$(CFLAGS)
endif

ifndef ONMACOSX
	AC_LIBS=-L$(QTDIR)/lib -lQtCore
else
	AC_LIBS = -F$(QTDIR) -framework QtCore
endif
ROOTVERSION := $(root-config --version)
ROOTDYNAMICLIBS := $(shell root-config --glibs)
ROOTDYNAMICLIBS += -lMinuit -lRooFitCore -lRooFit 

EXTRALIBS= -ldl -lTreePlayer -lProof -lProofPlayer
ifndef ONMACOSX
	EXTRALIBS += -lRFIO -lTMVA -lXMLIO -lMLP #-lNetx
else
	EXTRALIBS += -lTMVA -lXMLIO -lMLP
endif

EXTRALIBS=
# empty


INCLUDES=-I./ -I$(ROOTSYS)/include -I$(INC)

TARGETLIB := RooSpectra

SOURCES:= RooAsymmLorentzian.C RooBrokenPowerLaw.C RooBrokenPowerLaw2.C RooPowerLaw.C SlidingWindowFit.C
#IndexFitter.C

OBJECTS += $(SOURCES:%.C=./obj/%.o)

default: shared-lib static-lib root6 #test

#-----------------------------------------------

shared-lib: $(LIB)/lib$(TARGETLIB).so

$(LIB)/lib$(TARGETLIB).so: $(OBJECTS) $(OBJ)/Dict.o
	@echo "Creating shared library \"$@\" ..."
	@if ! [ -d $(LIB) ] ; then mkdir -p $(LIB); fi
	$(CXX) $(SOFLAGS) -o $(LIB)/lib$(TARGETLIB).so $(OBJECTS) $(OBJ)/Dict.o $(ROOTDYNAMICLIBS) $(EXTRALIBS) 

root6:	$(OBJ)/Dict_rdict.pcm
	@echo copying pcm file to $(LIB)
	cp -v $(OBJ)/Dict_rdict.pcm $(LIB)/Dict_rdict.pcm

static-lib: $(LIB)/lib$(TARGETLIB)a.a

$(LIB)/lib$(TARGETLIB)a.a: $(OBJECTS) $(OBJ)/Dict.o
	@echo "Creating static library \"$@\" ..."
	@if ! [ -d $(LIB) ] ; then mkdir -p $(LIB); fi
	ar rv $(LIB)/lib$(TARGETLIB)a.a $(OBJECTS) $(OBJ)/Dict.o 

#-----------------------------------------------

test: $(OBJ)/test.o
	@if ! [ -d $(EXE) ] ; then mkdir -p $(EXE); fi
	$(CXX) $(SRC)/test.C -o $(EXE)/$@ $(CFLAGS) $(INCLUDES) $(ROOTDYNAMICLIBS) $(EXTRALIBS)
#----------------------------------------------

./obj/%.o: $(SRC)/%.cpp
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(CXX) $(CFLAGS) $(ROOTDYNAMICLIBS) $(INCLUDES) -c $< -o $@

./obj/%.o: $(SRC)/%.cxx
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(CXX) $(CFLAGS) $(ROOTDYNAMICLIBS) $(INCLUDES) -c $< -o $@

./obj/%.o: $(SRC)/%.c
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(CXX) $(CFLAGS) $(ROOTDYNAMICLIBS) $(INCLUDES) -c $< -o $@

./obj/%.o: $(SRC)/%.C
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(CXX) $(CFLAGS) $(ROOTDYNAMICLIBS) $(INCLUDES) -c $< -o $@

#----------------------------------------------

$(OBJ)/Dict.o: $(OBJ)/Dict.cxx
	@echo Compiling  $< ...
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(CXX) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJ)/Dict.cxx: $(INC)/Roo*.h $(INC)/Sliding*.h $(INC)/LinkDef.h
	@echo Creating  $@ ...
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(CFLAGSFORCINT) $(INCLUDES) $^
#----------------------------------------------

clean: 
	rm -f $(OBJ)/*
	rm -f $(LIB)/lib$(TARGETLIB).so
	rm -f $(LIB)/lib$(TARGETLIB)a.a
	rm -f $(SRC)/*.so
	rm -f $(SRC)/*.d
	rm -f $(EXE)/*
	rm -f test

