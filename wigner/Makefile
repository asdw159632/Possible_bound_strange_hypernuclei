#
# setup control
#

#CXX           = g++
#CXXFLAGS      = -O -Wall -fPIC
#LD            = g++
#LDFLAGS       = -O
#SOFLAGS       = -shared


DEBUG = TRUE

TOP     := $(shell pwd)/
OBJ     := $(TOP)obj/
BIN     := $(TOP)bin/
LIB     := lib/
MYLIB   := ../mylib/*.o
MYINC		:= ../myinc
SRC     := $(TOP)src/
MAIN    := $(TOP)main/
INC     := $(TOP)src/
#STAF    := $(TOP)staf/
#STAFINC := $(STAF)inc/

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMathMore
ROOTGLIBS     = $(shell root-config --glibs)

#GFLIBS = -L/usr/local/gfortran/lib
#GFLIBS += -L/usr/local/gfortran/lib/i386
#GFLIBS += -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin14/4.9.1 -lgfortran

#GFLIBS = -L/Users/szhang/bin/gnuGCC/gcc-4.9.2-install/lib -lgfortran
#GFLIBS = -L/Users/szhang/bin/gnuGCC/gcc-5.3.0-install/lib -lgfortran

GFLIBS = -lgfortran 

ROOTINC := $(ROOTSYS)/include/

#
# set up dependencies
#
vpath %.h   $(ROOTINC):$(INC)
#vpath %.inc $(INC):$(STAFINC)
vpath %.cc  $(SRC):$(MAIN)
vpath %.cxx  $(SRC):$(MAIN)
vpath %.c   $(SRC):$(MAIN)
vpath %.C   $(SRC):$(MAIN)
vpath %.f   $(SRC)
vpath %.o   $(BIN)
#	
# set up compilers
#	

CPP     = g++
CPPFLAGS = -O -Wall -fPIC -I $(ROOTINC) -I$(INC) -I$(MYINC) -Wno-deprecated
CC      = gcc
ifeq ($(OPT),TRUE)	
#CPPFLAGS += -g 
endif	
ifeq ($(DEBUG), TRUE)
#CPPFLAGS += -g
endif
#
CFLAGS   = -w -I$(INC) $(MYINC)
#-I$(STAFINC) 
ifeq ($(OPT), TRUE)
#CFLAGS += -fast
endif
ifeq ($(DEBUG), TRUE)
#CFLAGS += -g	
endif	
#
FF      = gfortran 
#FFLAGS   = -c -e -xl -PIC -DSunOS -Dsun4os5 -I$(INC) -I$(STAFINC) 
FFLAGS   = -c -e -xl -DSunOS -Dsun4os5 -I$(INC)
#-I$(STAFINC) 
ifeq ($(OPT), TRUE)
FFLAGS += -fast
endif
ifeq ($(DEBUG), TRUE)
#FFLAGS += -g
endif
#
#CERNLIBS := $(shell $(CERN_ROOT)/bin/cernlib  pawlib graflib packlib mathlib kernlib)
#

ifeq ($(PURIFY), TRUE)
LD = purify $(CPP)
else
LD     = $(CPP)
endif
#LDFLAGS = -t -z muldefs -DSUN -DSOLARIS -Dsun 
ifeq ($(OPT), TRUE)
#LDFLAGS += -fast
endif
#
CSRCS := $(notdir $(wildcard $(SRC)*.c)) 
FSRCS := $(notdir $(wildcard $(SRC)*.f))
INCS  := $(notdir $(wildcard $(INC)*.h)) $(notdir $(wildcard $(INC)*.inc))
MAINS := $(notdir $(wildcard $(MAIN)*.c))
#OBJS  := $(FSRCS:.f=.o) $(CSRCS:.c=.o)
#
NOOPTOBJS :=  $(BIN)Raw.o $(BIN)TrkSort.o
#	
SOFLAGS = -G

SOURCE = $(wildcard $(SRC)*.cxx)
FSOURCE =$(wildcard $(SRC)*.f)
OBJS = $(patsubst $(SRC)%.cxx,$(OBJ)%.o,$(SOURCE)) $(patsubst $(SRC)%.f,$(OBJ)%.o,$(FSOURCE))

############## Make Executables ####################################
all: analysis
#all: ds kaka pika pipi lambdac charged d0rho thetapp l1520 thetapp2 d0


analysis : $(OBJS)
	$(LD) $(LDFLAGS) $(LDINC) $(XM_LIBS) $^ $(LDLIB) $(ROOTLIBS) $(ROOTLOADLIBS) $(MYLIB) $(GFLIBS) -o $(BIN)$(notdir $@)
	@echo


######################################################
$(OBJ)%.o : 	$(MAIN)%.c $(INCS)
	$(CC)     $(CFLAGS) -c $(MAIN)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(MAIN)%.cc $(INCS)
	$(CPP)     $(CPPFLAGS) -c $(MAIN)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(SRC)%.cxx $(INCS)
	$(CPP)     $(CPPFLAGS) -c $(SRC)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(MAIN)%.C $(INCS)
	$(CPP)     $(CPPFLAGS) -c $(MAIN)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(SRC)%.c $(INCS)
	$(CC)     $(CFLAGS) -c $(SRC)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(SRC)%.C $(INCS)
	$(CPP)     $(CPPFLAGS) -c $(SRC)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(SRC)%.f $(INCS)
	$(FF)     $(FFLAGS) -c $(SRC)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
foo:
	@echo $(CSRCS)
	@echo
	@echo $(FSRCS)
	@echo
	@echo $(MAINS)
	@echo
	@echo $(OBJS)
	@echo
	@echo $(INCS)
	@echo
clean:
	rm -f $(OBJ)*.o 











