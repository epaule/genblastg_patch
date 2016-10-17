#
# G++ makefile for gBlast
# enable -DDEBUG to compile with edge/node debugging output
#
# To build, it should simply be a matter of
#
#             make all
#
#

SRCDIRS := . #current directory
SRCEXTS := .c .cpp # C/C++ program

CXX	= g++
CC	= gcc
#####################################################
# uncomment this for regular build
CXXFLAGS  	:= -g -Wfatal-errors -O2 #-Wall
CFLAGS 		:= -g -Wfatal-errors -O2 #-Wall
#CCFLAGS	= -g -Wfatal-errors
#####################################################
# uncomment this for debug build
#CXXFLAGS  	:= -g -Wfatal-errors -DDEBUG#-Wall
#CFLAGS 		:= -g -Wfatal-errors -DDEBUG#-Wall
#CCFLAGS	= -g -Wfatal-errors -DDEBUG#-Wall
#####################################################
LD	:= g++
LDFLAGS	:= 

# use this command to erase files.
RM = /bin/rm -f

# program executable file name.
#####################################################
# uncomment this for regular build
PROGRAM  := genblast
#####################################################
# uncomment this for debug build
#PROGRAM  := genblastdebug
#####################################################

# list of generated object files.
#OBJS = gBlast.o graph.o data_manager.o edge.o scores.o build.o contin.o prune.o stats.o trees.o info.o

SOURCES = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))
OBJS = $(foreach x,$(SRCEXTS), $(patsubst %$(x),%.o,$(filter %$(x),$(SOURCES))))
DEPS = $(patsubst %.o,%.d,$(OBJS))

# Rules for creating the dependency files (.d).
#---------------------------------------------------
%.d : %.c
	@$(CC) -MM -MD $(CFLAGS) $<

%.d : %.cpp
	@$(CC) -MM -MD $(CXXFLAGS) $<

# Rules for producing the objects.
#---------------------------------------------------
objs : $(OBJS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $<

# Rules for producing the executable.
#----------------------------------------------
$(PROGRAM) : $(OBJS)
ifeq ($(strip $(SRCEXTS)), .c) # C file
	$(CC) -o $(PROGRAM) $(OBJS) $(LDFLAGS)
else # C++ file
	$(CXX) -o $(PROGRAM) $(OBJS) $(LDFLAGS)
endif

-include $(DEPS)

# top-level rule, to compile everything.
all: $(PROGRAM)

# rule for cleaning re-compilable files.
clean :
	@$(RM) $(OBJS) *.d $(PROGRAM)

