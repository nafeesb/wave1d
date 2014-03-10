CXX = g++
#CXX = icpc
CXXFLAGS  = -I../fltk-1.1 
CXXFLAGS += -I/System/Library/Frameworks/OpenGL.framework/Headers

LDFLAGS = -L../fltk-1.1/lib -lfltk_gl -lfltk -lpthread 
#LDFLAGS += -L/opt/intel/cc/9.1.029/lib
LDFLAGS += -framework OpenGL -framework AGL -framework Carbon -framework ApplicationServices
LDFLAGS += -lm
# LDFLAGS += -limf -lirc
REZT = /Developer/Tools/Rez -t APPL

EXECS = $(MAKECMDGOALS)
$(EXECS)_SRCS = $(EXECS).C

$(EXECS): $($(EXECS)_SRCS:.C=.o)
	  $(CXX) -o $(@) $($(@)_SRCS:.C=.o) $(LDFLAGS)
	  $(REZT) -o $(@) ../fltk-1.1/FL/mac.r
