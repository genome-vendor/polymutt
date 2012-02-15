CXXFLAGS = -O4 -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS

# C++ flags for different platform
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
   CXXFLAGS += 
endif

SOURCE := $(wildcard *.cpp)
OBJ = $(SOURCE:.cpp=.o)
TARGET = polymutt
LIB = ../core/libcore.a

all: $(TARGET)
$(TARGET): $(OBJ)
	g++ $(CXXFLAGS) -o $(TARGET) $(OBJ) $(LIB) -lm -lz -lgomp
%.o: %.cpp
	g++ -c $(CXXFLAGS) -fopenmp -I../core $<
clean:
	rm -f *.o $(TARGET)