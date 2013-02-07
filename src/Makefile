CXXFLAGS = -O4 -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS
CXXFLAGSDEBUG = -g -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS

# C++ flags for different platform
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
   CXXFLAGS += 
endif

SOURCE := $(wildcard *.cpp)
OBJ = $(SOURCE:.cpp=.o)
TARGET = polymutt
LIB = ../core/libcore.a ../base/lib-base.a ../libVcf/lib-vcf.a
TABIX_LIB = ../third/tabix/libtabix.a
PCRE_LIB = -lpcre

all: $(TARGET)
$(TARGET): $(OBJ)
	g++ $(CXXFLAGS) -o $(TARGET) $(OBJ) $(LIB) $(TABIX_LIB) $(PCRE_LIB) -lm -lz -lgomp -lbz2
%.o: %.cpp
	g++ -c $(CXXFLAGS) -fopenmp -I../core -I../base -I../libVcf -I../third/tabix -I../third/pcre/include $<

clean:
	rm -f *.o $(TARGET)

