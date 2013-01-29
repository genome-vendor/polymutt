LIB = lib-base.a
LIB_DBG = lib-dbg-base.a
BASE = Argument Exception IO OrderedMap Regex TypeConversion Utils Logger RangeList SimpleMatrix
OBJ = $(BASE:=.o)

DEFAULT_CXXFLAGS = -I. -I../third/tabix -I../third/pcre/include

.PHONY: all release debug
all: release

# 'make release' will: clean current build; build with -O4 flags; backup executable files to release/ folder
# NOTE: don't use -j flag!
release: CXXFLAGS = -O4 $(DEFAULT_CXXFLAGS)
release: $(LIB)

debug: CXXFLAGS = -Wall -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(LIB_DBG)

$(LIB): $(OBJ)
	ar rcs $(LIB) $(OBJ)
$(LIB_DBG): $(OBJ)
	ar rcs $(LIB_DBG) $(OBJ)

%: %.cpp # remove implicit Makefile rule of converting ABC.cpp to ABC
%.o: %.cpp %.h
	g++ -c $(CXXFLAGS) $<

clean:
	rm -rf $(LIB) $(LIB_DBG) *.o
