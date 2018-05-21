#
# makefile
#

PROGRAM=smooth
VERSION=1

# Files
OBJ=smooth.o decomment.o gauss_elimination.o create_matrix.o free_matrix.o
SRC=smooth.c decomment.c gauss_elimination.c create_matrix.c free_matrix.c
HEAD=decomment.h revision.h matrix.h gauss_elimination.h create_matrix.h free_matrix.h

# C compiler
CC = gcc

# Compiler optimalisation (-O0 for debug, -O2 for ussual opt.)
OPT = -O2 #-O0 

# For gdb, OPT must be -O0  (-g for gdb, -pd for gproof)
DBG = #-pg

# Ostatni parametry prekladace (-Wall -Wextra -pedantic)
CFLAGS = -Wall -Wextra #-pedantic

# Linkovane knihovny  libefence.a = -lefence
LIB = #-lefence

# Cilum build, install, uninstall, clean a dist neodpovida primo zadny soubor
# (predstirany '.PHONY' target)

.PHONY: build
.PHONY: install
.PHONY: uninstall
.PHONY: clean
.PHONY: dist

# List of valid suffixes through the use of the .SUFFIXES special target.
.SUFFIXES: .c .o

# Prvni cil je implicitni, neni treba volat 'make build', staci 'make'.
# Cil build nema zadnou akci, jen zavislost.

build: $(PROGRAM)

# install závisi na prelozeni projektu, volat ho muze jen root
install: build
	cp $(PROGRAM) /usr/local/bin

# uninstall ma jenom akci a zadnou zavislost, volat ho muze jen root
uninstall:
	rm -f /usr/local/bin/$(PROGRAM)

# Clean files
clean:
	rm -f *.o $(PROGRAM)

# Source package
dist:
	tar czf $(PROGRAM)-$(VERSION).tgz $(SRC) $(HEAD) makefile readme.txt test.dat

# Slinkovani
$(PROGRAM): $(OBJ) makefile
	$(CC) $(OBJ) $(DBG) $(LIB) -o $(PROGRAM)

#
# The target form '.s1.s2', where .s1 and .s2 are currently valid suffixes, then
# it defines a transformation from *.s1 to *.s2 (double suffix interference)
# If a target has the form '.s1', where .s1 is a currently valid suffix,
# then it defines a transformation from *.s1 to * (single suffix interference)
.c.o: $(HEAD) makefile
	$(CC) $(CFLAGS) $(OPT) $(DBG) -c $<
