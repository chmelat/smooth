#
# Makefile for smooth project
#

# Program name and version
PROGRAM = smooth
VERSION=$(shell grep '^#define VERSION' revision.h | sed 's/.*"\(.*\)".*/\1/')

# Installation directories
PREFIX = /usr/local
# Use this for user installation: $(HOME)
USER_PREFIX = $(HOME)

# Files
SRC = smooth.c decomment.c tikhonov.c polyfit.c savgol.c butterworth.c grid_analysis.c
OBJ = $(SRC:.c=.o)
HEAD = decomment.h revision.h tikhonov.h polyfit.h savgol.h butterworth.h grid_analysis.h

# C compiler and flags
CC = clang
# Standard optimization level (production)
OPT = -O2
# Warning flags
WFLAGS = -Wall -Wextra -pedantic

# Debug flags (uncomment for debugging)
# DBG = -g -O0

# Memory checking tools
# For memory leak detection (uncomment to use)
# MEMCHECK = -lefence
# Alternative: valgrind (used via 'make memcheck')

# Library settings
# Default library paths (can be overridden with environment variables)
LIBDIR ?= $(HOME)/lib
INCDIR ?= $(HOME)/include
BINDIR ?= $(HOME)/bin

# System-wide installation paths
SYS_LIBDIR = $(PREFIX)/lib
SYS_INCDIR = $(PREFIX)/include
SYS_BINDIR = $(PREFIX)/bin

# Libraries
LIBINCLUDE = -I$(INCDIR)
LIBPATH = -L$(LIBDIR)
LIB = -llapack -lblas $(MEMCHECK) -lm

# PHONY targets (not corresponding to files)
.PHONY: all build debug memcheck install install-user uninstall uninstall-user clean dist check help

# Default target
all: build

# Build target (production build)
build: CFLAGS = $(WFLAGS) $(OPT)
build: $(PROGRAM)

# Debug build
debug: CFLAGS = $(WFLAGS) -g -O0
debug: clean $(PROGRAM)
	@echo "Built with debug information"

# Memory check build (for use with valgrind)
memcheck: debug
	@echo "Run with: valgrind --leak-check=full ./$(PROGRAM) [options]"

# Install to system directories (requires root access)
install: build
	@echo "Installing to $(SYS_BINDIR)..."
	install -d $(SYS_BINDIR)
	install -m 755 $(PROGRAM) $(SYS_BINDIR)
	@echo "Installation complete"

# Install to user's home directory
install-user: build
	@echo "Installing to $(BINDIR)..."
	install -d $(BINDIR)
	install -m 755 $(PROGRAM) $(BINDIR)
	@echo "User installation complete"

# Uninstall from system
uninstall:
	@echo "Uninstalling from $(SYS_BINDIR)..."
	rm -f $(SYS_BINDIR)/$(PROGRAM)
	@echo "Uninstallation complete"

# Uninstall from user directory
uninstall-user:
	@echo "Uninstalling from $(BINDIR)..."
	rm -f $(BINDIR)/$(PROGRAM)
	@echo "User uninstallation complete"

# Clean generated files
clean:
	@echo "Cleaning up..."
	rm -f *.o $(PROGRAM)
	@echo "Clean complete"

# Create source package
dist:
	@echo "Creating distribution package..."
	mkdir -p $(PROGRAM)-$(VERSION)
	cp $(SRC) $(HEAD) Makefile README.md test_data.txt $(PROGRAM)-$(VERSION)/
	tar czf $(PROGRAM)-$(VERSION).tgz $(PROGRAM)-$(VERSION)
	rm -rf $(PROGRAM)-$(VERSION)
	@echo "Created $(PROGRAM)-$(VERSION).tgz"

# Help information
help:
	@echo "Makefile for $(PROGRAM) version $(VERSION)"
	@echo ""
	@echo "Available targets:"
	@echo "  all           : Default target, same as 'build'"
	@echo "  build         : Build the program with optimizations"
	@echo "  debug         : Build with debug information"
	@echo "  memcheck      : Build for memory checking with valgrind"
	@echo "  install       : Install program to $(SYS_BINDIR) (may require root)"
	@echo "  install-user  : Install program to $(BINDIR)"
	@echo "  uninstall     : Remove program from $(SYS_BINDIR)"
	@echo "  uninstall-user: Remove program from $(BINDIR)"
	@echo "  clean         : Remove object files and executable"
	@echo "  dist          : Create source distribution package"
	@echo "  help          : Display this help message"
	@echo ""
	@echo "Environment variables:"
	@echo "  LIBDIR        : Library directory (default: $(HOME)/lib)"
	@echo "  INCDIR        : Include directory (default: $(HOME)/include)"
	@echo "  BINDIR        : Binary installation directory (default: $(HOME)/bin)"
	@echo ""
	@echo "Example usage:"
	@echo "  make                  - Build the program"
	@echo "  make debug            - Build with debug information"
	@echo "  make install-user     - Install to user's bin directory"
	@echo "  LIBDIR=/opt/lib make  - Use custom library path"

# Link the program
$(PROGRAM): $(OBJ) Makefile
	@echo "Linking $(PROGRAM)..."
	$(CC) $(LIBPATH) $(OBJ) $(DBG) $(LIB) -o $(PROGRAM)
	@echo "Build complete"

# Compile source files
%.o: %.c $(HEAD) Makefile
	@echo "Compiling $<..."
	$(CC) $(LIBINCLUDE) $(CFLAGS) $(DBG) -c $<

# Explicit dependencies
smooth.o: smooth.c decomment.h revision.h Makefile
decomment.o: decomment.c decomment.h Makefile

#
# End
#

