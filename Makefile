#
# Makefile for smooth project
#

# Program name and version
PROGRAM = smooth
VERSION=$(shell grep '^#define VERSION' revision.h | sed 's/.*"\(.*\)".*/\1/')

# Installation directories
PREFIX = /usr/local

# Files
SRC = smooth.c decomment.c tikhonov.c polyfit.c savgol.c butterworth.c grid_analysis.c timestamp.c
OBJ = $(SRC:.c=.o)
HEAD = decomment.h revision.h tikhonov.h polyfit.h savgol.h butterworth.h grid_analysis.h timestamp.h

# Test files
TEST_DIR = tests
TEST_SRC = $(TEST_DIR)/test_main.c $(TEST_DIR)/test_grid_analysis.c $(TEST_DIR)/test_polyfit.c $(TEST_DIR)/test_savgol.c $(TEST_DIR)/test_tikhonov.c $(TEST_DIR)/test_butterworth.c $(TEST_DIR)/test_timestamp.c $(TEST_DIR)/test_helpers.c $(TEST_DIR)/grid_helpers.c $(TEST_DIR)/unity.c
TEST_OBJ = $(TEST_SRC:.c=.o)
TEST_MODULES = grid_analysis.o polyfit.o savgol.o tikhonov.o butterworth.o timestamp.o # Modules being tested (without main program)
TEST_RUNNER = $(TEST_DIR)/test_runner

# C compiler and flags
CC = clang
# Standard optimization level (production)
OPT = -O2
# Warning flags
WFLAGS = -Wall -Wextra -pedantic

# Memory checking tools
# For memory leak detection (uncomment to use)
# MEMCHECK = -lefence
# Alternative: valgrind (used via 'make memcheck')

# Library settings
# Default library paths (can be overridden with environment variables)
LIBDIR ?= $(HOME)/lib
INCDIR ?= $(HOME)/include
BINDIR ?= $(HOME)/bin

# System-wide installation path
SYS_BINDIR = $(PREFIX)/bin

# Libraries
LIBINCLUDE = -I$(INCDIR)
LIBPATH = -L$(LIBDIR)
LIB = -llapack -lblas $(MEMCHECK) -lm

# PHONY targets (not corresponding to files)
.PHONY: all build debug memcheck install install-user uninstall uninstall-user clean dist check help test test-valgrind

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

# Clean generated files (including test artifacts)
clean:
	@echo "Cleaning up..."
	rm -f *.o *.d $(PROGRAM)
	rm -f $(TEST_DIR)/*.o $(TEST_DIR)/*.d $(TEST_RUNNER)
	@echo "Clean complete"

# ============================================================================
# TESTING TARGETS
# ============================================================================

# Build and run unit tests
# Kompiluje testovací suite a spustí všechny testy
test: $(TEST_RUNNER)
	@echo ""
	@echo "Running unit tests..."
	@echo ""
	./$(TEST_RUNNER)
	@echo ""

# Build test runner
# Linkuje testovací program s testovanými moduly
$(TEST_RUNNER): $(TEST_OBJ) $(TEST_MODULES)
	@echo "Linking test runner..."
	$(CC) $(TEST_OBJ) $(TEST_MODULES) $(LIB) -o $(TEST_RUNNER)
	@echo "Test build complete"

# Compile test sources with automatic dependency generation
# Kompiluje testovací soubory s Unity frameworkem
$(TEST_DIR)/%.o: $(TEST_DIR)/%.c
	@echo "Compiling test $<..."
	$(CC) $(WFLAGS) -I. -I$(TEST_DIR) -DUNITY_INCLUDE_DOUBLE -MMD -MP -c $< -o $@

# Include test dependencies
-include $(TEST_OBJ:.o=.d)

# Run tests with Valgrind (memory leak detection)
# Spustí testy s kontrolou memory leaks
test-valgrind: $(TEST_RUNNER)
	@echo "Running tests with Valgrind..."
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./$(TEST_RUNNER)

# ============================================================================

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
	@echo "Testing targets:"
	@echo "  test          : Build and run unit tests"
	@echo "  test-valgrind : Run tests with Valgrind memory checking"
	@echo "  test-clean    : Remove test object files and test runner"
	@echo ""
	@echo "Environment variables:"
	@echo "  LIBDIR        : Library directory (default: $(HOME)/lib)"
	@echo "  INCDIR        : Include directory (default: $(HOME)/include)"
	@echo "  BINDIR        : Binary installation directory (default: $(HOME)/bin)"
	@echo ""
	@echo "Example usage:"
	@echo "  make                  - Build the program"
	@echo "  make debug            - Build with debug information"
	@echo "  make test             - Run unit tests"
	@echo "  make test-valgrind    - Check for memory leaks in tests"
	@echo "  make install-user     - Install to user's bin directory"
	@echo "  LIBDIR=/opt/lib make  - Use custom library path"

# Link the program
$(PROGRAM): $(OBJ) Makefile
	@echo "Linking $(PROGRAM)..."
	$(CC) $(CFLAGS) $(LIBPATH) $(OBJ) $(LIB) -o $(PROGRAM)
	@echo "Build complete"

# Compile source files with automatic dependency generation
%.o: %.c Makefile
	@echo "Compiling $<..."
	$(CC) $(LIBINCLUDE) $(CFLAGS) -MMD -MP -c $< -o $@

# Include automatically generated dependencies (if they exist)
# Note: Explicit dependencies are no longer needed - compiler generates them automatically in .d files
-include $(OBJ:.o=.d)

#
# End
#

