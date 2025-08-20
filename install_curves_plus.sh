#!/bin/bash
# CURVES+ Installation Script for BioStructBenchmark

set -e  # Exit on any error

echo "=== CURVES+ Installation Script ==="
echo "This script will guide you through installing CURVES+ for DNA analysis"
echo

# Check if running on macOS or Linux
OS="$(uname -s)"
case "${OS}" in
    Linux*)     MACHINE=Linux;;
    Darwin*)    MACHINE=Mac;;
    *)          MACHINE="UNKNOWN:${OS}"
esac

echo "Detected OS: ${MACHINE}"

# Create installation directory
INSTALL_DIR="$HOME/curves_plus"
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo "Installation directory: $INSTALL_DIR"

# Download CURVES+ (you'll need to get this manually from the website)
echo
echo "=== DOWNLOAD INSTRUCTIONS ==="
echo "1. Visit: https://bisi.ibcp.fr/tools/curves_plus/"
echo "2. Register and download the CURVES+ package"
echo "3. Extract the files to: $INSTALL_DIR"
echo "4. The package should contain source code and documentation"
echo

# Check if source files exist
if [ ! -f "Makefile" ] && [ ! -f "makefile" ]; then
    echo "ERROR: CURVES+ source files not found in $INSTALL_DIR"
    echo "Please download and extract CURVES+ package first"
    echo
    echo "Expected files:"
    echo "  - Makefile or makefile"
    echo "  - Source files (*.c, *.h)"
    echo "  - Documentation"
    exit 1
fi

echo "=== COMPILATION ==="
echo "Compiling CURVES+ from source..."

# Try different make commands
if [ -f "Makefile" ]; then
    make
elif [ -f "makefile" ]; then
    make -f makefile
else
    echo "ERROR: No Makefile found"
    exit 1
fi

# Check if compilation was successful
if [ ! -f "Cur+" ] && [ ! -f "curves+" ]; then
    echo "ERROR: Compilation failed - no executable found"
    exit 1
fi

# Find the executable
CURVES_EXE=""
if [ -f "Cur+" ]; then
    CURVES_EXE="Cur+"
elif [ -f "curves+" ]; then
    CURVES_EXE="curves+"
fi

echo "✓ Compilation successful: $CURVES_EXE"

# Make executable accessible
echo
echo "=== INSTALLATION ==="
echo "Adding CURVES+ to system PATH..."

# Create symlink in /usr/local/bin (requires sudo)
if [ -w "/usr/local/bin" ]; then
    ln -sf "$INSTALL_DIR/$CURVES_EXE" "/usr/local/bin/Cur+"
    echo "✓ CURVES+ installed to /usr/local/bin/Cur+"
else
    echo "Cannot write to /usr/local/bin"
    echo "Please run with sudo or add to PATH manually:"
    echo "  export PATH=\"$INSTALL_DIR:\$PATH\""
    echo "  echo 'export PATH=\"$INSTALL_DIR:\$PATH\"' >> ~/.bashrc"
fi

# Test installation
echo
echo "=== TESTING ==="
if command -v Cur+ &> /dev/null; then
    echo "✓ CURVES+ is available in PATH"
    Cur+ -h 2>&1 | head -5 || true
else
    echo "⚠ CURVES+ not in PATH. Add manually:"
    echo "  export PATH=\"$INSTALL_DIR:\$PATH\""
fi

echo
echo "=== SETUP COMPLETE ==="
echo "CURVES+ installation directory: $INSTALL_DIR"
echo "Executable: $INSTALL_DIR/$CURVES_EXE"
echo
echo "To use with BioStructBenchmark:"
echo "  from biostructbenchmark.analysis import CurvesAnalyzer"
echo "  analyzer = CurvesAnalyzer('$INSTALL_DIR/$CURVES_EXE')"
echo
echo "For system-wide access, add to your shell profile:"
echo "  echo 'export PATH=\"$INSTALL_DIR:\$PATH\"' >> ~/.bashrc"
echo "  source ~/.bashrc"