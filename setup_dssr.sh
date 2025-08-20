#!/bin/bash
# DSSR Installation and Setup Script for BioStructBenchmark

set -e  # Exit on any error

echo "=== DSSR (3DNA Suite) Setup Script ==="
echo "This script will help you set up DSSR for nucleic acid structure analysis"
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
INSTALL_DIR="$HOME/dssr"
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo "Installation directory: $INSTALL_DIR"

echo
echo "=== DSSR LICENSING INFORMATION ==="
echo "DSSR is distributed by Columbia Technology Ventures (CTV)"
echo "Thanks to NIH R24GM153869 grant, DSSR Basic is free for academic use"
echo
echo "Two versions available:"
echo "1. DSSR Basic (Free for Academic Use)"
echo "   - No modeling capabilities"
echo "   - Basic analysis/annotation features"
echo "   - Provided AS IS without support"
echo
echo "2. DSSR Pro (Commercial)"
echo "   - Full modeling capabilities"
echo "   - Complete analysis features"
echo "   - Includes manual and support"
echo

echo "=== DOWNLOAD INSTRUCTIONS ==="
echo "To obtain DSSR:"
echo "1. Visit: https://x3dna.org/"
echo "2. For Academic Users: Submit license request via 'Express Licensing'"
echo "3. For Commercial Users: Email techtransfer@columbia.edu"
echo "4. Download the appropriate binary for your system:"
echo "   - Linux: dssr-linux"
echo "   - macOS: dssr-macos"
echo "   - Windows: dssr-windows.exe"
echo

# Check if user already has DSSR binary
DSSR_BINARY=""
if [ -f "dssr" ]; then
    DSSR_BINARY="dssr"
elif [ -f "dssr-linux" ] && [ "$MACHINE" = "Linux" ]; then
    DSSR_BINARY="dssr-linux"
elif [ -f "dssr-macos" ] && [ "$MACHINE" = "Mac" ]; then
    DSSR_BINARY="dssr-macos"
elif [ -f "dssr-windows.exe" ]; then
    DSSR_BINARY="dssr-windows.exe"
fi

if [ -z "$DSSR_BINARY" ]; then
    echo "ERROR: DSSR binary not found in $INSTALL_DIR"
    echo
    echo "Please download DSSR binary and place it in: $INSTALL_DIR"
    echo "Expected filename based on your OS ($MACHINE):"
    if [ "$MACHINE" = "Linux" ]; then
        echo "  - dssr-linux"
    elif [ "$MACHINE" = "Mac" ]; then
        echo "  - dssr-macos"
    else
        echo "  - dssr (generic) or appropriate binary for your system"
    fi
    echo
    echo "After downloading, re-run this script"
    exit 1
fi

echo "✓ Found DSSR binary: $DSSR_BINARY"

# Rename to standard 'dssr' if needed
if [ "$DSSR_BINARY" != "dssr" ]; then
    cp "$DSSR_BINARY" "dssr"
    echo "✓ Renamed to standard 'dssr' executable"
fi

# Make executable
chmod +x dssr
echo "✓ Made executable"

# Test DSSR
echo
echo "=== TESTING DSSR ==="
if ./dssr --help > /dev/null 2>&1; then
    echo "✓ DSSR executable is working"
    
    # Show version information
    echo "DSSR Version Information:"
    ./dssr --version 2>&1 | head -5 || true
else
    echo "⚠ DSSR executable test failed"
    echo "This might be due to missing dependencies or incompatible binary"
fi

# Create symlink for system-wide access
echo
echo "=== SYSTEM INSTALLATION ==="
SYSTEM_BIN="/usr/local/bin"

if [ -w "$SYSTEM_BIN" ]; then
    ln -sf "$INSTALL_DIR/dssr" "$SYSTEM_BIN/dssr"
    echo "✓ DSSR installed to $SYSTEM_BIN/dssr"
else
    echo "Cannot write to $SYSTEM_BIN"
    echo "To install system-wide, run:"
    echo "  sudo ln -sf $INSTALL_DIR/dssr $SYSTEM_BIN/dssr"
fi

# Add to PATH for current session
export PATH="$INSTALL_DIR:$PATH"

# Update shell profile
SHELL_PROFILE=""
if [ -f "$HOME/.bashrc" ]; then
    SHELL_PROFILE="$HOME/.bashrc"
elif [ -f "$HOME/.zshrc" ]; then
    SHELL_PROFILE="$HOME/.zshrc"
elif [ -f "$HOME/.profile" ]; then
    SHELL_PROFILE="$HOME/.profile"
fi

if [ -n "$SHELL_PROFILE" ]; then
    if ! grep -q "$INSTALL_DIR" "$SHELL_PROFILE" 2>/dev/null; then
        echo >> "$SHELL_PROFILE"
        echo "# DSSR (3DNA suite) - Added by BioStructBenchmark setup" >> "$SHELL_PROFILE"
        echo "export PATH=\"$INSTALL_DIR:\$PATH\"" >> "$SHELL_PROFILE"
        echo "✓ Added DSSR to PATH in $SHELL_PROFILE"
    else
        echo "✓ DSSR already in PATH in $SHELL_PROFILE"
    fi
else
    echo "⚠ Could not determine shell profile to update PATH"
    echo "Please manually add to your shell profile:"
    echo "  export PATH=\"$INSTALL_DIR:\$PATH\""
fi

# Test system installation
echo
echo "=== FINAL VERIFICATION ==="
if command -v dssr &> /dev/null; then
    echo "✓ DSSR is available in PATH"
    
    # Test with a simple command
    echo "Testing DSSR functionality..."
    dssr --version 2>&1 | head -3 || true
else
    echo "⚠ DSSR not found in PATH"
    echo "You may need to:"
    echo "1. Restart your terminal"
    echo "2. Or run: source $SHELL_PROFILE"
    echo "3. Or manually add to PATH: export PATH=\"$INSTALL_DIR:\$PATH\""
fi

echo
echo "=== SETUP COMPLETE ==="
echo "DSSR installation directory: $INSTALL_DIR"
echo "Executable: $INSTALL_DIR/dssr"
echo
echo "To use with BioStructBenchmark:"
echo "  from biostructbenchmark.analysis import DSSRAnalyzer"
echo "  analyzer = DSSRAnalyzer()  # Auto-finds DSSR"
echo "  # OR specify path explicitly:"
echo "  analyzer = DSSRAnalyzer('$INSTALL_DIR/dssr')"
echo
echo "Example usage:"
echo "  results = analyzer.analyze_structure('dna_structure.pdb')"
echo "  interface = analyzer.analyze_protein_dna_interface('complex.pdb')"
echo
echo "Documentation: https://x3dna.org/"
echo "Support: 3DNA Forum at https://x3dna.org/"

# Create a simple test script
cat > "$INSTALL_DIR/test_dssr.py" << 'EOF'
#!/usr/bin/env python3
"""
Simple test script for DSSR installation
"""
import sys
from pathlib import Path

# Add BioStructBenchmark to path (adjust as needed)
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from biostructbenchmark.analysis.dssr import DSSRAnalyzer
    
    print("Testing DSSR integration...")
    
    # Test initialization
    analyzer = DSSRAnalyzer()
    print(f"✓ DSSR found at: {analyzer.dssr_exe}")
    
    print("\nDSSR integration test successful!")
    print("Ready for nucleic acid structure analysis")
    
except Exception as e:
    print(f"✗ DSSR integration test failed: {e}")
    sys.exit(1)
EOF

chmod +x "$INSTALL_DIR/test_dssr.py"
echo
echo "Created test script: $INSTALL_DIR/test_dssr.py"
echo "Run it to verify BioStructBenchmark integration"