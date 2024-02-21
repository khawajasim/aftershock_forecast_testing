#!/bin/bash

# Command file for Sphinx documentation

# Change directory to the directory of the script
cd "$(dirname "$0")"

# Set variables
SPHINXBUILD="sphinx-build"
SOURCEDIR="."
BUILDDIR="_build"

# Check if SPHINXBUILD is set, otherwise set to default
if [ -z "$SPHINXBUILD" ]; then
    SPHINXBUILD="sphinx-build"
fi

# Check if sphinx-build is available
command -v $SPHINXBUILD >/dev/null 2>&1
if [ $? -ne 0 ]; then
    echo ""
    echo "The 'sphinx-build' command was not found. Make sure you have Sphinx"
    echo "installed, then set the SPHINXBUILD environment variable to point"
    echo "to the full path of the 'sphinx-build' executable. Alternatively you"
    echo "may add the Sphinx directory to PATH."
    echo ""
    echo "If you don't have Sphinx installed, grab it from"
    echo "https://www.sphinx-doc.org/"
    echo ""
    exit 1
fi

# Check if command argument is provided
if [ -z "$1" ]; then
    $SPHINXBUILD -M help $SOURCEDIR $BUILDDIR $SPHINXOPTS $O
else
    $SPHINXBUILD -M $1 $SOURCEDIR $BUILDDIR $SPHINXOPTS $O
fi

