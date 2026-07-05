#!/bin/bash
set -e

echo "Downloading jemalloc 5.3.0..."
wget -qnc https://github.com/jemalloc/jemalloc/releases/download/5.3.0/jemalloc-5.3.0.tar.bz2 || true

echo "Extracting jemalloc..."
tar -xf jemalloc-5.3.0.tar.bz2
cd jemalloc-5.3.0

echo "Configuring jemalloc..."
./configure --prefix=$(pwd)/../jemalloc_install

echo "Building and installing jemalloc..."
make -j$(nproc)
make install

echo "Done! You can now run your CMake build."
