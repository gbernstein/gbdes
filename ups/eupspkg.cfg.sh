build() {
mkdir build
cd build
cmake -S .. -B .
cmake --build . --config Release
}
