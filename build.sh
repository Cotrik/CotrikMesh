if [ ! -d "Release" ]; then
    mkdir "Release"
fi
cd Release
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 8
cd ../

if [ ! -d "Debug" ]; then
    mkdir "Debug"
fi
cd Debug
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j 8

