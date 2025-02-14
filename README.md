### Build
For building sequential tests
```
cmake -S . -B build -DWITH_MPI=0 -DBUILD_TESTING=1
```
For building parallel implementation only
```
cmake -S . -B build -DWITH_MPI=1
```
Following the invocation of cmake:
```
cd build
make
```
### Run
The size of the matrix needs to be specified as a commandline argument, e.g.

```
./EVP 10000
```

### Testing
```
GTEST_COLOR=1 ctest --test-dir build --output-on-failure -j12
```

Updating submodules (e.g. googletest, GitHub [https://github.com/google/googletest](https://github.com/google/googletest)) code source:
[https://github.com/cpp-for-yourself/lectures-and-homeworks/blob/main/lectures/googletest.md#update-the-submodules-automagically](https://github.com/cpp-for-yourself/lectures-and-homeworks/blob/main/lectures/googletest.md#update-the-submodules-automagically).
