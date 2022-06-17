# Internal Benchmarks

These files benchmark clad's error estimation versus gradient generation times. To run any of these files, execute the following commands.

```bash
clang -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang /path/to/clad.so -Ipath/to/clad/include -x c++ -lstdc++ -lm -isystem path/to/Google/benchmark/include  -Lpath/to/benchmark/build/src -O2 ./SimpsonsGradVsEE.cpp -lbenchmark -lpthread
```
