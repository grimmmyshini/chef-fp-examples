# Arclength Benchmark

### Build

1. Build the Print Model plugin which is at the base of the repo

2. Run this command to compile the benchmark:
    ```Nix
    clang++ -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang "PATH_TO_CLAD.SO" \
    -Xclang -plugin-arg-clad -Xclang -fcustom-estimation-model                       \
    -Xclang -plugin-arg-clad -Xclang ../PrintModel/libPrintModel.so                  \
    -isystem "PATH_TO_GOOGLE-BENCHMARK_INCLUDE"                                      \
    -isystem "PATH_TO_CoDiPack-1.9_INCLUDE"                                          \
    -isystem "PATH_TO_ADAPT-FP"                                                      \
    -isystem "PATH_TO_CLAD_INCLUDE"                                                  \
    -lstdc++ -lm -O3 -std=c++11                                                      \
    benchmark.cpp -o benchmark.out                                                   \
    -L "PATH_TO_BENCHMARK_BUILD/src"                                                 \
    -lpthread -lbenchmark -DCODI_ZeroAdjointReverse=0
    ```

### Running the benchmark

Execute `benchmark.out`

```Nix
./benchmark.out
```

### Results

Benchmark                                | Time      |  CPU    | Iterations
-----------------------------------------|-----------|---------|-----------
ErrorEstimateArcLenClad/iterations:10    | 0.761 s   | 0.761 s | 10
ErrorEstimateArcLenAdapt/iterations:10   | 10.9 s    |  10.9 s | 10