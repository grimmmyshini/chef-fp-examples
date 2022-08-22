#include "memory-manager.hpp"

static void BM_malloc_free(benchmark::State& state) {
    int *ptr[10];
    for (auto _ : state) {
        for (int i =0; i < 10; i++) {
            benchmark::DoNotOptimize(ptr[i] = (int *) malloc(10 * sizeof(int *)));
        }
        for (int i =0; i < 10; i++) {
            benchmark::DoNotOptimize(ptr[i]);
            free(ptr[i]);
        }
    }
}

static void BM_new_delete(benchmark::State& state) {
    int *ptr[10];
    for (auto _ : state) {
        for (int i =0; i < 10; i++) {
            benchmark::DoNotOptimize(ptr[i] = new int());
        }
        for (int i =0; i < 10; i++) {
            benchmark::DoNotOptimize(ptr[i]);
            delete ptr[i];
        }
    }
}

static void BM_new_delete_arr(benchmark::State& state) {
    int *ptr[10];
    for (auto _ : state) {
        for (int i =0; i < 10; i++) {
            benchmark::DoNotOptimize(ptr[i] = new int[10]);
        }
        for (int i =0; i < 10; i++) {
            benchmark::DoNotOptimize(ptr[i]);
            delete[] ptr[i];
        }
    }
}

static void BM_empty(benchmark::State& state) {
    int *ptr[10];
    for (auto _ : state) {
        for (int i =0; i < 10; i++) {
            benchmark::DoNotOptimize(ptr[i]);
        }
        for (int i =0; i < 10; i++) {
            benchmark::DoNotOptimize(ptr[i]);
        }
    }
}

BENCHMARK(BM_malloc_free)->Unit(benchmark::kMillisecond)->Iterations(17);
BENCHMARK(BM_new_delete)->Unit(benchmark::kMillisecond)->Iterations(17);
BENCHMARK(BM_new_delete_arr)->Unit(benchmark::kMillisecond)->Iterations(17);
BENCHMARK(BM_empty)->Unit(benchmark::kMillisecond)->Iterations(17);

//BENCHMARK_MAIN();
int main(int argc, char** argv)
{
    ::benchmark::RegisterMemoryManager(mm.get());
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::RegisterMemoryManager(nullptr);
}