#ifndef __MEMORY_MANAGER_HPP
#define __MEMORY_MANAGER_HPP

#include <memory>
#include <benchmark/benchmark.h>

class CustomMemoryManager : public benchmark::MemoryManager
{
public:
    int64_t num_allocs;
    int64_t max_bytes_used;
    int64_t total_allocated_bytes;
    int64_t net_heap_growth;

    void Start() BENCHMARK_OVERRIDE
    {
        num_allocs = 0;
        max_bytes_used = 0;
        total_allocated_bytes = 0;
        net_heap_growth = 0;
    }

    void Stop(Result *result) BENCHMARK_OVERRIDE
    {
        result->num_allocs = num_allocs;
        result->max_bytes_used = max_bytes_used;
        result->total_allocated_bytes = total_allocated_bytes;
        result->net_heap_growth = net_heap_growth;
    }
};

std::unique_ptr<CustomMemoryManager> mm(new CustomMemoryManager());

#ifdef MEMORY_PROFILER
void *custom_malloc(size_t size)
{
    size_t *p = (size_t *)malloc(sizeof(size_t) + size);
    auto mm_ptr = mm.get();
    if (mm_ptr)
    {
        mm_ptr->num_allocs += 1;
        mm_ptr->total_allocated_bytes += size;
        int64_t &net_heap_growth = mm_ptr->net_heap_growth;
        int64_t &max_bytes_used = mm_ptr->max_bytes_used;
        net_heap_growth += size;
        max_bytes_used = net_heap_growth > max_bytes_used ? net_heap_growth : max_bytes_used;
        p[0] = size;
    }
    return (void *)(p + 1);
}

void custom_free(void *ptr)
{
    auto mm_ptr = mm.get();
    if (mm_ptr)
    {
        size_t *p = static_cast<size_t *>(ptr) - 1;
        size_t p_size = p[0];
        mm_ptr->net_heap_growth -= p_size;
        ptr = (void *) p;
    }
    free(ptr);
}
#define malloc(size) custom_malloc(size)
#define free(ptr) custom_free(ptr)

void *operator new(size_t size)
{
    void *p = malloc(size);
    return p;
}

void operator delete(void *p) noexcept
{
    free(p);
}

void *operator new[](size_t size)
{
    void *p = malloc(size);
    return p;
}

void operator delete[](void *p) noexcept
{
    free(p);
}
#endif

#endif // __MEMORY_MANAGER_HPP