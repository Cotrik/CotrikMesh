#include <algorithm>
#include <thread>
#include <functional>
#include <vector>

#define NUM_THREADS std::thread::hardware_concurrency() == 0 ? 8 : std::thread::hardware_concurrency()

static void parallelFor(int n, std::function<void (int start, int end)> func) {
    int nThreads = NUM_THREADS;
    int batchSize = n / nThreads;
    std::vector<std::thread> threads(nThreads);
    for (int i = 0; i < nThreads; i++) {
        int start = i * batchSize;
        int end = (start + batchSize) > n ? n - 1 : (start + batchSize);
        threads[i] = std::thread(func, start, end);
    }
    for (auto& t: threads) t.join();
}

#define PARALLEL_FOR_BEGIN(n) parallelFor(n, [&](int start, int end){for (int i = start; i < end; i++)
#define PARALLEL_FOR_END()})
