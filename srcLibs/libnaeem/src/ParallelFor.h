#include <iostream>
#include <algorithm>
#include <thread>
#include <functional>
#include <vector>
#include <mutex>
#include <math.h>

#define NUM_THREADS std::thread::hardware_concurrency() == 0 ? 8 : std::thread::hardware_concurrency() * 2

static void parallelFor(int n, std::function<void (int start, int end)> func) {
    int nThreads = NUM_THREADS;
    // int batchSize = ceil((float) n / (float) nThreads);
    int batchSize = (n / nThreads) + 1;
    if (n <= nThreads) {
        nThreads = 1;
        batchSize = n;
    }
    // std::vector<std::thread> threads(nThreads);
    std::vector<std::thread> threads;
    for (int i = 0; i < nThreads; i++) {
        int start = i * batchSize;
        int end = (start + batchSize) > n ? n : (start + batchSize);
        threads.push_back(std::thread(func, start, end));
        if (end >= n) break;
        // threads[i] = std::thread(func, start, end);
    }
    for (auto& t: threads) t.join();
}

#define PARALLEL_FOR_BEGIN(n) parallelFor(n, [&](int start, int end){for (int i = start; i < end; i++)
#define PARALLEL_FOR_END()})