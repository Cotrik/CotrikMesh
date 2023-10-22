#include <iostream>
#include <algorithm>
#include <thread>
#include <functional>
#include <vector>
#include <mutex>
#include <math.h>

#define NUM_THREADS std::thread::hardware_concurrency() == 0 ? 4 : std::thread::hardware_concurrency()

static void parallelFor(int nStart, int nEnd, std::function<void (int start, int end)> func) {
    int nThreads = NUM_THREADS;
    int n = nEnd;
    int batchSize = (n / nThreads) + 1;
    if (n <= nThreads) {
        nThreads = 1;
        batchSize = n;
    }
    std::cout << "num threads: " << nThreads << std::endl;
    std::cout << "batch size: " << batchSize << std::endl;
    std::vector<std::thread> threads;
    for (int i = 0; i < nThreads; i++) {
        int start = (i * batchSize) + nStart;
        int end = (start + batchSize) > n ? n : (start + batchSize);
        threads.push_back(std::thread(func, start, end));
        if (end >= n) break;
    }
    for (auto& t: threads) t.join();
}

#define PARALLEL_FOR_BEGIN(nStart, nEnd) parallelFor(nStart, nEnd, [&](int start, int end){for (int i = start; i < end; i++)
#define PARALLEL_FOR_END()})
