#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <iostream>
#include <vector>
#include <queue>
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>

#define NUM_THREADS std::thread::hardware_concurrency() == 0 ? 4 : 2 * std::thread::hardware_concurrency()

class ThreadPool {
public:
    ThreadPool() : running(true) {
        int n = NUM_THREADS;
        std::cout << "NUM_THREADS: " << n << std::endl;
        for (int i = 0; i < n; ++i) {
            threads.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock, [this] { return !this->running || !this->queue.empty(); });
                        if (!this->running && this->queue.empty()) return;
                        task = std::move(this->queue.front());
                        this->queue.pop();
                    }
                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            running = false;
        }

        condition.notify_all();

        for (auto& t : threads) {
            t.join();
        }
        std::cout << "Destroyed all threads" << std::endl;
    }	
    
    template <typename Func, typename... Args>
    void submit(Func&& func, Args&&... args) {
        auto boundTask = std::bind(std::forward<Func>(func), std::forward<Args>(args)...);
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            queue.push(boundTask);
        }
        condition.notify_one();
    }

    bool available() {
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            return queue.size() < maxQueueSize;
        }
    }
private:

    std::vector< std::thread > threads;
    std::queue< std::function<void()>> queue;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool running;
    int maxQueueSize = 1000;
};

#endif // THREADPOOL_H
