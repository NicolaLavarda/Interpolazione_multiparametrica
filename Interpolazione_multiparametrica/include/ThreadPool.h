#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>
#include <future>

class ThreadPool {
public:
    ThreadPool();

    template <class F, class... Args>       // Per le classi template devono essere dichiarate ed implementate nello stesso file
    std::future<typename std::result_of<F(Args...)>::type> enqueue(F&& f, Args&&... args) {
        using return_type = typename std::result_of<F(Args...)>::type;

        // Crea una promessa per restituire il risultato del task
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        // Creazione di una future che verrà usata per ottenere il risultato
        std::future<return_type> res = task->get_future();

        {
            std::unique_lock<std::mutex> lock(queueMutex);

            // Accodo il nuovo task
            tasks.emplace([task]() { (*task)(); });
        }

        condition.notify_one(); // Notifico un thread che c'è un nuovo task

        return res; // Restituisco il futuro
    }

    ~ThreadPool();

private:
    unsigned int numThreads = std::thread::hardware_concurrency();

    std::vector<std::thread> workers;                // Vettore di thread
    std::queue<std::function<void()>> tasks;         // Coda di task
    std::mutex queueMutex;                           // Mutex per proteggere la coda
    std::condition_variable condition;               // Variabile di condizione per sincronizzazione
    std::atomic<bool> stop;                          // Flag per terminare il pool
};

#endif
