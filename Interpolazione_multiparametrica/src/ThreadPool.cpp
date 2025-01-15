#include "ThreadPool.h"

#include <iostream>

ThreadPool::ThreadPool() : stop(false) {
    // Creazione dei thread
    for (unsigned int i = 0; i < numThreads; ++i) {
        workers.emplace_back([this] {
            while (true) {
                std::function<void()> task;

                { // Blocco per sincronizzare l'accesso alla coda
                    std::unique_lock<std::mutex> lock(this->queueMutex);

                    // Attendo che ci siano task disponibili o che il pool venga chiuso
                    this->condition.wait(lock, [this] {
                        return this->stop || !this->tasks.empty();
                        });

                    // Se il pool è stato chiuso e non ci sono più task, termino il thread
                    if (this->stop && this->tasks.empty()) return;

                    // Recupero un task dalla coda
                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }

                // Eseguo il task
                task();

                // Incremento il contatore delle task completate
                tasksCompleted++;
                //std::cout << "->" << tasksCompleted << std::endl;
            }
            });
    }
}


// Funzione per ottenere il numero di task completati
int ThreadPool::getTasksCompleted() const {
    return tasksCompleted.load(); // Carica il valore in modo thread-safe
}


ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queueMutex);      // accesso thread-safe
        stop = true;

        while (!tasks.empty()) {
            tasks.pop();
        }
    }
    condition.notify_all(); // Sveglio tutti i thread per terminare

    for (std::thread& worker : workers) {
        worker.join(); // Aspetto la terminazione dei thread
    }
}
