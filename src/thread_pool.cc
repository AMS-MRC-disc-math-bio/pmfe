// Copyright (c) 2015 Andrew Gainer-Dewar

#include "thread_pool.h"

#include <boost/asio/io_service.hpp>
#include <boost/thread.hpp>
#include <boost/function.hpp>

namespace pmfe {
    SimpleThreadPool::SimpleThreadPool(size_t num_threads):
        io_service(),
        work(io_service),
        thread_group(),
        incomplete_threads(0)
    {
        if (num_threads == 0) {
            this->num_threads = boost::thread::hardware_concurrency();
        } else {
            this->num_threads = num_threads;
        }

        make_threads();
    };

    void SimpleThreadPool::post(boost::function<void()> job) {
        ++incomplete_threads;

        boost::function<void()> job_with_counter = [this, job] () {
            job();
            --incomplete_threads;
            cv.notify_all();
        };

        io_service.post(job_with_counter);;
    };

    void SimpleThreadPool::wait_for_all_jobs() {
        boost::unique_lock<boost::mutex> lock(mutex);
        cv.wait(lock, [this] () {return incomplete_threads == 0;});
    };

    void SimpleThreadPool::make_threads() {
        for (size_t i = 0; i < num_threads; ++i) {
            thread_group.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
        }
    };
}
