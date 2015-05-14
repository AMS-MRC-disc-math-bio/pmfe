// Copyright (c) 2015 Andrew Gainer-Dewar

#include "thread_pool.h"

#include <boost/asio/io_service.hpp>
#include <boost/thread.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace pmfe {
    SimpleThreadPool::SimpleThreadPool(size_t num_threads):
        io_service(),
        work(io_service),
        thread_group()
    {
        if (num_threads == 0) {
            this->num_threads = boost::thread::hardware_concurrency();
        } else {
            this->num_threads = num_threads;
        }

        make_threads();
    };

    void SimpleThreadPool::post(boost::function<void()> job) {
        io_service.post(job);
    };

    void SimpleThreadPool::make_threads() {
        for (size_t i = 0; i < num_threads; ++i) {
            thread_group.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
        }
    };

    SimpleJobGroup::SimpleJobGroup(SimpleThreadPool& thread_pool):
        thread_pool(thread_pool) {};

    void SimpleJobGroup::post(boost::function<void()> job) {
        typedef boost::packaged_task<void> task_t;

        boost::shared_ptr<task_t> task = boost::make_shared<task_t>(job);
        boost::unique_future<void> future = task->get_future();
        pending_jobs.push_back(boost::move(future));

        thread_pool.post(boost::bind(&task_t::operator(), task));
    };

    void SimpleJobGroup::wait_for_all_jobs() {
        boost::wait_for_all(pending_jobs.begin(), pending_jobs.end());
    };
}
