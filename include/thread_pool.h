// Copyright (c) 2015 Andrew Gainer-Dewar

#ifndef _THREAD_POOL_H
#define _THREAD_POOL_H

#include <boost/asio/io_service.hpp>
#include <boost/thread.hpp>
#include <boost/function.hpp>
#include <boost/atomic.hpp>
#include <deque>

namespace pmfe {
    class SimpleThreadPool {
    public:
        SimpleThreadPool(size_t num_threads = 0);

        void post(boost::function<void()> job);

        void wait_for_all_jobs();

    protected:
        size_t num_threads;

    private:
        void make_threads();

        boost::asio::io_service io_service;
        boost::asio::io_service::work work;
        boost::thread_group thread_group;

        boost::atomic<unsigned int> incomplete_threads;
        boost::condition_variable cv;
        boost::mutex mutex;
    };

    template<typename T>
        class SimpleThreadSafeDeque {
    public:
        void push_back(const T& data){
            boost::mutex::scoped_lock lock(m_mutex);
            m_deque.push_back(data);
        };

        void push_front(const T& data){
            boost::mutex::scoped_lock lock(m_mutex);
            m_deque.push_front(data);
        };

        bool empty() const {
            boost::mutex::scoped_lock lock(m_mutex);
            return m_deque.empty();
        };

        size_t size() const {
            boost::mutex::scoped_lock lock(m_mutex);
            return m_deque.size();
        };

        T front() const {
            boost::mutex::scoped_lock lock(m_mutex);
            return m_deque.front();
        };

        void pop_front() {
            boost::mutex::scoped_lock lock(m_mutex);
            m_deque.pop_front();
        }

        T pop_and_return_front() {
            boost::mutex::scoped_lock lock(m_mutex);
            T result = m_deque.front();
            m_deque.pop_front();
            return result;
        };

        T back() const {
            boost::mutex::scoped_lock lock(m_mutex);
            return m_deque.back();
        };

        void pop_back() {
            boost::mutex::scoped_lock lock(m_mutex);
            m_deque.pop_back();
        }

        T pop_and_return_back() {
            boost::mutex::scoped_lock lock(m_mutex);
            T result = m_deque.back();
            m_deque.pop_back();
            return result;
        };

    private:
        std::deque<T> m_deque;
        mutable boost::mutex m_mutex;
    };
}

#endif
