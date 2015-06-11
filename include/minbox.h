// Copyright (c) 2015 Andrew Gainer-Dewar

#ifndef MINBOX_H
#define MINBOX_H

#include <stdexcept>

namespace pmfe {
    template <typename F>
        class MinBox {
        // Specialized container that records the minimum of the insert()ed objects
    public:
        MinBox(){};

        void insert(const F& elt) {
            if (not initialized) {
                initialized = true;
                m_minimum = elt;
            } else {
                if (elt < m_minimum) {
                    m_minimum = elt;
                }
            }
        };

        F minimum() const {
            if (initialized) {
                return m_minimum;
            } else {
                throw std::logic_error("Called minimum() on empty container.");
            }
        };

    protected:
        F m_minimum;
        bool initialized = false;
    };
}

#endif
