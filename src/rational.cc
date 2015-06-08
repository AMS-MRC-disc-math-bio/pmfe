// Copyright (c) 2015 Andrew Gainer-Dewar

#include <stdexcept>
#include <iostream>

#include <gmpxx.h>
#include <CGAL/Gmpq.h>

#include "rational.h"

namespace pmfe {
    Rational Rational::infinity() {
        Rational result(0);
        result.m_isFinite = false;
        return result;
    }

    bool Rational::isFinite() const {
        return m_isFinite;
    };

    double Rational::get_d() const {
        if (isFinite()) {
            return m_value.get_d();
        } else {
            return std::numeric_limits<double>::infinity();
        }
    };

    void Rational::canonicalize() {
        m_value.canonicalize();
    }

    std::string Rational::get_str(int base) const {
        if (isFinite()) {
            return m_value.get_str(base);
        } else {
            return "∞";
        }
    }

    Rational& Rational::operator+=(const Rational& rhs) {
        if (isFinite() and rhs.isFinite()) {
            m_value += rhs.m_value;
            return *this;
        } else if (not isFinite()) {
            return *this;
        } else {
            m_value = rhs.m_value;
            m_isFinite = rhs.m_isFinite;
            return *this;
        }
    };

    Rational operator+(const Rational& lhs, const Rational& rhs) {
        Rational result = lhs;
        result += rhs;
        return result;
    }

    Rational& Rational::operator-=(const Rational& rhs) {
        if (isFinite() and rhs.isFinite()) {
            m_value -= rhs.m_value;
            return *this;
        } else if (not isFinite() and rhs.isFinite()) {
            return *this;
        } else {
            throw std::logic_error("Invalid arithmetic with ∞.");
        }
    };

    Rational operator-(const Rational& lhs, const Rational& rhs) {
        Rational result = lhs;
        result -= rhs;
        return result;
    }


    Rational& Rational::operator*=(const Rational& rhs) {
        if (isFinite() and rhs.isFinite()) {
            m_value *= rhs.m_value;
            return *this;
        } else if (not isFinite() and rhs != 0) {
            return *this;
        } else if (not rhs.isFinite() and this != 0) {
            m_value = 0;
            m_isFinite = rhs.m_isFinite;
            return *this;
        } else {
            throw std::logic_error("Invalid arithmetic with ∞.");
        }
    };

    Rational operator*(const Rational& lhs, const Rational& rhs) {
        Rational result = lhs;
        result *= rhs;
        return result;
    }

    Rational& Rational::operator/=(const Rational& rhs) {
        if (isFinite() and rhs.isFinite()) {
            m_value /= rhs.m_value;
            return *this;
        } else if (not isFinite() and rhs != 0) {
            return *this;
        } else if (isFinite() and this != 0 and not rhs.isFinite()) {
            m_value = 0;
            return *this;
        } else {
            throw std::logic_error("Invalid arithmetic with ∞.");
        }
    };

    Rational operator/(const Rational& lhs, const Rational& rhs) {
        Rational result = lhs;
        result /= rhs;
        return result;
    };

    bool operator==(const Rational& lhs, const Rational& rhs) {
        if (lhs.isFinite() and rhs.isFinite() and lhs.m_value == rhs.m_value) {
            return true;
        } else if (not lhs.isFinite() and not rhs.isFinite()) {
            return true;
        } else {
            return false;
        }
    };

    bool operator!=(const Rational& lhs, const Rational& rhs) {
        return not (lhs == rhs);
    };

    bool operator<(const Rational& lhs, const Rational& rhs) {
        if (lhs.isFinite() and rhs.isFinite() and lhs.m_value < rhs.m_value) {
            return true;
        } else if (lhs.isFinite() and not rhs.isFinite()) {
            return true;
        } else {
            return false;
        }
    };

    bool operator<=(const Rational& lhs, const Rational& rhs) {
        return (lhs < rhs or lhs == rhs);
    }

    bool operator>(const Rational& lhs, const Rational& rhs) {
        return (rhs < lhs);
    }

    bool operator>=(const Rational& lhs, const Rational& rhs) {
        return (rhs <= lhs);
    }

    std::ostream& operator<<(std::ostream& os, const Rational& energy) {
        if (energy.isFinite())  {
            os << energy.m_value;
        } else {
            os << "∞";
        }

        return os;
    }

    Rational::operator mpq_class() const {
        if (isFinite()) {
            return m_value;
        } else {
            throw std::logic_error("Cannot convert ∞ to rational.");
        }
    }

    Rational::operator CGAL::Gmpq() const {
        if (isFinite()) {
            return CGAL::Gmpq(m_value.get_mpq_t());
        } else {
            throw std::logic_error("Cannot convert ∞ to rational.");
        }
    }
};
