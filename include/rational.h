// Copyright (c) 2015 Andrew Gainer-Dewar

#ifndef DPENERGY_H
#define DPENERGY_H

#include <iostream>

#include <gmpxx.h>
#include <CGAL/Gmpq.h>

namespace pmfe {
    typedef mpz_class Integer;

    class Rational {
    public:
    Rational():
        m_value(0),
            m_isFinite(true)
            {};

    Rational(mpq_class energy):
        m_value(energy),
            m_isFinite(true)
            {};

    Rational(Integer num, Integer den = 1):
        m_value(num, den),
            m_isFinite(true)
            {};

    Rational(std::string word):
        m_value(word),
            m_isFinite(true)
            {};

    Rational(int val):
        m_value(val),
            m_isFinite(true)
            {};

        static Rational infinity();

        bool isFinite() const;
        double get_d() const;
        void canonicalize();
        std::string get_str(int base = 10) const;

        Rational& operator+=(const Rational& rhs);
        friend Rational operator+(const Rational& lhs, const Rational& rhs);

        Rational& operator-=(const Rational& rhs);
        friend Rational operator-(const Rational& lhs, const Rational& rhs);

        Rational& operator*=(const Rational& rhs);
        friend Rational operator*(const Rational& lhs, const Rational& rhs);

        Rational& operator/=(const Rational& rhs);
        friend Rational operator/(const Rational& lhs, const Rational& rhs);

        friend bool operator==(const Rational& lhs, const Rational& rhs);
        friend bool operator!=(const Rational& lhs, const Rational& rhs);
        friend bool operator<(const Rational& lhs, const Rational& rhs);
        friend bool operator<=(const Rational& lhs, const Rational& rhs);
        friend bool operator>(const Rational& lhs, const Rational& rhs);
        friend bool operator>=(const Rational& lhs, const Rational& rhs);

        friend std::ostream& operator<<(std::ostream& os, const Rational& energy);

        operator mpq_class() const;
        operator CGAL::Gmpq() const;

    protected:
        mpq_class m_value;
        bool m_isFinite;
    };
}
#endif
