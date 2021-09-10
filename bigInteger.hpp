//
// Created by Freddy Don on 2019/5/26.
//

#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <cstdint>
#include <list>
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <stdexcept>
#include <random>
#include <compare>
#include <map>

class bigInteger final {
private:
    std::list<std::uint32_t> arr;
    bool negative;

    inline bool absEquals(const bigInteger &o) const {
        return arr == o.arr;
    }

    inline bool absLessThan(const bigInteger &o) const {
        if (arr.size() != o.arr.size())
            return arr.size() < o.arr.size();
        for (auto iter = arr.cbegin(), oiter = o.arr.cbegin(); arr.cend() != iter; ++iter, ++oiter)
            if (*iter != *oiter)
                return *iter < *oiter;
        return false;
    }

    inline bool absGreaterThan(const bigInteger &o) const {
        if (arr.size() != o.arr.size())
            return arr.size() > o.arr.size();
        for (auto iter = arr.cbegin(), oiter = o.arr.cbegin(); arr.cend() != iter; ++iter, ++oiter)
            if (*iter != *oiter)
                return *iter > *oiter;
        return false;
    }

    inline void absAdd(const bigInteger &o) {
        if (o.isZero())
            return;
        while (arr.size() < o.arr.size())
            arr.push_front(0);
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        auto iter = arr.rbegin();
        for (auto oiter = o.arr.crbegin(); o.arr.crend() != oiter; ++iter, ++oiter)
            *iter = static_cast<std::uint32_t>(u64 = static_cast<std::uint64_t>(*iter) + *oiter + u32[1]);
        for (; arr.rend() != iter && u32[1] > 0; ++iter)
            *iter = static_cast<std::uint32_t>(u64 = static_cast<std::uint64_t>(*iter) + u32[1]);
        if (u32[1] > 0)
            arr.push_front(u32[1]);
    }

    inline void absSub(const bigInteger &o) {
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        auto iter = arr.rbegin();
        for (auto oiter = o.arr.crbegin(); o.arr.rend() != oiter; ++iter, ++oiter)
            *iter = static_cast<std::uint32_t>(u64 = static_cast<std::uint64_t>(*iter) - *oiter - (-u32[1]));
        for (; arr.rend() != iter && u32[1] != 0; ++iter)
            *iter = static_cast<std::uint32_t>(u64 = static_cast<std::uint64_t>(*iter) - (-u32[1]));
        while (!arr.empty() && arr.front() == 0)
            arr.pop_front();
        negative &= !arr.empty();
    }

    inline void absMul(const std::uint32_t i) {
        if (arr.empty())
            return;
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        for (auto iter = arr.rbegin(); iter != arr.rend(); ++iter)
            *iter = static_cast<std::uint32_t>(u64 = static_cast<std::uint64_t>(*iter) * i + u32[1]);
        if (0 != u32[1])
            arr.push_front(u32[1]);
        negative &= !arr.empty();
    }

    inline void absBitAndSelf(const bigInteger &o, std::uint32_t reverse = 0) {
        while (arr.size() > o.arr.size())
            arr.pop_front();
        auto iter = arr.rbegin();
        for (auto oiter = o.arr.crbegin(); arr.rend() != iter; *iter++ &= *oiter++ ^ reverse);
        while (!arr.empty() && 0 == arr.front())
            arr.pop_front();
        negative &= !arr.empty();
    }

    inline void absBitOrSelf(const bigInteger &o) {
        while (arr.size() < o.arr.size())
            arr.push_front(0);
        auto iter = arr.rbegin();
        for (auto oiter = o.arr.crbegin(); o.arr.rend() != oiter; *iter++ |= *oiter++);
    }

    inline void divisionRightShiftBooster(std::size_t step = 1) {
        bool carryBit = false;
        union {
            std::uint64_t u64;
            std::uint32_t u32[2];
        };
        ++step;
        auto iter = arr.begin();
        while (iter != arr.end()) {
            u64 = 0;
            u32[1] = *iter;
            u64 >>= 1;
            *iter = u32[1];
            if (carryBit)
                *iter |= 0x80000000;
            carryBit = u32[0] > 0;
            ++iter;
            --step;
            if (0 == step)
                break;
        }
        if (!arr.empty() && 0 == arr.front())
            arr.pop_front();
        negative &= !isZero();
    }

    inline void divisionSubBooster(const bigInteger &o, std::size_t step) {
        auto iter = arr.begin();
        auto oiter = o.arr.begin();
        if (arr.size() > o.arr.size())
            ++iter;
        while (step--) {
            ++iter;
            ++oiter;
        }
        if (iter != arr.end()) {
            ++iter;
            ++oiter;
        }
        union {
            std::uint64_t u64;
            std::uint32_t u32[2];
        };
        u64 = 0x0000000100000000;
        do {
            --iter;
            --oiter;
            if (u32[1] == 0)
                u64 = 0x00000000ffffffff;
            else
                u64 = 0x0000000100000000;
            u64 += (*iter);
            u64 -= *oiter;
            *iter = u32[0];
        } while (oiter != o.arr.begin());
        if (u32[1] == 0) {
            --iter;
            while (0 == *iter) {
                *iter = 0xffffffff;
                --iter;
            }
            (*iter)--;
        }
        while (!isZero() && arr.front() == 0)
            arr.pop_front();
    }

public:
    bigInteger() = default;

    ~bigInteger() = default;

    bigInteger(const bigInteger &) = default;

    bigInteger(bigInteger &&) noexcept(std::is_nothrow_move_constructible_v<decltype(arr)>) = default;
    
    template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
    bigInteger(T t) : arr(), negative(false) {
        if constexpr (std::is_same_v<std::remove_cv_t<T>, bool>) {
            if (t)
                arr.push_back(1);
        } else {
            T temp = t;
            std::size_t digits = std::numeric_limits<T>::digits;
            if constexpr (std::numeric_limits<T>::is_signed) {
                digits += 1;
                if (t < 0) {
                    negative = true;
                    temp = -temp;
                }
            }
            std::size_t i = 32;
            for (; i <= digits; i += 32) {
                arr.push_front((uint32_t) temp);
                if constexpr (std::numeric_limits<T>::digits >= 32) {
                    temp >>= 31;
                    temp >>= 1;
                }
            }
            if (digits % 32) {
                uint32_t inserted = 1 << (digits % 32 - 1);
                inserted |= inserted >> 1;
                inserted |= inserted >> 2;
                inserted |= inserted >> 4;
                inserted |= inserted >> 8;
                inserted |= inserted >> 16;
                inserted &= (uint32_t) temp;
                arr.push_front(inserted);
            }
            while (!arr.empty() && 0 == arr.front())
                arr.pop_front();
        }
    }
    
    template<typename charT = char, typename Traits = std::char_traits<charT>, typename Allocator = std::allocator<charT>>
    bigInteger(std::basic_stringstream<charT, Traits, Allocator> in) : arr(), negative(false) {
        if (!in.good())
            return;
        bool isNegative = false;
        while (' ' == in.peek() || '\t' == in.peek()) in.get();
        if ('-' == in.peek()) {
            in.get();
            isNegative = true;
        } else if ('+' == in.peek()) {
            in.get();
        }
        if ('0' == in.peek()) {
            in.get();
            std::list<bool> l;
            if ('x' == in.peek() || 'X' == in.peek()) {
                in.get();
                while ('_' == in.peek() || ('0' <= in.peek() && '9' >= in.peek()) || ('a' <= in.peek() && 'f' >= in.peek()) || ('A' <= in.peek() && 'F' >= in.peek())) {
                    switch (in.peek()) {
                        case '0':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '1':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case '2':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case '3':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        case '4':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '5':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case '6':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case '7':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        case '8':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '9':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case 'a':
                        case 'A':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case 'b':
                        case 'B':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        case 'c':
                        case 'C':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case 'd':
                        case 'D':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case 'e':
                        case 'E':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case 'f':
                        case 'F':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        default:
                            break;
                    }
                    in.get();
                }
            } else {
                while ('_' == in.peek() || ('0' <= in.peek() && '7' >= in.peek())) {
                    switch (in.peek()) {
                        case '0':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '1':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case '2':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case '3':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        case '4':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '5':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case '6':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case '7':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        default:
                            break;
                    }
                    in.get();
                }
            }
            while (!l.empty() && false == l.front()) l.pop_front();
            while (0 != (l.size() & 31)) l.push_front(false);
            while (!l.empty()) {
                std::uint32_t temp = 0;
                for (int i = 31; i >= 0; --i) {
                    if (l.front())
                        temp |= 1 << i;
                    l.pop_front();
                }
                arr.push_back(temp);
            }
        } else {
            while ('_' == in.peek() || ('0' <= in.peek() && '9' >= in.peek())) {
                if ('_' == in.peek())
                    continue;
                absMul(10);
                absAdd(in.peek() - '0');
                in.get();
            }
        }
        negative = isNegative && !isZero();
    }
    
    template<typename charT = char, typename Traits = std::char_traits<charT>, typename Allocator = std::allocator<charT>>
    bigInteger(const std::basic_string<charT, Traits, Allocator> &s) : bigInteger(std::basic_stringstream<charT, Traits, Allocator>(s)) {}
    
    template<typename charT = char, typename Traits = std::char_traits<charT>, typename Allocator = std::allocator<charT>>
    bigInteger(const charT * s) : bigInteger(std::basic_string<charT, Traits, Allocator>(s)) {}

    bigInteger &operator=(const bigInteger &) & = default;

    bigInteger &operator=(bigInteger &&) & noexcept(std::is_nothrow_move_assignable_v<decltype(arr)>) = default;

    friend bigInteger &operator+=(bigInteger &left, const bigInteger &right) {
        if (left.negative ^ right.negative) {
            if (left.absLessThan(right)) {
                bigInteger newBigInteger(right);
                newBigInteger.absSub(left);
                left = std::move(newBigInteger);
            } else
                left.absSub(right);
        } else
            left.absAdd(right);
        return left;
    }

    friend bigInteger &operator+=(bigInteger &left, bigInteger &&right) {
        if (left.negative ^ right.negative) {
            if (left.absLessThan(right)) {
                right.absSub(left);
                left = std::move(right);
            } else
                left.absSub(right);
        } else
            left.absAdd(right);
        return left;
    }
    
    friend bigInteger &&operator+(bigInteger &&left, const bigInteger &right) {
        if (left.negative ^ right.negative) {
            if (left.absLessThan(right)) {
                bigInteger newBigInteger(right);
                newBigInteger.absSub(left);
                left = std::move(newBigInteger);
            } else
            left.absSub(right);
        }
        else
            left.absAdd(right);
        return static_cast<bigInteger &&>(left);
    }

    friend bigInteger &&operator+(bigInteger &&left, bigInteger &&right) {
        if (left.negative ^ right.negative) {
            if (left.absLessThan(right)) {
                right.absSub(left);
                left = std::move(right);
            } else
                left.absSub(right);
        } else
            left.absAdd(right);
        return static_cast<bigInteger &&>(left);
    }

    friend bigInteger operator+(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) + right;
    }

    friend bigInteger operator+(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) + std::move(right);
    }

    friend bigInteger &operator-=(bigInteger &left, const bigInteger &right) {
        if (&left == &right) {
            left.negative = false;
            left.arr.clear();
            return left;
        }
        if (left.negative ^ right.negative)
            left.absAdd(right);
        else if (left.absLessThan(right)) {
            bigInteger newBigInteger(right);
            newBigInteger.absSub(left);
            left = std::move(newBigInteger);
            left.negative ^= true;
        } else
            left.absSub(right);
        return left;
    }
    
    friend bigInteger &operator-=(bigInteger &left, bigInteger &&right) {
        if (left.negative ^ right.negative)
            left.absAdd(right);
        else if (left.absLessThan(right)) {
            right.absSub(left);
            left = std::move(right);
            left.negative ^= true;
        } else
            left.absSub(right);
        return left;
    }

    friend bigInteger &&operator-(bigInteger &&left, const bigInteger &right) {
        if (left.negative ^ right.negative)
            left.absAdd(right);
        else if (left.absLessThan(right)) {
            bigInteger newBigInteger(right);
            newBigInteger.absSub(left);
            left = std::move(newBigInteger);
            left.negative ^= true;
        } else
            left.absSub(right);
        return static_cast<bigInteger &&>(left);
    }

    friend bigInteger &&operator-(bigInteger &&left, bigInteger &&right) {
        if (left.negative ^ right.negative)
            left.absAdd(right);
        else if (left.absLessThan(right)) {
            right.absSub(left);
            left = std::move(right);
            left.negative ^= true;
        } else
            left.absSub(right);
        return static_cast<bigInteger &&>(left);
    }
    
    friend bigInteger operator-(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) - right;
    }
    
    friend bigInteger operator-(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) - std::move(right);
    }

    friend bigInteger operator-(const bigInteger &o) {
        return -bigInteger(o);
    }

    friend bigInteger &&operator-(bigInteger &&o) {
        o.negative ^= !o.arr.empty();
        return static_cast<bigInteger &&>(o);
    }

    friend const bigInteger &operator+(const bigInteger &o) {
        return o;
    }

    friend bigInteger &&operator+(bigInteger &&o) {
        return static_cast<bigInteger &&>(o);
    }

    friend bigInteger &operator*=(bigInteger &left, const bigInteger &right) {
        bigInteger ret, temp;
        if (right.arr.empty() || left.arr.empty()) {
            left.negative = false;
            left.arr.clear();
            return left;
        }
        for (std::uint32_t num : right.arr) {
            ret.arr.push_back(0);
            temp = left;
            temp.absMul(num);
            ret.absAdd(temp);
        }
        ret.negative = left.negative ^ right.negative;
        return left = std::move(ret);
    }

    friend bigInteger &operator*=(bigInteger &left, bigInteger &&right) {
        bigInteger ret, temp;
        if (right.arr.empty() || left.arr.empty()) {
            left.negative = false;
            left.arr.clear();
            return left;
        }
        for (std::uint32_t num : right.arr) {
            ret.arr.push_back(0);
            temp = left;
            temp.absMul(num);
            ret.absAdd(temp);
        }
        ret.negative = left.negative ^ right.negative;
        return left = std::move(ret);
    }
    
    friend bigInteger &&operator*(bigInteger &&left, const bigInteger &right) {
        bigInteger ret, temp;
        if (right.arr.empty() || left.arr.empty()) {
            left.negative = false;
            left.arr.clear();
            return static_cast<bigInteger &&>(left);
        }
        for (std::uint32_t num : right.arr) {
            ret.arr.push_back(0);
            temp = left;
            temp.absMul(num);
            ret.absAdd(temp);
        }
        ret.negative = left.negative ^ right.negative;
        return static_cast<bigInteger &&>(left = std::move(ret));
    }

    friend bigInteger &&operator*(bigInteger &&left, bigInteger &&right) {
        bigInteger ret, temp;
        if (right.arr.empty() || left.arr.empty()) {
            left.negative = false;
            left.arr.clear();
            return static_cast<bigInteger &&>(left);
        }
        for (std::uint32_t num : right.arr) {
            ret.arr.push_back(0);
            temp = left;
            temp.absMul(num);
            ret.absAdd(temp);
        }
        ret.negative = left.negative ^ right.negative;
        return static_cast<bigInteger &&>(left = std::move(ret));
    }
    
    friend bigInteger operator*(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) * right;
    }
    
    friend bigInteger operator*(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) * std::move(right);
    }

    friend bigInteger &operator/=(bigInteger &left, const bigInteger &right) {
        if (right.arr.empty())
            throw std::overflow_error("divided zero");
        if (&left == &right) {
            left.negative = false;
            left.arr.clear();
            left.arr.push_front(1);
            return left;
        }
        left.negative ^= right.negative;
        auto step = right.arr.size();
        bigInteger pos(1), divisor(right.abs()), divided(left.abs());
        left.arr.clear();
        while (divisor.absLessThan(divided)) {
            divisor.arr.push_back(0);
            pos.arr.push_back(0);
        }
        while (left.arr.size() < pos.arr.size())
            left.arr.push_front(0);
        auto it = left.arr.begin();
        for (; !pos.arr.empty(); pos.divisionRightShiftBooster(), divisor.divisionRightShiftBooster(step)) {
            if (!divisor.absGreaterThan(divided)) {
                *it |= pos.arr.front();
                divided.divisionSubBooster(divisor, step);
            }
            if (1 == pos.arr.front())
                ++it;
        }
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        left.negative &= !left.isZero();
        return left;
    }

    friend bigInteger &operator/=(bigInteger &left, bigInteger &&right) {
        if (right.arr.empty())
            throw std::overflow_error("divided zero");
        left.negative ^= right.negative;
        std::size_t step = right.arr.size();
        bigInteger pos(1), divisor(static_cast<bigInteger &&>(right).abs()), divided(left.abs());
        left.arr.clear();
        while (divisor.absLessThan(divided)) {
            divisor.arr.push_back(0);
            pos.arr.push_back(0);
        }
        while (left.arr.size() < pos.arr.size())
            left.arr.push_front(0);
        auto it = left.arr.begin();
        for (; !pos.arr.empty(); pos.divisionRightShiftBooster(), divisor.divisionRightShiftBooster(step)) {
            if (!divisor.absGreaterThan(divided)) {
                *it |= pos.arr.front();
                divided.divisionSubBooster(divisor, step);
            }
            if (1 == pos.arr.front())
                ++it;
        }
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        left.negative &= !left.isZero();
        return left;
    }
    
    friend bigInteger &&operator/(bigInteger &&left, const bigInteger &right) {
        if (right.arr.empty())
            throw std::overflow_error("divided zero");
        left.negative ^= right.negative;
        std::size_t step = right.arr.size();
        bigInteger pos(1), divisor(right.abs()), divided(static_cast<bigInteger &&>(left).abs());
        left.arr.clear();
        while (divisor.absLessThan(divided)) {
            divisor.arr.push_back(0);
            pos.arr.push_back(0);
        }
        while (left.arr.size() < pos.arr.size())
            left.arr.push_front(0);
            auto it = left.arr.begin();
            for (; !pos.arr.empty(); pos.divisionRightShiftBooster(), divisor.divisionRightShiftBooster(step)) {
                if (!divisor.absGreaterThan(divided)) {
                    *it |= pos.arr.front();
                    divided.divisionSubBooster(divisor, step);
                }
                if (1 == pos.arr.front())
                    ++it;
            }
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
            left.negative &= !left.isZero();
            return static_cast<bigInteger &&>(left);
    }

    friend bigInteger &&operator/(bigInteger &&left, bigInteger &&right) {
        if (right.arr.empty())
            throw std::overflow_error("divided zero");
        left.negative ^= right.negative;
        std::size_t step = right.arr.size();
        bigInteger pos(1), divisor(static_cast<bigInteger &&>(right).abs()), divided(static_cast<bigInteger &&>(left).abs());
        left.arr.clear();
        while (divisor.absLessThan(divided)) {
            divisor.arr.push_back(0);
            pos.arr.push_back(0);
        }
        while (left.arr.size() < pos.arr.size())
            left.arr.push_front(0);
        auto it = left.arr.begin();
        for (; !pos.arr.empty(); pos.divisionRightShiftBooster(), divisor.divisionRightShiftBooster(step)) {
            if (!divisor.absGreaterThan(divided)) {
                *it |= pos.arr.front();
                divided.divisionSubBooster(divisor, step);
            }
            if (1 == pos.arr.front())
                ++it;
        }
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        left.negative &= !left.isZero();
        return static_cast<bigInteger &&>(left);
    }
    
    friend bigInteger operator/(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) / right;
    }
    
    friend bigInteger operator/(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) / std::move(right);
    }

    friend bigInteger &operator%=(bigInteger &left, const bigInteger &right) {
        if (right.arr.empty())
            throw std::overflow_error("divided zero");
        if (&left == &right) {
            left.negative = false;
            left.arr.clear();
            return left;
        }
        std::size_t step = right.arr.size();
        bigInteger pos(1), divisor(right.abs());
        while (divisor.absLessThan(left)) {
            divisor.arr.push_back(0);
            pos.arr.push_back(0);
        }
        for (; !pos.arr.empty(); pos.divisionRightShiftBooster(), divisor.divisionRightShiftBooster(step))
            if (!left.absLessThan(divisor))
                left.divisionSubBooster(divisor, step);
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        left.negative &= !left.isZero();
        return left;
    }

    friend bigInteger &operator%=(bigInteger &left, bigInteger &&right) {
        if (right.arr.empty())
            throw std::overflow_error("divided zero");
        std::size_t step = right.arr.size();
        bigInteger pos(1);
        while (right.absLessThan(left)) {
            right.arr.push_back(0);
            pos.arr.push_back(0);
        }
        for (; !pos.arr.empty(); pos.divisionRightShiftBooster(), right.divisionRightShiftBooster(step))
            if (!left.absLessThan(right))
                left.divisionSubBooster(right, step);
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        left.negative &= !left.isZero();
        return left;
    }
    
    friend bigInteger &&operator%(bigInteger &&left, const bigInteger &right) {
        if (right.arr.empty())
            throw std::overflow_error("divided zero");
        std::size_t step = right.arr.size();
        bigInteger pos(1), divisor(right.abs());
        while (divisor.absLessThan(left)) {
            divisor.arr.push_back(0);
            pos.arr.push_back(0);
        }
        for (; !pos.arr.empty(); pos.divisionRightShiftBooster(), divisor.divisionRightShiftBooster(step))
            if (!left.absLessThan(divisor))
                left.divisionSubBooster(divisor, step);
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        left.negative &= !left.isZero();
        return static_cast<bigInteger &&>(left);
    }
    
    friend bigInteger &&operator%(bigInteger &&left, bigInteger &&right) {
        if (right.arr.empty())
            throw std::overflow_error("divided zero");
        std::size_t step = right.arr.size();
        bigInteger pos(1);
        while (right.absLessThan(left)) {
            right.arr.push_back(0);
            pos.arr.push_back(0);
        }
        for (; !pos.arr.empty(); pos.divisionRightShiftBooster(), right.divisionRightShiftBooster(step))
            if (!left.absLessThan(right))
                left.divisionSubBooster(right, step);
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        left.negative &= !left.isZero();
        return static_cast<bigInteger &&>(left);
    }
    
    friend bigInteger operator%(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) % right;
    }
    
    friend bigInteger operator%(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) % std::move(right);
    }

    friend bigInteger &operator>>=(bigInteger &left, const bigInteger &right) {
        if (right.arr.empty() || left.arr.empty())
            return left;
        if (right.negative)
            throw std::logic_error("left shift count is a negative number");
        if (left.negative)
            left.absSub(1);
        bigInteger times32(right);
        for (; !times32.absLessThan(32) && !left.isZero(); times32.absSub(32), left.arr.pop_back());
        if (times32.isZero() || left.isZero()) {
            if (left.negative)
                left.absAdd(1);
            return left;
        }
        std::uint32_t times1 = times32.arr.front() & 31;
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        for (std::uint32_t &num : left.arr) {
            u64 = (static_cast<std::uint64_t>(num) << (32 - times1)) | (static_cast<std::uint64_t>(u32[0]) << 32);
            num = u32[1];
        }
        if (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        if (left.negative)
            left.absAdd(1);
        return left;
    }

    friend bigInteger &operator>>=(bigInteger &left, bigInteger &&right) {
        if (right.arr.empty() || left.arr.empty())
            return left;
        if (right.negative)
            throw std::logic_error("left shift count is a negative number");
        if (left.negative)
            left.absSub(1);
        for (; !right.absLessThan(32) && !left.isZero(); right.absSub(32), left.arr.pop_back());
        if (right.isZero() || left.isZero()) {
            if (left.negative)
                left.absAdd(1);
            return left;
        }
        std::uint32_t times1 = *right.arr.begin() & 31;
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        for (std::uint32_t &num : left.arr) {
            u64 = (static_cast<std::uint64_t>(num) << (32 - times1)) | (static_cast<std::uint64_t>(u32[0]) << 32);
            num = u32[1];
        }
        if (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        if (left.negative)
            left.absAdd(1);
        return left;
    }
    
    friend bigInteger &&operator>>(bigInteger &&left, const bigInteger &right) {
        if (right.arr.empty() || left.arr.empty())
            return static_cast<bigInteger &&>(left);
        if (right.negative)
            throw std::logic_error("left shift count is a negative number");
        if (left.negative)
            left.absSub(1);
            bigInteger times32(right);
            for (; !times32.absLessThan(32) && !left.isZero(); times32.absSub(32), left.arr.pop_back());
        if (times32.isZero() || left.isZero()) {
            if (left.negative)
                left.absAdd(1);
            return static_cast<bigInteger &&>(left);
        }
        std::uint32_t times1 = *times32.arr.begin() & 31;
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        for (std::uint32_t &num : left.arr) {
            u64 = (static_cast<std::uint64_t>(num) << (32 - times1)) | (static_cast<std::uint64_t>(u32[0]) << 32);
            num = u32[1];
        }
        if (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        if (left.negative)
            left.absAdd(1);
        return static_cast<bigInteger &&>(left);
    }

    friend bigInteger &&operator>>(bigInteger &&left, bigInteger &&right) {
        if (right.arr.empty() || left.arr.empty())
            return static_cast<bigInteger &&>(left);
        if (right.negative)
            throw std::logic_error("left shift count is a negative number");
        if (left.negative)
            left.absSub(1);
        for (; !right.absLessThan(32) && !left.isZero(); right.absSub(32), left.arr.pop_back());
        if (right.isZero() || left.isZero()) {
            if (left.negative)
                left.absAdd(1);
            return static_cast<bigInteger &&>(left);
        }
        std::uint32_t times1 = *right.arr.begin() & 31;
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        for (std::uint32_t &num : left.arr) {
            u64 = (static_cast<std::uint64_t>(num) << (32 - times1)) | (static_cast<std::uint64_t>(u32[0]) << 32);
            num = u32[1];
        }
        if (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        if (left.negative)
            left.absAdd(1);
        return static_cast<bigInteger &&>(left);
    }
    
    friend bigInteger operator>>(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) >> right;
    }
    
    friend bigInteger operator>>(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) >> std::move(right);
    }

    friend bigInteger &operator<<=(bigInteger &left, const bigInteger &right) {
        if (right.arr.empty() || left.arr.empty())
            return left;
        if (right.negative)
            throw std::logic_error("left shift count is a negative number");
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        std::uint32_t times1 = right.arr.back() & 31;
        for (auto iter = left.arr.rbegin(); iter != left.arr.rend(); *iter++ = static_cast<std::uint32_t>(u64 = static_cast<std::uint64_t>(*iter) << times1 | u32[1]));
        if (0 != u32[1])
            left.arr.push_front(u32[1]);
        for (bigInteger times32(right); !times32.absLessThan(32); times32.absSub(32), left.arr.push_back(0));
        return left;
    }

    friend bigInteger &operator<<=(bigInteger &left, bigInteger &&right) {
        if (right.arr.empty() || left.arr.empty())
            return left;
        if (right.negative)
            throw std::logic_error("left shift count is a negative number");
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        std::uint32_t times1 = right.arr.back() & 31;
        for (auto iter = left.arr.rbegin(); iter != left.arr.rend(); *iter++ = static_cast<std::uint32_t>(u64 = static_cast<std::uint64_t>(*iter) << times1 | u32[1]));
        if (0 != u32[1])
            left.arr.push_front(u32[1]);
        for (; !right.absLessThan(32); right.absSub(32), left.arr.push_back(0));
        return left;
    }
    
    friend bigInteger &&operator<<(bigInteger &&left, const bigInteger &right) {
        if (right.arr.empty() || left.arr.empty())
            return static_cast<bigInteger &&>(left);
        if (right.negative)
            throw std::logic_error("left shift count is a negative number");
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        std::uint32_t times1 = right.arr.back() & 31;
        for (auto iter = left.arr.rbegin(); iter != left.arr.rend(); *iter++ = static_cast<std::uint32_t>(u64 = static_cast<std::uint64_t>(*iter) << times1 | u32[1]));
        if (0 != u32[1])
            left.arr.push_front(u32[1]);
            for (bigInteger times32(right); !times32.absLessThan(32); times32.absSub(32), left.arr.push_back(0));
        return static_cast<bigInteger &&>(left);
    }

    friend bigInteger &&operator<<(bigInteger &&left, bigInteger &&right) {
        if (right.arr.empty() || left.arr.empty())
            return static_cast<bigInteger &&>(left);
        if (right.negative)
            throw std::logic_error("left shift count is a negative number");
        union {
            std::uint64_t u64 = 0;
            std::uint32_t u32[2];
        };
        std::uint32_t times1 = right.arr.back() & 31;
        for (auto iter = left.arr.rbegin(); iter != left.arr.rend(); *iter++ = static_cast<std::uint32_t>(u64 = static_cast<std::uint64_t>(*iter) << times1 | u32[1]));
        if (0 != u32[1])
            left.arr.push_front(u32[1]);
        for (; !right.absLessThan(32); right.absSub(32), left.arr.push_back(0));
        return static_cast<bigInteger &&>(left);
    }
    
    friend bigInteger operator<<(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) << right;
    }
    
    friend bigInteger operator<<(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) << std::move(right);
    }

    friend bigInteger &operator&=(bigInteger &left, const bigInteger &right) {
        bigInteger temp(right);
        if (left.negative) {
            left.absSub(1);
            if (temp.negative) {
                temp.absSub(1);
                left.absBitOrSelf(temp);
                left.absAdd(1);
            } else {
                temp.absBitAndSelf(left, 0xffffffff);
                left = std::move(temp);
            }
        } else if (temp.negative) {
            temp.absSub(1);
            left.absBitAndSelf(temp, 0xffffffff);
        } else
            left.absBitAndSelf(temp);
        return left;
    }

    friend bigInteger &operator&=(bigInteger &left, bigInteger &&right) {
        if (left.negative) {
            left.absSub(1);
            if (right.negative) {
                right.absSub(1);
                left.absBitOrSelf(right);
                left.absAdd(1);
            } else {
                right.absBitAndSelf(left, 0xffffffff);
                left = std::move(right);
            }
        } else if (right.negative) {
            right.absSub(1);
            left.absBitAndSelf(right, 0xffffffff);
        } else
            left.absBitAndSelf(right);
        return left;
    }
    
    friend bigInteger &&operator&(bigInteger &&left, const bigInteger &right) {
        bigInteger temp(right);
        if (left.negative) {
            left.absSub(1);
            if (temp.negative) {
                temp.absSub(1);
                left.absBitOrSelf(temp);
                left.absAdd(1);
            } else {
                temp.absBitAndSelf(left, 0xffffffff);
                left = std::move(temp);
            }
        } else if (temp.negative) {
            temp.absSub(1);
            left.absBitAndSelf(temp, 0xffffffff);
        } else
            left.absBitAndSelf(temp);
            return static_cast<bigInteger &&>(left);
    }

    friend bigInteger &&operator&(bigInteger &&left, bigInteger &&right) {
        if (left.negative) {
            left.absSub(1);
            if (right.negative) {
                right.absSub(1);
                left.absBitOrSelf(right);
                left.absAdd(1);
            } else {
                right.absBitAndSelf(left, 0xffffffff);
                left = std::move(right);
            }
        } else if (right.negative) {
            right.absSub(1);
            left.absBitAndSelf(right, 0xffffffff);
        } else
            left.absBitAndSelf(right);
        return static_cast<bigInteger &&>(left);
    }
    
    friend bigInteger operator&(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) & right;
    }
    
    friend bigInteger operator&(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) & std::move(right);
    }

    friend bigInteger &operator|=(bigInteger &left, const bigInteger &right) {
        bigInteger temp(right);
        if (left.negative) {
            left.absSub(1);
            if (temp.negative) {
                temp.absSub(1);
                left.absBitAndSelf(temp);
                left.negative = true;
                left.absAdd(1);
            } else {
                left.absBitAndSelf(temp, 0xffffffff);
                left.negative = true;
                left.absAdd(1);
            }
        } else if (temp.negative) {
            temp.absSub(1);
            temp.absBitAndSelf(left, 0xffffffff);
            left = std::move(temp);
            left.negative = true;
            left.absAdd(1);
        } else
            left.absBitOrSelf(temp);
        return left;
    }

    friend bigInteger &operator|=(bigInteger &left, bigInteger &&right) {
        if (left.negative) {
            left.absSub(1);
            if (right.negative) {
                right.absSub(1);
                left.absBitAndSelf(right);
                left.negative = true;
                left.absAdd(1);
            } else {
                left.absBitAndSelf(right, 0xffffffff);
                left.negative = true;
                left.absAdd(1);
            }
        } else if (right.negative) {
            right.absSub(1);
            right.absBitAndSelf(left, 0xffffffff);
            left = std::move(right);
            left.negative = true;
            left.absAdd(1);
        } else
            left.absBitOrSelf(right);
        return left;
    }
    
    friend bigInteger &&operator|(bigInteger &&left, const bigInteger &right) {
        bigInteger temp(right);
        if (left.negative) {
            left.absSub(1);
            if (temp.negative) {
                temp.absSub(1);
                left.absBitAndSelf(temp);
                left.negative = true;
                left.absAdd(1);
            } else {
                left.absBitAndSelf(temp, 0xffffffff);
                left.negative = true;
                left.absAdd(1);
            }
        } else if (temp.negative) {
            temp.absSub(1);
            temp.absBitAndSelf(left, 0xffffffff);
            left = std::move(temp);
            left.negative = true;
            left.absAdd(1);
        } else
            left.absBitOrSelf(temp);
        return static_cast<bigInteger &&>(left);
    }

    friend bigInteger &&operator|(bigInteger &&left, bigInteger &&right) {
        if (left.negative) {
            left.absSub(1);
            if (right.negative) {
                right.absSub(1);
                left.absBitAndSelf(right);
                left.negative = true;
                left.absAdd(1);
            } else {
                left.absBitAndSelf(right, 0xffffffff);
                left.negative = true;
                left.absAdd(1);
            }
        } else if (right.negative) {
            right.absSub(1);
            right.absBitAndSelf(left, 0xffffffff);
            left = std::move(right);
            left.negative = true;
            left.absAdd(1);
        } else
            left.absBitOrSelf(right);
        return static_cast<bigInteger &&>(left);
    }
    
    friend bigInteger operator|(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) | right;
    }
    
    friend bigInteger operator|(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) | std::move(right);
    }

    friend bigInteger &operator^=(bigInteger &left, const bigInteger &right) {
        bigInteger o(right);
        bool neg = left.negative ^ o.negative;
        if (left.negative)
            left.absSub(1);
        if (o.negative)
            o.absSub(1);
        left.negative ^= neg;
        auto oiter = o.arr.crbegin();
        for (auto iter = left.arr.rbegin(); left.arr.rend() != iter && o.arr.crend() != oiter; *iter++ ^= *oiter++);
        for (; o.arr.crend() != oiter; left.arr.push_back(*oiter++));
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        if (left.negative)
            left.absAdd(1);
        return left;
    }

    friend bigInteger &operator^=(bigInteger &left, bigInteger &&right) {
        bool neg = left.negative ^ right.negative;
        if (left.negative)
            left.absSub(1);
        if (right.negative)
            right.absSub(1);
        left.negative = neg;
        auto oiter = right.arr.rbegin();
        for (auto iter = left.arr.rbegin(); iter != left.arr.rend() && oiter != right.arr.rend(); *iter++ ^= *oiter++);
        for (; oiter != right.arr.rend(); left.arr.push_back(*oiter++));
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        if (left.negative)
            left.absAdd(1);
        return left;
    }
    
    friend bigInteger &&operator^(bigInteger &&left, const bigInteger &right) {
        bigInteger o(right);
        bool neg = left.negative ^ o.negative;
        if (left.negative)
            left.absSub(1);
        if (o.negative)
            o.absSub(1);
        left.negative = neg;
        auto oiter = o.arr.rbegin();
        for (auto iter = left.arr.rbegin(); iter != left.arr.rend() && oiter != o.arr.rend(); *iter++ ^= *oiter++);
        for (; oiter != o.arr.rend(); left.arr.push_back(*oiter++));
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        if (left.negative)
            left.absAdd(1);
        return static_cast<bigInteger &&>(left);
    }

    friend bigInteger &&operator^(bigInteger &&left, bigInteger &&right) {
        bool neg = left.negative ^ right.negative;
        if (left.negative)
            left.absSub(1);
        if (right.negative)
            right.absSub(1);
        left.negative = neg;
        auto oiter = right.arr.rbegin();
        for (auto iter = left.arr.rbegin(); iter != left.arr.rend() && oiter != right.arr.rend(); *iter++ ^= *oiter++);
        for (; oiter != right.arr.rend(); left.arr.push_back(*oiter++));
        while (!left.arr.empty() && 0 == left.arr.front())
            left.arr.pop_front();
        if (left.negative)
            left.absAdd(1);
        return static_cast<bigInteger &&>(left);
    }
    
    friend bigInteger operator^(const bigInteger &left, const bigInteger &right) {
        return bigInteger(left) ^ right;
    }
    
    friend bigInteger operator^(const bigInteger &left, bigInteger &&right) {
        return bigInteger(left) ^ std::move(right);
    }

    friend bigInteger operator~(const bigInteger &o) {
        bigInteger ret(o);
        if (ret.negative)
            ret.absSub(1);
        else
            ret.absAdd(1);
        ret.negative ^= true;
        return ret;
    }

    friend bigInteger &&operator~(bigInteger &&o) {
        if (o.negative)
            o.absSub(1);
        else
            o.absAdd(1);
        o.negative ^= true;
        return static_cast<bigInteger &&>(o);
    }
    
    friend bigInteger &operator++(bigInteger &o) {
        return o += 1;
    }
    
    friend bigInteger &operator--(bigInteger &o) {
        return o -= 1;
    }
    
    friend bigInteger &&operator++(bigInteger &&o) {
        return static_cast<bigInteger &&>(o += 1);
    }
    
    friend bigInteger &&operator--(bigInteger &&o) {
        return static_cast<bigInteger &&>(o -= 1);
    }

    friend bigInteger operator++(bigInteger &o, int) {
        bigInteger ret(o);
        o += 1;
        return ret;
    }
    
    friend bigInteger operator--(bigInteger &o, int) {
        bigInteger ret(o);
        o -= 1;
        return ret;
    }
    
    bigInteger gcd() const {
        return abs();
    }
    
    bigInteger gcd(const bigInteger &o) const {
        if (o.isZero())
            return abs();
        auto temp1(this->abs()), temp2(o.abs());
        while (temp1 %= temp2)
            std::swap(temp1, temp2);
        return temp2;
    }
    
    template<typename... Tp>
    bigInteger gcd(const bigInteger &o, const Tp &... args) const {
        if (o.isZero())
            return abs().gcd(args...);
        auto temp1(this->abs()), temp2(o.abs());
        while (temp1 %= temp2)
            std::swap(temp1, temp2);
        return temp2.gcd(args...);
    }
    
    bigInteger lcm() const {
        return abs();
    }
    
    bigInteger lcm(const bigInteger &o) const {
        if (isZero() || o.isZero())
            return 0;
        auto temp1(this->abs()), temp2(o.abs());
        return temp1 * temp2 / temp1.gcd(temp2);
    }
    
    template<typename... Tp>
    bigInteger lcm(const bigInteger &o, const Tp &... args) const {
        if (isZero() || o.isZero())
            return 0;
        auto temp1(this->abs()), temp2(o.abs());
        return (temp1 * temp2 / temp1.gcd(temp2)).lcm(args...);
    }
    
    template<typename... Tp>
    static bigInteger GCD(const bigInteger &o, const Tp &... args) {
        return o.gcd(args...);
    }
    
    template<typename... Tp>
    static bigInteger LCM(const bigInteger &o, const Tp &... args) {
        return o.lcm(args...);
    }

    friend bool operator==(const bigInteger &left, const bigInteger &right) {
        return left.absEquals(right) && left.negative == right.negative;
    }

    friend bool operator!=(const bigInteger &left, const bigInteger &right) {
        return left.negative != right.negative || !left.absEquals(right);
    }

    friend bool operator<=(const bigInteger &left, const bigInteger &right) {
        return left.negative != right.negative ? left.negative : (left.negative ? !left.absLessThan(right) : !left.absGreaterThan(right));
    }

    friend bool operator<(const bigInteger &left, const bigInteger &right) {
        return left.negative != right.negative ? left.negative : (left.negative ? left.absGreaterThan(right) : left.absLessThan(right));
    }

    friend bool operator>=(const bigInteger &left, const bigInteger &right) {
        return left.negative != right.negative ? right.negative : (left.negative ? !left.absGreaterThan(right) : !left.absLessThan(right));
    }

    friend bool operator>(const bigInteger &left, const bigInteger &right) {
        return left.negative != right.negative ? right.negative : (left.negative ? left.absLessThan(right) : left.absGreaterThan(right));
    }

    friend std::strong_ordering operator<=>(const bigInteger &left, const bigInteger &right) {
        if (left.negative != right.negative)
            return left.negative ? std::strong_ordering::less : std::strong_ordering::greater;
        if (left.absEquals(right))
            return std::strong_ordering::equal;
        return left.absLessThan(right) != left.negative ? std::strong_ordering::less : std::strong_ordering::greater;
    }

    explicit operator bool() const noexcept {
        return !isZero();
    }
    
    explicit operator std::list<bool>() const {
        std::list<bool> ret;
        for (auto i : arr)
            for (std::uint32_t j = 0x80000000; j > 0; j >>= 1)
                ret.push_back(0 != (j & i));
        while (!ret.empty() && !ret.front())
            ret.pop_front();
        return ret;
    }
    
    template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
    explicit operator T() const {
        T ret = 0;
        std::size_t digits = std::numeric_limits<T>::digits;
        if constexpr (std::numeric_limits<T>::is_signed)
            ++digits;
        if (arr.size() < (std::numeric_limits<T>::digits + 31) / 32) {
            for (uint32_t item : arr) {
                if constexpr (std::numeric_limits<T>::digits >= 32) {
                    ret <<= 31;
                    ret <<= 1;
                }
                ret |= item;
            }
        } else {
            auto iter = arr.end();
            for (std::size_t i = 0; i < (std::numeric_limits<T>::digits + 31) / 32; ++i)
                --iter;
            for (; iter != arr.end(); ++iter) {
                if constexpr (std::numeric_limits<T>::digits >= 32) {
                    ret <<= 31;
                    ret <<= 1;
                }
                ret |= *iter;
            }
        }
        if (negative)
            ret = -ret;
        return ret;
    }

    constexpr inline bool isZero() const noexcept {
        return arr.empty();
    }

    constexpr inline bool isNegative() const noexcept {
        return negative;
    }

    constexpr inline bool isPositive() const noexcept {
        return !arr.empty() && !negative;
    }

    bigInteger abs() const & {
        bigInteger ret(*this);
        ret.negative = false;
        return ret;
    }

    bigInteger &&abs() && {
        negative = false;
        return static_cast<bigInteger &&>(*this);
    }

    bigInteger pow(const bigInteger &o) const {
        bigInteger temp(1);
        if (o.isNegative() || isZero())
            return 0;
        for (bool b:(std::list<bool>) o) {
            temp *= temp;
            if (b)
                temp *= (*this);
        }
        return temp;
    }

    bigInteger powMod(const bigInteger &o, const bigInteger &m) const {
        bigInteger temp(1), temp2(*this % m);
        if (o.isNegative() || temp2.isZero() || m.isZero()) {
            temp.arr.clear();
            return temp;
        }
        for (bool b:(std::list<bool>) o) {
            temp *= temp;
            temp %= m;
            if (b) {
                temp *= temp2;
                temp %= m;
            }
        }
        return temp;
    }

    template<typename charT = char, typename Traits = std::char_traits<charT>, typename Allocator = std::allocator<charT>>
    std::basic_string<charT, Traits, Allocator> getString() const {
        if (isZero()) {
            auto cc = new charT[2];
            cc[0] = (charT) '0';
            std::basic_string<charT, Traits, Allocator> ret(cc);
            delete[] cc;
            return ret;
        }
        std::list<charT> ll;
        for (bigInteger i(this->abs()); !i.isZero(); i /= bigInteger(10))
            ll.push_front((charT) (int) ((i % 10) + '0'));
        if (negative)
            ll.push_front('-');
        auto cc = new charT[ll.size() + 1];
        cc[ll.size()] = '\0';
        for (std::size_t i = 0; !ll.empty(); ll.pop_front())
            cc[i++] = ll.front();
        std::basic_string<charT, Traits, Allocator> ret(cc);
        delete[] cc;
        return ret;
    }

    template<typename charT = char, typename Traits = std::char_traits<charT>, typename Allocator = std::allocator<charT>>
    std::basic_string<charT, Traits, Allocator> getOctString() const {
        if (isZero()) {
            auto cc = new charT[2];
            cc[0] = (charT) '0';
            std::basic_string<charT, Traits, Allocator> ret(cc);
            delete[] cc;
            return ret;
        }
        auto bits = (std::list<bool>) *this;
        while (0 != bits.size() % 3) bits.push_front(false);
        std::list<charT> ll;
        while (!bits.empty()) {
            int temp = 0;
            if (bits.front()) temp += 4;
            bits.pop_front();
            if (bits.front()) temp += 2;
            bits.pop_front();
            if (bits.front()) temp += 1;
            bits.pop_front();
            ll.push_back(temp + '0');
        }
        ll.push_front('0');
        if (negative)
            ll.push_front('-');
        auto *cc = new charT[ll.size() + 1];
        cc[ll.size()] = '\0';
        for (std::size_t i = 0; !ll.empty(); ll.pop_front())
            cc[i++] = ll.front();
        std::basic_string<charT, Traits, Allocator> ret(cc);
        delete[] cc;
        return ret;
    }

    template<typename charT = char, typename Traits = std::char_traits<charT>, typename Allocator = std::allocator<charT>>
    std::basic_string<charT, Traits, Allocator> getHexLowerString() const {
        if (isZero()) {
            auto cc = new charT[2];
            cc[0] = (charT) '0';
            std::basic_string<charT, Traits, Allocator> ret(cc);
            delete[] cc;
            return ret;
        }
        auto bits = (std::list<bool>) *this;
        while (0 != (bits.size() & 3)) bits.push_front(false);
        std::list<charT> ll;
        while (!bits.empty()) {
            int temp = 0;
            if (bits.front()) temp += 8;
            bits.pop_front();
            if (bits.front()) temp += 4;
            bits.pop_front();
            if (bits.front()) temp += 2;
            bits.pop_front();
            if (bits.front()) temp += 1;
            bits.pop_front();
            if (temp <= 9) ll.push_back(temp + '0'); else ll.push_back(temp - 10 + 'a');
        }
        ll.push_front('x');
        ll.push_front('0');
        if (negative)
            ll.push_front('-');
        auto *cc = new charT[ll.size() + 1];
        cc[ll.size()] = '\0';
        for (std::size_t i = 0; !ll.empty(); ll.pop_front())
            cc[i++] = ll.front();
        std::basic_string<charT, Traits, Allocator> ret(cc);
        delete[] cc;
        return ret;
    }
            
    template<typename charT = char, typename Traits = std::char_traits<charT>, typename Allocator = std::allocator<charT>>
    std::basic_string<charT, Traits, Allocator> getHexUpperCaseString() const {
        if (isZero()) {
            auto cc = new charT[2];
            cc[0] = (charT) '0';
            std::basic_string<charT, Traits, Allocator> ret(cc);
            delete[] cc;
            return ret;
        }
        auto bits = (std::list<bool>) *this;
        while (0 != (bits.size() & 3)) bits.push_front(false);
        std::list<charT> ll;
        while (!bits.empty()) {
            int temp = 0;
            if (bits.front()) temp += 8;
            bits.pop_front();
            if (bits.front()) temp += 4;
            bits.pop_front();
            if (bits.front()) temp += 2;
            bits.pop_front();
            if (bits.front()) temp += 1;
            bits.pop_front();
            if (temp <= 9) ll.push_back(temp + '0'); else ll.push_back(temp - 10 + 'A');
        }
        ll.push_front('X');
        ll.push_front('0');
        if (negative)
            ll.push_front('-');
        auto *cc = new charT[ll.size() + 1];
        cc[ll.size()] = '\0';
        for (std::size_t i = 0; !ll.empty(); ll.pop_front())
            cc[i++] = ll.front();
        std::basic_string<charT, Traits, Allocator> ret(cc);
        delete[] cc;
        return ret;
    }

    template<typename charT = char, typename Traits = std::char_traits<charT>>
    friend std::basic_istream<charT, Traits> &operator>>(std::basic_istream<charT, Traits> &in, bigInteger &bi) {
        if (!in.good())
            return in;
        bool isNegative = false;
        bi.negative = false;
        bi.arr.clear();
        while (' ' == in.peek() || '\t' == in.peek()) in.get();
        if ('-' == in.peek()) {
            in.get();
            isNegative = true;
        } else if ('+' == in.peek()) {
            in.get();
        }
        if (std::ios_base::dec != (in.flags() & std::ios_base::basefield)) {
            in.get();
            std::list<bool> l;
            if (std::ios_base::hex == (in.flags() & std::ios_base::basefield)) {
                in.get();
                while ('_' == in.peek() || ('0' <= in.peek() && '9' >= in.peek()) || ('a' <= in.peek() && 'f' >= in.peek()) || ('A' <= in.peek() && 'F' >= in.peek())) {
                    switch (in.peek()) {
                        case '0':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '1':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case '2':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case '3':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        case '4':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '5':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case '6':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case '7':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        case '8':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '9':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case 'a':
                        case 'A':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case 'b':
                        case 'B':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        case 'c':
                        case 'C':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case 'd':
                        case 'D':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case 'e':
                        case 'E':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case 'f':
                        case 'F':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        default:
                            break;
                    }
                    in.get();
                }
            } else {
                while ('_' == in.peek() || ('0' <= in.peek() && '7' >= in.peek())) {
                    switch (in.peek()) {
                        case '0':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '1':
                            l.push_back(false);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case '2':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case '3':
                            l.push_back(false);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        case '4':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(false);
                            break;
                        case '5':
                            l.push_back(true);
                            l.push_back(false);
                            l.push_back(true);
                            break;
                        case '6':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(false);
                            break;
                        case '7':
                            l.push_back(true);
                            l.push_back(true);
                            l.push_back(true);
                            break;
                        default:
                            break;
                    }
                    in.get();
                }
            }
            while (!l.empty() && false == l.front()) l.pop_front();
            while (0 != (l.size() & 31)) l.push_front(false);
            while (!l.empty()) {
                std::uint32_t temp = 0;
                for (int i = 31; i >= 0; --i) {
                    if (l.front())
                        temp |= 1 << i;
                    l.pop_front();
                }
                bi.arr.push_back(temp);
            }
        } else {
            while ('_' == in.peek() || ('0' <= in.peek() && '9' >= in.peek())) {
                if ('_' == in.peek())
                    continue;
                bi *= 10;
                bi += (int) in.peek() - '0';
                in.get();
            }
        }
        bi.negative = isNegative && !bi.isZero();
        return in;
    }

    template<typename charT = char, typename Traits = std::char_traits<charT>>
    friend std::basic_ostream<charT, Traits> &operator<<(std::basic_ostream<charT, Traits> &os, const bigInteger &bi) {
        switch(os.flags() & std::ios_base::basefield) {
            case std::ios_base::hex:
                if (0 == (os.flags() & std::ios_base::uppercase))
                    os << bi.getHexLowerString<charT>();
                else
                    os << bi.getHexUpperCaseString<charT>();
                break;
            case std::ios_base::oct:
                os << bi.getOctString<charT>();
                break;
            case std::ios_base::dec:
                os << bi.getString<charT>();
                break;
            default:
                std::cerr << "unsupported base when print a bigInteger" << std::endl;
        }
        return os;
    }
            
    static bigInteger rand(const bigInteger &o) {
        static std::default_random_engine e;
        e.seed(static_cast<unsigned>(clock()));
        if(!o.isPositive())
            return bigInteger(0);
        bigInteger ret(o.arr.front());
        while(ret.arr.size() < o.arr.size())
            ret.arr.push_back(e());
        if(ret.absGreaterThan(o)){
            std::uniform_int_distribution<std::uint32_t> u(0, o.arr.front() - 1);
            ret.arr.front() = u(e);
        }else{
            std::uniform_int_distribution<std::uint32_t> u(0,o.arr.front());
            ret.arr.front() = u(e);
        }
        while (!ret.isZero() && 0 == ret.arr.front())
            ret.arr.pop_front();
        return ret;
    }
            
    static bigInteger rand(const bigInteger &o, const bigInteger &i) {
        return o > i ? i + rand(o - i) : o + rand(i - o);
    }

};

bigInteger operator"" _b(const char * const str) {
    return bigInteger(std::string(str));
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator+=(T &left, const bigInteger &right) {
    return left += (T)right;
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator-=(T &left, const bigInteger &right) {
    return left -= (T)right;
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator*=(T &left, const bigInteger &right) {
    return left *= (T)right;
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator/=(T &left, const bigInteger &right) {
    if (left < right)
        return left = 0;
    return left /= (T)right;
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator%=(T &left, const bigInteger &right) {
    return left < right ? left : left %= (T)right;
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator<<=(T &left, const bigInteger &right) {
    return left <<= (T)right;
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator>>=(T &left, const bigInteger &right) {
    return left >>= (T)right;
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator&=(T &left, const bigInteger &right) {
    return left &= (T)right;
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator|=(T &left, const bigInteger &right) {
    return left |= (T)right;
}

template<typename T, std::enable_if_t<std::is_integral_v<T>, T*> ptr = nullptr>
T &operator^=(T &left, const bigInteger &right) {
    return left ^= (T)right;
}

template<>
struct std::hash<bigInteger> : public unary_function<bigInteger, std::size_t> {
    std::size_t operator()(const bigInteger &val) const noexcept {
        return std::hash<std::string>()(val.getHexLowerString());
    };
};

#endif //BIGINTEGER_H
