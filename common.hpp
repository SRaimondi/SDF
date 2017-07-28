//
// Created by Simone Raimondi on 22.07.17.
//

#ifndef SDF_COMMON_HPP
#define SDF_COMMON_HPP

#include <cstdlib>

// Return number if it's positive, else return zero
template<typename T>
static inline T Positive(T x) {
    return x > (T) 0 ? x : (T) 0;
}

// Return squared number if positive, else return zero
template<typename T>
static inline T PositiveSQ(T x) {
    return x > (T) 0 ? x * x : (T) 0;
}

// Return number if it's negative, else return zero
template<typename T>
static inline T Negative(T x) {
    return x < (T) 0 ? x : (T) 0;
}

// Return squared number if negative, else return zero
template<typename T>
static inline T NegativeSQ(T x) {
    return x < (T) 0 ? x * x : (T) 0;
}

// Return maximum between two numbers
template<typename T>
static inline T Max(T x, T y) {
    return x > y ? x : y;
}

// Return minimum between two numbers
template<typename T>
static inline T Min(T x, T y) {
    return x < y ? x : y;
}

// Lerp
template<typename T>
static inline T Lerp(T t, T a, T b) {
    return ((T) 1 - t) * a + t * b;
}

#endif //SDF_COMMON_HPP
