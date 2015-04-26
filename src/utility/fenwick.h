/**
 * Implementation of the Fenwick tree data structure,
 * supporting the following operations:
 *  a) get(x): gets the maximum element in the interval [0,x] 
 *  b) update(x, v): sets the x-th element to v, if v
 *                   is bigger than the previously stored value
 *                   on position x, otherwise nothing happens
 *  Note: x >= 0
 * @author: Filip Pavetic (fpavetic@gmail.com)
 */

#ifndef LCSKPP_FENWICK
#define LCSKPP_FENWICK

#include <cassert>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

template<class T>
class FenwickMax {
 public:
  FenwickMax(size_t n) {
    elements_ = std::vector<T> (n+1);
  }
  
  void update(size_t pos, const T& val) {
    ++pos;
    for ( ; pos < elements_.size(); pos += lobit(pos)) {
      elements_[pos] = std::max(elements_[pos], val);
    }
  }

  T get(size_t pos) {
    ++pos;
    T ret = T();
    for ( ; pos > 0; pos -= lobit(pos)) {
      ret = std::max(ret, elements_[pos]);
    }
    return ret;
  }

 private:
  size_t lobit(const size_t& a) { return a&-a; }
  
 private:
  std::vector<T> elements_;
};

#endif
