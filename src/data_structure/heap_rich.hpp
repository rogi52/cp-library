#include "../../src/cp-template.hpp"

namespace std_ex {

template < class T > T pop(std::queue < T >& c) { T x = c.front(); c.pop(); return x; }
template < class T > T pop(std::stack < T >& c) { T x = c.top();   c.pop(); return x; }
template < class T > T pop(std::vector< T >& c) { T x = c.top();   c.pop(); return x; } 
template < class T > T pop(   heap_min< T >& c) { T x = c.top();   c.pop(); return x; }
template < class T > T pop(   heap_max< T >& c) { T x = c.top();   c.pop(); return x; }

}

template < class container >
struct erasable {
    using T = typename container::value_type;
    container c, rm_c;

    erasable() {}
    erasable(std::vector< T >& a) : c(a.begin(), a.end()) {}
    void push(T x) { c.push(x); }
    int size() const { return c.size() - rm_c.size(); }
    int empty() { return size() == 0; }
    T pop() { set(); return std_ex::pop(c); }
    T top() { set(); return c.top(); }
    void erase(T x) { rm_c.push(x); }

  private:
    void set() { while(not rm_c.empty() and rm_c.top() == c.top()) rm_c.pop(), c.pop(); }
};


template < class T >
struct heap_minmax {
    erasable< heap_min< T > > h_min;
    erasable< heap_max< T > > h_max;

    heap_minmax() {}
    heap_minmax(std::vector< T >& a) : h_min(a), h_max(a) {}
    void push(T x) {
        h_min.push(x), h_max.push(x);
    }
    int size() { return h_min.size(); }
    int empty() { return size() == 0; }
    T min() { return h_min.top(); }
    T max() { return h_max.top(); }
    T pop_min() {
        T x = h_min.top(); h_min.pop();
        h_max.erase(x);
        return x;
    }
    T pop_max() {
        T x = h_max.top(); h_max.pop();
        h_min.erase(x);
        return x;
    }
    void erase(T x) {
        h_min.erase(x);
        h_max.erase(x);
    }
};
