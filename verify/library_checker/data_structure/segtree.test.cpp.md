---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: src/algebra/affine.hpp
    title: src/algebra/affine.hpp
  - icon: ':question:'
    path: src/cp-template.hpp
    title: src/cp-template.hpp
  - icon: ':heavy_check_mark:'
    path: src/data_structure/segtree.hpp
    title: src/data_structure/segtree.hpp
  - icon: ':heavy_check_mark:'
    path: src/number/modint.hpp
    title: modint
  - icon: ':question:'
    path: src/utility/io.hpp
    title: src/utility/io.hpp
  - icon: ':question:'
    path: src/utility/rep_itr.hpp
    title: src/utility/rep_itr.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/point_set_range_composite
    links:
    - https://judge.yosupo.jp/problem/point_set_range_composite
  bundledCode: "#line 1 \"verify/library_checker/data_structure/segtree.test.cpp\"\
    \n#define PROBLEM \"https://judge.yosupo.jp/problem/point_set_range_composite\"\
    \n\n#line 2 \"src/cp-template.hpp\"\n#include <bits/stdc++.h>\nusing namespace\
    \ std;\nusing ll = long long;\nusing ld = long double;\nusing uint = unsigned\
    \ int;\nusing ull  = unsigned long long;\nusing i128 = __int128_t;\ntemplate <\
    \ class T > bool chmin(T& a, T b) { if(a > b) { a = b; return true; } return false;\
    \ }\ntemplate < class T > bool chmax(T& a, T b) { if(a < b) { a = b; return true;\
    \ } return false; }\n\n#line 2 \"src/utility/rep_itr.hpp\"\ntemplate < class T\
    \ > struct itr {\n    T i, d;\n    constexpr itr(const T i) noexcept : i(i), d(1)\
    \ {}\n    constexpr itr(const T i, const T d) noexcept : i(i), d(d) {}\n    void\
    \ operator++() noexcept { i += d; }\n    constexpr int operator*() const noexcept\
    \ { return i; }\n    constexpr bool operator!=(const itr x) const noexcept {\n\
    \        return d > 0 ? i < x.i : i > x.i;\n    }\n};\n\ntemplate < class T >\
    \ struct rep {\n    const itr< T > s, t;\n    constexpr rep(const T t) noexcept\
    \ : s(0), t(t) {}\n    constexpr rep(const T s, const T t) noexcept : s(s), t(t)\
    \ {}\n    constexpr rep(const T s, const T t, const T d) noexcept : s(s, d), t(t,\
    \ d) {}\n    constexpr auto begin() const noexcept { return s; }\n    constexpr\
    \ auto end() const noexcept { return t; }\n};\n\ntemplate < class T > struct revrep\
    \ {\n    const itr < T > s, t;\n    constexpr revrep(const T t) noexcept : s(t\
    \ - 1, -1), t(-1, -1) {}\n    constexpr revrep(const T s, const T t) noexcept\
    \ : s(t - 1, -1), t(s - 1, -1) {}\n    constexpr revrep(const T s, const T t,\
    \ const T d) noexcept : s(t - 1, -d), t(s - 1, -d) {}\n    constexpr auto begin()\
    \ const noexcept { return s; }\n    constexpr auto end() const noexcept { return\
    \ t; }\n};\n#line 2 \"src/utility/io.hpp\"\nnamespace scanner {\n    struct sca\
    \ {\n        template < class T > operator T() {\n            T s; cin >> s; return\
    \ s;\n        }\n    };\n    struct vec {\n        int n;\n        vec(int n)\
    \ : n(n) {}\n        template < class T > operator vector< T >() {\n         \
    \   vector< T > v(n);\n            for(T& x : v) cin >> x;\n            return\
    \ v;\n        }\n    };\n    struct mat {\n        int h,w;\n        mat(int h,\
    \ int w) : h(h), w(w) {}\n        template < class T > operator vector< vector<\
    \ T > >() {\n            vector m(h, vector< T >(w));\n            for(vector<\
    \ T >& v : m) for(T& x : v) cin >> x;\n            return m;\n        }\n    };\n\
    \    struct speedup {\n        speedup() {\n            cin.tie(0);\n        \
    \    ios::sync_with_stdio(0);\n        }\n    } su;\n}\nscanner::sca in() { return\
    \ scanner::sca(); }\nscanner::vec in(int n) { return scanner::vec(n); }\nscanner::mat\
    \ in(int h, int w) { return scanner::mat(h, w); }\n\nnamespace printer {\n   \
    \ void precision(int d) {\n        cout << fixed << setprecision(d);\n    }\n\
    \    void flush() {\n        cout.flush();\n    }\n}\nint print() { cout << '\\\
    n'; return 0; }\ntemplate < class head, class... tail > int print(head&& h, tail&&...\
    \ t) {\n    cout << h; if(sizeof...(tail)) cout << ' ';\n    return print(forward<tail>(t)...);\n\
    }\ntemplate < class T > int print(vector< T > a, char sep = ' ') {\n    int n\
    \ = a.size();\n    for(int i : rep(n)) cout << a[i] << (i != n - 1 ? sep : '\\\
    n');\n    return 0;\n}\ntemplate < class T > int print(vector< vector< T > > a)\
    \ {\n    if(a.empty()) return 0;\n    int h = a.size(), w = a[0].size();\n   \
    \ for(int i : rep(h)) for(int j : rep(w)) cout << a[i][j] << (j != w - 1 ? ' '\
    \ : '\\n');\n    return 0;\n}\n#line 1 \"src/data_structure/segtree.hpp\"\ntemplate\
    \ < class monoid > struct segtree {\n    using S = typename monoid::set;\n\n \
    \   segtree() : segtree(0) {}\n    segtree(int n) : segtree(vector< S >(n, monoid::id))\
    \ {}\n    segtree(const vector< S >& v) : _n(int(v.size())) {\n        log = ceil_pow2(_n);\n\
    \        size = 1 << log;\n        d = vector< S >(2 * size, monoid::id);\n  \
    \      for(int i = 0; i < _n; i++) d[size + i] = v[i];\n        for(int i = size\
    \ - 1; i >= 1; i--) update(i);\n    }\n    // a[i] <- x\n    void set(int i, S\
    \ x) {\n        assert(0 <= i && i < _n);\n        i += size;\n        d[i] =\
    \ x;\n        for(int p = 1; p <= log; p++) update(i >> p);\n    }\n    // a[i]\n\
    \    S get(int i) {\n        assert(0 <= i && i < _n);\n        return d[i + size];\n\
    \    }\n    // [l, r)\n    S prod(int l, int r) {\n        assert(0 <= l && l\
    \ <= r && r <= _n);\n        S sml = monoid::id, smr = monoid::id;\n        l\
    \ += size, r += size;\n        while(l < r) {\n            if(l & 1) sml = monoid::op(sml,\
    \ d[l++]);\n            if(r & 1) smr = monoid::op(d[--r], smr);\n           \
    \ l >>= 1, r >>= 1;\n        }\n        return monoid::op(sml, smr);\n    }\n\
    \    S all_prod() { return d[1]; }\n    template < class func > int max_right(int\
    \ l, func f) {\n        assert(0 <= l && l <= _n);\n        assert(f(monoid::id));\n\
    \        if(l == _n) return _n;\n        l += size;\n        S sm = monoid::id;\n\
    \        do {\n            while(l % 2 == 0) l >>= 1;\n            if(!f(monoid::op(sm,\
    \ d[l]))) {\n                while(l < size) {\n                    l = 2 * l;\n\
    \                    if(f(monoid::op(sm, d[l]))) {\n                        sm\
    \ = monoid::op(sm, d[l]);\n                        l++;\n                    }\n\
    \                }\n                return l - size;\n            }\n        \
    \    sm = monoid::op(sm, d[l]);\n            l++;\n        } while((l & -l) !=\
    \ l);\n        return _n;\n    }\n    template < class func > int min_left(int\
    \ r, func f) {\n        assert(0 <= r && r <= _n);\n        assert(f(monoid::id));\n\
    \        if(r == 0) return 0;\n        r += size;\n        S sm = monoid::id;\n\
    \        do {\n            r--;\n            while(r > 1 && (r % 2)) r >>= 1;\n\
    \            if(!f(monoid::op(d[r], sm))) {\n                while(r < size) {\n\
    \                    r = 2 * r + 1;\n                    if(f(monoid::op(d[r],\
    \ sm))) {\n                        sm = monoid::op(d[r], sm);\n              \
    \          r--;\n                    }\n                }\n                return\
    \ r + 1 - size;\n            }\n            sm = monoid::op(d[r], sm);\n     \
    \   } while((r & -r) != r);\n        return 0;\n    }\n\n  private:\n    int _n,\
    \ size, log;\n    vector< S > d;\n    int ceil_pow2(int n) { int x = 0; while((1U\
    \ << x) < uint(n)) x++; return x; }\n    void update(int k) { d[k] = monoid::op(d[2\
    \ * k], d[2 * k + 1]); }\n};\n#line 1 \"src/algebra/affine.hpp\"\ntemplate < class\
    \ T > class affine {\n  public:\n    T a, b;\n    constexpr affine() = default;\n\
    \    constexpr affine(const T &a, const T &b) : a(a), b(b) {}\n    constexpr T\
    \ eval(const T &x) const { return x * a + b; }\n    constexpr affine composite(const\
    \ affine &r) const {\n        return affine(a * r.a, b * r.a + r.b);\n    }\n\
    \    static constexpr affine id() {\n        return affine(T(1), T(0));\n    }\n\
    };\n\ntemplate < class T > class affine_composite_monoid {\n  public:\n    using\
    \ F = affine< T >;\n    using set = F;\n    static constexpr F op(const F &l,\
    \ const F &r) { return l.composite(r); }\n    static constexpr F id = F::id();\n\
    };\ntemplate < class T > constexpr affine< T > affine_composite_monoid< T >::id;\n\
    #line 1 \"src/number/modint.hpp\"\nstruct modinfo { uint mod, root, isprime; };\n\
    template < modinfo const &ref >\nstruct modint {\n    static constexpr uint const\
    \ &mod = ref.mod;\n    static constexpr uint const &root = ref.root;\n    static\
    \ constexpr uint const &isprime = ref.isprime;\n    uint v = 0;\n    constexpr\
    \ modint& s(uint v) { this->v = v < mod ? v : v - mod; return *this; }\n    constexpr\
    \ modint(ll v = 0) { s(v % mod + mod); }\n    modint operator-() const { return\
    \ modint() - *this; }\n    modint& operator+=(const modint& rhs) { return s(v\
    \ + rhs.v); }\n    modint& operator-=(const modint& rhs) { return s(v + mod -\
    \ rhs.v); }\n    modint& operator*=(const modint& rhs) { v = ull(v) * rhs.v %\
    \ mod; return *this; }\n    modint& operator/=(const modint& rhs) { return *this\
    \ *= inv(rhs); }\n    modint operator+(const modint& rhs) const { return modint(*this)\
    \ += rhs; }\n    modint operator-(const modint& rhs) const { return modint(*this)\
    \ -= rhs; }\n    modint operator*(const modint& rhs) const { return modint(*this)\
    \ *= rhs; }\n    modint operator/(const modint& rhs) const { return modint(*this)\
    \ /= rhs; }\n    friend modint pow(modint x, ll n) { modint res(1); while(n >\
    \ 0) { if(n & 1) res *= x; x *= x; n >>= 1; } return res; }\n    friend modint\
    \ inv(modint v) {\n        if(isprime) {\n            return pow(v, mod - 2);\n\
    \        } else {\n            ll a = v.v, b = modint::mod, x = 1, y = 0, t;\n\
    \            while(b > 0) { t = a / b; swap(a -= t * b, b); swap(x -= t * y, y);\
    \ }\n            return modint(x);\n        }\n    }\n    friend modint operator+(int\
    \ x, const modint& y) { return modint(x) + y; }\n    friend modint operator-(int\
    \ x, const modint& y) { return modint(x) - y; }\n    friend modint operator*(int\
    \ x, const modint& y) { return modint(x) * y; }\n    friend modint operator/(int\
    \ x, const modint& y) { return modint(x) / y; }\n    friend istream& operator>>(istream&\
    \ is, modint& m) { ll x; is >> x; m = modint(x); return is; }\n    friend ostream&\
    \ operator<<(ostream& os, const modint& m) { return os << m.v; }\n    bool operator==(const\
    \ modint& r) const { return v == r.v; }\n    bool operator!=(const modint& r)\
    \ const { return v != r.v; }\n    static uint get_mod() { return mod; }\n};\n\
    constexpr modinfo base998244353 { 998244353, 3, 1 };\nconstexpr modinfo base1000000007\
    \ { 1000000007, 0, 1 };\nusing mint998244353 = modint< base998244353 >;\nusing\
    \ mint1000000007 = modint< base1000000007 >;\n#line 7 \"verify/library_checker/data_structure/segtree.test.cpp\"\
    \n\nint main(){\n    cin.tie(0);\n    ios::sync_with_stdio(0);\n    \n    int\
    \ N,Q; cin >> N >> Q;\n    using mint = mint998244353;\n    using F = affine<mint>;\n\
    \    vector< F > f(N);\n    for(int i : rep(N)) {\n        mint a,b; cin >> a\
    \ >> b;\n        f[i] = F(a, b);\n    }\n\n    segtree< affine_composite_monoid<mint>\
    \ > seg(f);\n    for(int _ : rep(Q)) {\n        int t; cin >> t;\n        if(t\
    \ == 0) {\n            int p; mint c,d; cin >> p >> c >> d;\n            seg.set(p,\
    \ F(c, d));\n        }\n\n        if(t == 1) {\n            int l,r; mint x; cin\
    \ >> l >> r >> x;\n            cout << seg.prod(l, r).eval(x) << \"\\n\";\n  \
    \      }\n    }\n}\n"
  code: "#define PROBLEM \"https://judge.yosupo.jp/problem/point_set_range_composite\"\
    \n\n#include \"src/cp-template.hpp\"\n#include \"src/data_structure/segtree.hpp\"\
    \n#include \"src/algebra/affine.hpp\"\n#include \"src/number/modint.hpp\"\n\n\
    int main(){\n    cin.tie(0);\n    ios::sync_with_stdio(0);\n    \n    int N,Q;\
    \ cin >> N >> Q;\n    using mint = mint998244353;\n    using F = affine<mint>;\n\
    \    vector< F > f(N);\n    for(int i : rep(N)) {\n        mint a,b; cin >> a\
    \ >> b;\n        f[i] = F(a, b);\n    }\n\n    segtree< affine_composite_monoid<mint>\
    \ > seg(f);\n    for(int _ : rep(Q)) {\n        int t; cin >> t;\n        if(t\
    \ == 0) {\n            int p; mint c,d; cin >> p >> c >> d;\n            seg.set(p,\
    \ F(c, d));\n        }\n\n        if(t == 1) {\n            int l,r; mint x; cin\
    \ >> l >> r >> x;\n            cout << seg.prod(l, r).eval(x) << \"\\n\";\n  \
    \      }\n    }\n}\n"
  dependsOn:
  - src/cp-template.hpp
  - src/utility/rep_itr.hpp
  - src/utility/io.hpp
  - src/data_structure/segtree.hpp
  - src/algebra/affine.hpp
  - src/number/modint.hpp
  isVerificationFile: true
  path: verify/library_checker/data_structure/segtree.test.cpp
  requiredBy: []
  timestamp: '2023-05-10 11:30:59+09:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/library_checker/data_structure/segtree.test.cpp
layout: document
redirect_from:
- /verify/verify/library_checker/data_structure/segtree.test.cpp
- /verify/verify/library_checker/data_structure/segtree.test.cpp.html
title: verify/library_checker/data_structure/segtree.test.cpp
---