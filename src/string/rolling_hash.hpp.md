---
data:
  _extendedDependsOn:
  - icon: ':question:'
    path: src/algorithm/argsort.hpp
    title: src/algorithm/argsort.hpp
  - icon: ':question:'
    path: src/algorithm/bin_search.hpp
    title: src/algorithm/bin_search.hpp
  - icon: ':question:'
    path: src/cp-template.hpp
    title: src/cp-template.hpp
  - icon: ':heavy_check_mark:'
    path: src/utility/hash.hpp
    title: src/utility/hash.hpp
  - icon: ':question:'
    path: src/utility/heap.hpp
    title: src/utility/heap.hpp
  - icon: ':question:'
    path: src/utility/io.hpp
    title: src/utility/io.hpp
  - icon: ':question:'
    path: src/utility/key_val.hpp
    title: src/utility/key_val.hpp
  - icon: ':question:'
    path: src/utility/random.hpp
    title: src/utility/random.hpp
  - icon: ':question:'
    path: src/utility/rep_itr.hpp
    title: src/utility/rep_itr.hpp
  - icon: ':question:'
    path: src/utility/vec_op.hpp
    title: src/utility/vec_op.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/library_checker/string/rolling_hash.test.cpp
    title: verify/library_checker/string/rolling_hash.test.cpp
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 2 \"src/cp-template.hpp\"\n#include <bits/stdc++.h>\nusing namespace\
    \ std;\nusing ll = long long;\nusing ld = long double;\nusing uint = unsigned\
    \ int;\nusing ull  = unsigned long long;\nusing i32 = int;\nusing u32 = unsigned\
    \ int;\nusing i64 = long long;\nusing u64 = unsigned long long;\nusing i128 =\
    \ __int128_t;\ntemplate < class T > bool chmin(T& a, T b) { if(a > b) { a = b;\
    \ return true; } return false; }\ntemplate < class T > bool chmax(T& a, T b) {\
    \ if(a < b) { a = b; return true; } return false; }\n\n#line 2 \"src/utility/rep_itr.hpp\"\
    \ntemplate < class T > struct itr {\n    T i, d;\n    constexpr itr(const T i)\
    \ noexcept : i(i), d(1) {}\n    constexpr itr(const T i, const T d) noexcept :\
    \ i(i), d(d) {}\n    void operator++() noexcept { i += d; }\n    constexpr int\
    \ operator*() const noexcept { return i; }\n    constexpr bool operator!=(const\
    \ itr x) const noexcept {\n        return d > 0 ? i < x.i : i > x.i;\n    }\n\
    };\n\ntemplate < class T > struct rep {\n    const itr< T > s, t;\n    constexpr\
    \ rep(const T t) noexcept : s(0), t(t) {}\n    constexpr rep(const T s, const\
    \ T t) noexcept : s(s), t(t) {}\n    constexpr rep(const T s, const T t, const\
    \ T d) noexcept : s(s, d), t(t, d) {}\n    constexpr auto begin() const noexcept\
    \ { return s; }\n    constexpr auto end() const noexcept { return t; }\n};\n\n\
    template < class T > struct revrep {\n    const itr < T > s, t;\n    constexpr\
    \ revrep(const T t) noexcept : s(t - 1, -1), t(-1, -1) {}\n    constexpr revrep(const\
    \ T s, const T t) noexcept : s(t - 1, -1), t(s - 1, -1) {}\n    constexpr revrep(const\
    \ T s, const T t, const T d) noexcept : s(t - 1, -d), t(s - 1, -d) {}\n    constexpr\
    \ auto begin() const noexcept { return s; }\n    constexpr auto end() const noexcept\
    \ { return t; }\n};\n#line 2 \"src/utility/io.hpp\"\nnamespace scanner {\n   \
    \ struct sca {\n        template < class T > operator T() {\n            T s;\
    \ cin >> s; return s;\n        }\n    };\n    struct vec {\n        int n;\n \
    \       vec(int n) : n(n) {}\n        template < class T > operator vector< T\
    \ >() {\n            vector< T > v(n);\n            for(T& x : v) cin >> x;\n\
    \            return v;\n        }\n    };\n    struct mat {\n        int h,w;\n\
    \        mat(int h, int w) : h(h), w(w) {}\n        template < class T > operator\
    \ vector< vector< T > >() {\n            vector m(h, vector< T >(w));\n      \
    \      for(vector< T >& v : m) for(T& x : v) cin >> x;\n            return m;\n\
    \        }\n    };\n    struct speedup {\n        speedup() {\n            cin.tie(0);\n\
    \            ios::sync_with_stdio(0);\n        }\n    } speedup_instance;\n}\n\
    scanner::sca in() { return scanner::sca(); }\nscanner::vec in(int n) { return\
    \ scanner::vec(n); }\nscanner::mat in(int h, int w) { return scanner::mat(h, w);\
    \ }\n\nnamespace printer {\n    void precision(int d) {\n        cout << fixed\
    \ << setprecision(d);\n    }\n    void flush() {\n        cout.flush();\n    }\n\
    }\n\ntemplate < class T >\nostream& operator<<(ostream& os, const std::vector<\
    \ T > a) {\n    int n = a.size();\n    for(int i : rep(n)) { os << a[i]; if(i\
    \ != n - 1) os << ' '; }\n    return os;\n}\n\nint print() { cout << '\\n'; return\
    \ 0; }\ntemplate < class head, class... tail > int print(head&& h, tail&&... t)\
    \ {\n    cout << h; if(sizeof...(tail)) cout << ' ';\n    return print(forward<tail>(t)...);\n\
    }\n\ntemplate < class T > int print_n(const std::vector< T > a) {\n    int n =\
    \ a.size();\n    for(int i : rep(n)) cout << a[i] << \"\\n\";\n    return 0;\n\
    }\n#line 2 \"src/utility/key_val.hpp\"\ntemplate < class K, class V >\nstruct\
    \ key_val {\n    K key; V val;\n    key_val() {}\n    key_val(K key, V val) :\
    \ key(key), val(val) {}\n};\n#line 2 \"src/utility/vec_op.hpp\"\ntemplate < class\
    \ T >\nkey_val< int, T > max_of(const vector< T >& a) {\n    int i = max_element(a.begin(),\
    \ a.end()) - a.begin();\n    return {i, a[i]};\n}\n\ntemplate < class T >\nkey_val<\
    \ int, T > min_of(const vector< T >& a) {\n    int i = min_element(a.begin(),\
    \ a.end()) - a.begin();\n    return {i, a[i]};\n}\n\ntemplate < class T >\nT sum_of(const\
    \ vector< T >& a) {\n    T sum = 0;\n    for(const T x : a) sum += x;\n    return\
    \ sum;\n}\n\ntemplate < class T >\nvector<int> freq_of(const vector< T >& a, T\
    \ L, T R) {\n    vector<int> res(R - L);\n    for(const T x : a) res[x - L]++;\n\
    \    return res;\n}\n\ntemplate < class T >\nvector<int> freq_of(const vector<\
    \ T >& a, T R) {\n    return freq_of(a, T(0), R);\n}\n\ntemplate < class T >\n\
    struct prefix_sum {\n    vector< T > s;\n    prefix_sum(const vector< T >& a)\
    \ : s(a) {\n        s.insert(s.begin(), T(0));\n        for(int i : rep(a.size()))\
    \ s[i + 1] += s[i];\n    }\n    // [L, R)\n    T sum(int L, int R) {\n       \
    \ return s[R] - s[L];\n    }\n};\n#line 3 \"src/utility/heap.hpp\"\n\ntemplate\
    \ < class T > using heap_min = std::priority_queue< T, std::vector< T >, std::greater<\
    \ T > >;\ntemplate < class T > using heap_max = std::priority_queue< T, std::vector<\
    \ T >, std::less< T > >;\n\n#line 21 \"src/cp-template.hpp\"\n\n#line 1 \"src/algorithm/bin_search.hpp\"\
    \ntemplate < class T, class F >\nT bin_search(T ok, T ng, F& f) {\n    while(abs(ok\
    \ - ng) > 1) {\n        T mid = (ok + ng) / 2;\n        (f(mid) ? ok : ng) = mid;\n\
    \    }\n    return ok;\n}\n\ntemplate < class T, class F >\nT bin_search_real(T\
    \ ok, T ng, F& f, int step = 80) {\n    while(step--) {\n        T mid = (ok +\
    \ ng) / 2;\n        (f(mid) ? ok : ng) = mid;\n    }\n    return ok;\n}\n#line\
    \ 2 \"src/algorithm/argsort.hpp\"\n\ntemplate < class T > std::vector< int > argsort(const\
    \ std::vector< T > &a) {\n    std::vector< int > ids((int)a.size());\n    std::iota(ids.begin(),\
    \ ids.end(), 0);\n    std::sort(ids.begin(), ids.end(), [&](int i, int j) {\n\
    \        return a[i] < a[j] || (a[i] == a[j] && i < j);\n    });\n    return ids;\n\
    }\n#line 3 \"src/utility/random.hpp\"\n\nnamespace randnum {\n\nstatic uint seed;\n\
    static std::mt19937 mt;\nstruct gen_seed {\n    gen_seed() {\n        seed = std::random_device()();\n\
    \        mt = std::mt19937(seed);\n    }\n} gs;\n\n// [L, R)\ntemplate < class\
    \ T >\nT gen_int(T L, T R) {\n    return std::uniform_int_distribution< T >(L,\
    \ R - 1)(mt);\n}\n\ntemplate < class T >\nT get_real(T L, T R) {\n    return std::uniform_real_distribution<\
    \ T >(L, R)(mt);\n}\n\n}\n#line 3 \"src/utility/hash.hpp\"\n\ntemplate < int num_of_mod\
    \ = 2 >\nstruct hash_vector : public array<ll, num_of_mod> {\n    static constexpr\
    \ ll MODS[] = {999999937, 1000000007, 1000000009, 1000000021};\n    static_assert(1\
    \ <= num_of_mod and num_of_mod <= 4);\n    using array<ll, num_of_mod>::operator[];\n\
    \    using H = hash_vector;\n    static constexpr int n = num_of_mod;\n    hash_vector()\
    \ : array<ll,n>() {}\n    hash_vector(ll x) : H() { for(int i : rep(n)) (*this)[i]\
    \ = x % MODS[i]; }\n    H& operator+=(const H& rhs) { for(int i : rep(n)) if(((*this)[i]\
    \ += rhs[i]) >= MODS[i]) (*this)[i] -= MODS[i]; return *this; }\n    H& operator-=(const\
    \ H& rhs) { for(int i : rep(n)) if(((*this)[i] += MODS[i] - rhs[i]) >= MODS[i])\
    \ (*this)[i] -= MODS[i]; return *this; }\n    H& operator*=(const H& rhs) { for(int\
    \ i : rep(n)) (*this)[i] = (*this)[i] * rhs[i] % MODS[i]; return *this; }\n  \
    \  H& operator+=(const ll rhs) { for(int i : rep(n)) if(((*this)[i] += rhs % MODS[i])\
    \ >= MODS[i]) (*this)[i] -= MODS[i]; return *this; }\n    H& operator-=(const\
    \ ll rhs) { for(int i : rep(n)) if(((*this)[i] += MODS[i] - rhs % MODS[i]) >=\
    \ MODS[i]) (*this)[i] -= MODS[i]; return *this; }\n    H& operator*=(const ll\
    \ rhs) { for(int i : rep(n)) (*this)[i] = (*this)[i] * (rhs % MODS[i]) % MODS[i];\
    \ return *this; }\n    H operator+(const H& rhs) const { return H(*this) += rhs;\
    \ }\n    H operator-(const H& rhs) const { return H(*this) -= rhs; }\n    H operator*(const\
    \ H& rhs) const { return H(*this) *= rhs; }\n    H operator+(const ll rhs) const\
    \ { return H(*this) += rhs; }\n    H operator-(const ll rhs) const { return H(*this)\
    \ -= rhs; }\n    H operator*(const ll rhs) const { return H(*this) *= rhs; }\n\
    \    H operator-() const { return H().fill(0) - *this; }\n    friend H operator+(ll\
    \ x, const H& y) { return H(x) + y; }\n    friend H operator-(ll x, const H& y)\
    \ { return H(x) + y; }\n    friend H operator*(ll x, const H& y) { return H(x)\
    \ * y; }\n    bool operator==(const H& rhs) { for(int i : rep(n)) if((*this)[i]\
    \ != rhs[i]) return false; return true ; }\n    bool operator!=(const H& rhs)\
    \ { for(int i : rep(n)) if((*this)[i] != rhs[i]) return true ; return false; }\n\
    };\n#line 4 \"src/string/rolling_hash.hpp\"\n\ntemplate< int num_of_mod = 2 >\n\
    struct rolling_hash {\n\n    static const ll BASE;\n\n    vector< hash_vector<\
    \ num_of_mod > > pb, hs;\n    rolling_hash() {}\n    rolling_hash(const string&\
    \ s) {\n        int n = s.size();\n        hs.resize(n + 1); hs[0].fill(0);\n\
    \        pb.resize(n + 1); pb[0].fill(1);\n        for(int i : rep(n)) {\n   \
    \         hs[i + 1] = hs[i] * BASE + s[i];\n            pb[i + 1] = pb[i] * BASE;\n\
    \        }\n    }\n\n    // [l, r)\n    hash_vector< num_of_mod > get(int l, int\
    \ r) const {\n        return hs[r] - hs[l] * pb[r - l];\n    }\n\n    hash_vector<\
    \ num_of_mod > concat(hash_vector< num_of_mod > h1, hash_vector< num_of_mod >\
    \ h2, int h2_len) {\n        assert(0 <= h2_len and h2_len < int(pb.size()));\n\
    \        return h1 * pb[h2_len] + h2;\n    }\n\n    template < int n >\n    static\
    \ int lcp(const rolling_hash< n >& rh1, int l1, int r1, const rolling_hash< n\
    \ >& rh2, int l2, int r2) {\n        int lo = -1, hi = min(r1 - l1, r2 - l2) +\
    \ 1;\n        while(hi - lo > 1) {\n            int mid = (lo + hi) / 2;\n   \
    \         (rh1.get(l1, l1 + mid) == rh2.get(l2, l2 + mid) ? lo : hi) = mid;\n\
    \        }\n        return lo;\n    }\n\n    template < int n >\n    static int\
    \ cmp(const string& s1, const rolling_hash< n >& rh1, int l1, int r1,\n      \
    \             const string& s2, const rolling_hash< n >& rh2, int l2, int r2)\
    \ {\n        int len = lcp(rh1, l1, r1, rh2, l2, r2);\n        if(len == r1 -\
    \ l1 && len == r2 - l2) return 0;\n        if(len == r1 - l1) return -1;\n   \
    \     if(len == r2 - l2) return +1;\n        return (s1[l1 + len] < s2[l2 + len]\
    \ ? -1 : +1);\n    }\n};\n\ntemplate < int num_of_mod >\nconst ll rolling_hash<\
    \ num_of_mod >::BASE = randnum::gen_int<ll>(ll(0), hash_vector< num_of_mod >::MODS[0]);\n"
  code: "#include \"../../src/cp-template.hpp\"\n#include \"../../src/utility/random.hpp\"\
    \n#include \"../../src/utility/hash.hpp\"\n\ntemplate< int num_of_mod = 2 >\n\
    struct rolling_hash {\n\n    static const ll BASE;\n\n    vector< hash_vector<\
    \ num_of_mod > > pb, hs;\n    rolling_hash() {}\n    rolling_hash(const string&\
    \ s) {\n        int n = s.size();\n        hs.resize(n + 1); hs[0].fill(0);\n\
    \        pb.resize(n + 1); pb[0].fill(1);\n        for(int i : rep(n)) {\n   \
    \         hs[i + 1] = hs[i] * BASE + s[i];\n            pb[i + 1] = pb[i] * BASE;\n\
    \        }\n    }\n\n    // [l, r)\n    hash_vector< num_of_mod > get(int l, int\
    \ r) const {\n        return hs[r] - hs[l] * pb[r - l];\n    }\n\n    hash_vector<\
    \ num_of_mod > concat(hash_vector< num_of_mod > h1, hash_vector< num_of_mod >\
    \ h2, int h2_len) {\n        assert(0 <= h2_len and h2_len < int(pb.size()));\n\
    \        return h1 * pb[h2_len] + h2;\n    }\n\n    template < int n >\n    static\
    \ int lcp(const rolling_hash< n >& rh1, int l1, int r1, const rolling_hash< n\
    \ >& rh2, int l2, int r2) {\n        int lo = -1, hi = min(r1 - l1, r2 - l2) +\
    \ 1;\n        while(hi - lo > 1) {\n            int mid = (lo + hi) / 2;\n   \
    \         (rh1.get(l1, l1 + mid) == rh2.get(l2, l2 + mid) ? lo : hi) = mid;\n\
    \        }\n        return lo;\n    }\n\n    template < int n >\n    static int\
    \ cmp(const string& s1, const rolling_hash< n >& rh1, int l1, int r1,\n      \
    \             const string& s2, const rolling_hash< n >& rh2, int l2, int r2)\
    \ {\n        int len = lcp(rh1, l1, r1, rh2, l2, r2);\n        if(len == r1 -\
    \ l1 && len == r2 - l2) return 0;\n        if(len == r1 - l1) return -1;\n   \
    \     if(len == r2 - l2) return +1;\n        return (s1[l1 + len] < s2[l2 + len]\
    \ ? -1 : +1);\n    }\n};\n\ntemplate < int num_of_mod >\nconst ll rolling_hash<\
    \ num_of_mod >::BASE = randnum::gen_int<ll>(ll(0), hash_vector< num_of_mod >::MODS[0]);\n"
  dependsOn:
  - src/cp-template.hpp
  - src/utility/rep_itr.hpp
  - src/utility/io.hpp
  - src/utility/key_val.hpp
  - src/utility/vec_op.hpp
  - src/utility/heap.hpp
  - src/algorithm/bin_search.hpp
  - src/algorithm/argsort.hpp
  - src/utility/random.hpp
  - src/utility/hash.hpp
  isVerificationFile: false
  path: src/string/rolling_hash.hpp
  requiredBy: []
  timestamp: '2023-10-24 04:26:14+09:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/library_checker/string/rolling_hash.test.cpp
documentation_of: src/string/rolling_hash.hpp
layout: document
title: Rolling Hash
---

- $\mathrm{rolling}$ _ $\mathrm{hash}(s)$ : 文字列 $s$ のハッシュテーブルを構築する。 $O(\|S\|)$。
- $\mathrm{get}(l, r)$ : $s_{[l, r)}$ のハッシュ値を返す。 $O(1)$。  

長さが $O(N)$ の文字列 $s_1, s_2$ のハッシュテーブルをそれぞれ ${rh}_1, {rh}_2$ とする。
-  $\mathrm{lcp}({rh}_ 1, l_1, r_1, {rh}_ 2, l_2, r_2)$ : $s_{1[l_1, r_1)}$ と $s_{2[l_2, r_2)}$の最長共通接頭辞の長さを返す。 $O(\log N)$。
- $\mathrm{cmp}(s_1, {rh}_ 1, l_1, r_1, s_2, {rh}_ 2, l_2, r_2)$  :   文字列を辞書順比較し、次の値を返す。 $O(\log N)$。
	- $s_{1[l_1, r_1)} = s_{2[l_2, r_2)} \implies 0$
	- $s_{1[l_1, r_1)} < s_{2[l_2, r_2)} \implies +1$
	- $s_{1[l_1, r_1)} > s_{2[l_2, r_2)} \implies - 1$