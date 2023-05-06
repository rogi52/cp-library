---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: src/cp-template.hpp
    title: src/cp-template.hpp
  - icon: ':heavy_check_mark:'
    path: src/graph/tree/cartesian_tree.hpp
    title: src/graph/tree/cartesian_tree.hpp
  - icon: ':heavy_check_mark:'
    path: src/utility/io.hpp
    title: src/utility/io.hpp
  - icon: ':heavy_check_mark:'
    path: src/utility/rep_itr.hpp
    title: src/utility/rep_itr.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/cartesian_tree
    links:
    - https://judge.yosupo.jp/problem/cartesian_tree
  bundledCode: "#line 1 \"verify/library_checker/graph/tree/cartesian_tree.test.cpp\"\
    \n#define PROBLEM \"https://judge.yosupo.jp/problem/cartesian_tree\"\n\n#line\
    \ 1 \"src/cp-template.hpp\"\n#include <bits/stdc++.h>\nusing namespace std;\n\
    using ll = long long;\nusing ld = long double;\nusing uint = unsigned int;\nusing\
    \ ull  = unsigned long long;\nusing i128 = __int128_t;\ntemplate < class T > bool\
    \ chmin(T& a, T b) { if(a > b) { a = b; return true; } return false; }\ntemplate\
    \ < class T > bool chmax(T& a, T b) { if(a < b) { a = b; return true; } return\
    \ false; }\n\n#line 1 \"src/utility/rep_itr.hpp\"\ntemplate < class T > struct\
    \ itr {\n    T i, d;\n    constexpr itr(const T i) noexcept : i(i), d(1) {}\n\
    \    constexpr itr(const T i, const T d) noexcept : i(i), d(d) {}\n    void operator++()\
    \ noexcept { i += d; }\n    constexpr int operator*() const noexcept { return\
    \ i; }\n    constexpr bool operator!=(const itr x) const noexcept {\n        return\
    \ d > 0 ? i < x.i : i > x.i;\n    }\n};\n\ntemplate < class T > struct rep {\n\
    \    const itr< T > s, t;\n    constexpr rep(const T t) noexcept : s(0), t(t)\
    \ {}\n    constexpr rep(const T s, const T t) noexcept : s(s), t(t) {}\n    constexpr\
    \ rep(const T s, const T t, const T d) noexcept : s(s, d), t(t, d) {}\n    constexpr\
    \ auto begin() const noexcept { return s; }\n    constexpr auto end() const noexcept\
    \ { return t; }\n};\n\ntemplate < class T > struct revrep {\n    const itr < T\
    \ > s, t;\n    constexpr revrep(const T t) noexcept : s(t - 1, -1), t(-1, -1)\
    \ {}\n    constexpr revrep(const T s, const T t) noexcept : s(t - 1, -1), t(s\
    \ - 1, -1) {}\n    constexpr revrep(const T s, const T t, const T d) noexcept\
    \ : s(t - 1, -d), t(s - 1, -d) {}\n    constexpr auto begin() const noexcept {\
    \ return s; }\n    constexpr auto end() const noexcept { return t; }\n};\n#line\
    \ 1 \"src/utility/io.hpp\"\nnamespace scanner {\n    struct sca {\n        template\
    \ < class T > operator T() {\n            T s; cin >> s; return s;\n        }\n\
    \    };\n    struct vec {\n        int n;\n        vec(int n) : n(n) {}\n    \
    \    template < class T > operator vector< T >() {\n            vector< T > v(n);\n\
    \            for(T& x : v) cin >> x;\n            return v;\n        }\n    };\n\
    \    struct mat {\n        int h,w;\n        mat(int h, int w) : h(h), w(w) {}\n\
    \        template < class T > operator vector< vector< T > >() {\n           \
    \ vector m(h, vector< T >(w));\n            for(vector< T >& v : m) for(T& x :\
    \ v) cin >> x;\n            return m;\n        }\n    };\n    struct speedup {\n\
    \        speedup() {\n            cin.tie(0);\n            ios::sync_with_stdio(0);\n\
    \        }\n    } su;\n}\nscanner::sca in() { return scanner::sca(); }\nscanner::vec\
    \ in(int n) { return scanner::vec(n); }\nscanner::mat in(int h, int w) { return\
    \ scanner::mat(h, w); }\n\nnamespace printer {\n    void precision(int d) {\n\
    \        cout << fixed << setprecision(d);\n    }\n    void flush() {\n      \
    \  cout.flush();\n    }\n}\nint print() { cout << '\\n'; return 0; }\ntemplate\
    \ < class head, class... tail > int print(head&& h, tail&&... t) {\n    cout <<\
    \ h; if(sizeof...(tail)) cout << ' ';\n    return print(forward<tail>(t)...);\n\
    }\ntemplate < class T > int print(vector< T >& a, char sep = ' ') {\n    int n\
    \ = a.size();\n    for(int i : rep(n)) cout << a[i] << (i != n - 1 ? sep : '\\\
    n');\n    return 0;\n}\ntemplate < class T > int print(vector< vector< T > >&\
    \ a) {\n    if(a.empty()) return 0;\n    int h = a.size(), w = a[0].size();\n\
    \    for(int i : rep(h)) for(int j : rep(w)) cout << a[i][j] << (j != w - 1 ?\
    \ ' ' : '\\n');\n    return 0;\n}\n#line 2 \"src/graph/tree/cartesian_tree.hpp\"\
    \n\ntemplate < class T >\npair< vector<vector<int>>, int > cartesian_tree(vector<\
    \ T > &a) {\n    int n = a.size();\n    vector<vector<int>> tree(n);\n    vector<int>\
    \ parent(n, -1);\n    stack<int> st;\n    for(int i = 0; i < n; i++) {\n     \
    \   int prev = -1;\n        while(!st.empty() && a[i] < a[st.top()]) { prev =\
    \ st.top(); st.pop(); }\n        if(prev != -1) parent[prev] = i;\n        if(!st.empty())\
    \ parent[i] = st.top();\n        st.push(i);\n    }\n\n    int root = -1;\n  \
    \  for(int i = 0; i < n; i++) {\n        if(parent[i] == -1) root = i;\n     \
    \   else tree[parent[i]].push_back(i);\n    }\n    return {tree, root};\n}\n#line\
    \ 5 \"verify/library_checker/graph/tree/cartesian_tree.test.cpp\"\n\nint main(){\n\
    \    cin.tie(0);\n    ios::sync_with_stdio(0);\n    \n    int n; cin >> n;\n \
    \   vector<int> a = in(n);\n    auto [G, r] = cartesian_tree<int>(a);\n    vector<int>\
    \ ans(n);\n    ans[r] = r;\n    for(int i : rep(n)) for(int v : G[i]) ans[v] =\
    \ i;\n    for(int p : ans) cout << p << \" \"; cout << endl;\n};\n"
  code: "#define PROBLEM \"https://judge.yosupo.jp/problem/cartesian_tree\"\n\n#include\
    \ \"src/cp-template.hpp\"\n#include \"src/graph/tree/cartesian_tree.hpp\"\n\n\
    int main(){\n    cin.tie(0);\n    ios::sync_with_stdio(0);\n    \n    int n; cin\
    \ >> n;\n    vector<int> a = in(n);\n    auto [G, r] = cartesian_tree<int>(a);\n\
    \    vector<int> ans(n);\n    ans[r] = r;\n    for(int i : rep(n)) for(int v :\
    \ G[i]) ans[v] = i;\n    for(int p : ans) cout << p << \" \"; cout << endl;\n\
    };\n"
  dependsOn:
  - src/cp-template.hpp
  - src/utility/rep_itr.hpp
  - src/utility/io.hpp
  - src/graph/tree/cartesian_tree.hpp
  isVerificationFile: true
  path: verify/library_checker/graph/tree/cartesian_tree.test.cpp
  requiredBy: []
  timestamp: '2023-05-06 10:51:50+09:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/library_checker/graph/tree/cartesian_tree.test.cpp
layout: document
redirect_from:
- /verify/verify/library_checker/graph/tree/cartesian_tree.test.cpp
- /verify/verify/library_checker/graph/tree/cartesian_tree.test.cpp.html
title: verify/library_checker/graph/tree/cartesian_tree.test.cpp
---