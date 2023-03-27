---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: src/cp-template.hpp
    title: src/cp-template.hpp
  - icon: ':heavy_check_mark:'
    path: src/data_structure/cht_offline_get_min.hpp
    title: src/data_structure/cht_offline_get_min.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/segment_add_get_min
    links:
    - https://judge.yosupo.jp/problem/segment_add_get_min
  bundledCode: "#line 1 \"verify/library_checker/data_structure/segment_add_get_min.test.cpp\"\
    \n#define PROBLEM \"https://judge.yosupo.jp/problem/segment_add_get_min\"\n\n\
    #line 1 \"src/cp-template.hpp\"\n#include <bits/stdc++.h>\n#define rep(i,n) for(int\
    \ i = 0; i < (n); i++)\nusing namespace std;\nusing ll = long long;\nusing ld\
    \ = long double;\nusing uint = unsigned int;\nusing ull  = unsigned long long;\n\
    #line 1 \"src/data_structure/cht_offline_get_min.hpp\"\ntemplate < class T > class\
    \ CHT_offline_get_min {\n  private:\n    struct Line {\n        T a, b;\n    \
    \    Line(T a, T b) : a(a), b(b) {}\n        T eval(T x) { return a * x + b; }\n\
    \    };\n\n    T sgn;\n    int n;\n    vector< Line > ls;\n    vector< T > xs;\n\
    \n  public:\n    T inf = numeric_limits< T >::max();\n\n    CHT_offline_get_min(vector<\
    \ T > &x, bool is_min = true) : xs(x) {\n        sort(xs.begin(), xs.end());\n\
    \        xs.erase(unique(xs.begin(), xs.end()), xs.end());\n        n = xs.size();\n\
    \        ls.assign(n << 1, Line(0, inf));\n        sgn = is_min ? +1 : -1;\n \
    \   }\n\n    void add_line(T a, T b) { update(a, b, 0, n); }\n\n    void add_segment(T\
    \ a, T b, T l, T r) {\n        int xl = distance(xs.begin(), lower_bound(xs.begin(),\
    \ xs.end(), l));\n        int xr = distance(xs.begin(), lower_bound(xs.begin(),\
    \ xs.end(), r));\n        update(a, b, xl, xr);\n    }\n\n    T query(T x) {\n\
    \        int i = distance(xs.begin(), lower_bound(xs.begin(), xs.end(), x));\n\
    \        assert(i != n && x == xs[i]);\n        T v = inf;\n        for(i += n;\
    \ i > 0; i >>= 1) v = min(v, ls[i].eval(x));\n        return sgn * v;\n    }\n\
    \n  private:\n    void update(T a, T b, int l, int r) {\n        a *= sgn, b *=\
    \ sgn;\n        Line f(a, b);\n        for(l += n, r += n; l < r; l >>= 1, r >>=\
    \ 1) {\n            if(l & 1) descend(f, l++);\n            if(r & 1) descend(f,\
    \ --r);\n        }\n    }\n\n    void descend(Line g, int i) {\n        int l\
    \ = i, r = i + 1;\n        while(l < n) l <<= 1, r <<= 1;\n        while(l < r)\
    \ {\n            int m = (l + r) >> 1;\n            T xl = xs[l - n], xm = xs[m\
    \ - n], xr = xs[r - 1 - n];\n            Line &f = ls[i];\n            if(f.eval(xl)\
    \ <= g.eval(xl) && f.eval(xr) <= g.eval(xr)) return;\n            if(f.eval(xl)\
    \ >= g.eval(xl) && f.eval(xr) >= g.eval(xr)) { f = g; return; }\n            if(f.eval(xm)\
    \ >  g.eval(xm)) swap(f, g);\n            if(f.eval(xl) >  g.eval(xl)) i = i <<\
    \ 1 | 0, r = m; else i = i << 1 | 1, l = m;\n        }\n    }\n};\n#line 5 \"\
    verify/library_checker/data_structure/segment_add_get_min.test.cpp\"\n\nint main(){\n\
    \    cin.tie(0);\n    ios::sync_with_stdio(0);\n\n    int N,Q; cin >> N >> Q;\n\
    \    using Line = tuple<ll,ll,ll,ll>;\n    vector< Line > lines(N);\n    for(auto\
    \ &[l, r, a, b] : lines) cin >> l >> r >> a >> b;\n\n    vector< pair< int, Line\
    \ > > query(Q);\n    vector< ll > xs;\n    rep(i,Q) {\n        int t; cin >> t;\n\
    \        if(t == 0) {\n            ll l, r, a, b; cin >> l >> r >> a >> b;\n \
    \           query[i] = {t, {l, r, a, b}};\n        } else {\n            ll x;\
    \ cin >> x;\n            query[i] = {t, {x, 0, 0, 0}};\n            xs.push_back(x);\n\
    \        }\n    }\n\n    CHT_offline_get_min<ll> cht(xs);\n    for(auto [l, r,\
    \ a, b] : lines) cht.add_segment(a, b, l, r);\n    rep(i,Q) {\n        int t =\
    \ query[i].first;\n        if(t == 0) {\n            auto [l, r, a, b] = query[i].second;\n\
    \            cht.add_segment(a, b, l, r);\n        } else {\n            auto\
    \ [x, _, __, ___] = query[i].second;\n            ll ans = cht.query(x);\n   \
    \         if(ans == cht.inf) {\n                cout << \"INFINITY\" << '\\n';\n\
    \            } else {\n                cout << ans << '\\n';\n            }\n\
    \        }\n    }\n}\n"
  code: "#define PROBLEM \"https://judge.yosupo.jp/problem/segment_add_get_min\"\n\
    \n#include \"src/cp-template.hpp\"\n#include \"src/data_structure/cht_offline_get_min.hpp\"\
    \n\nint main(){\n    cin.tie(0);\n    ios::sync_with_stdio(0);\n\n    int N,Q;\
    \ cin >> N >> Q;\n    using Line = tuple<ll,ll,ll,ll>;\n    vector< Line > lines(N);\n\
    \    for(auto &[l, r, a, b] : lines) cin >> l >> r >> a >> b;\n\n    vector< pair<\
    \ int, Line > > query(Q);\n    vector< ll > xs;\n    rep(i,Q) {\n        int t;\
    \ cin >> t;\n        if(t == 0) {\n            ll l, r, a, b; cin >> l >> r >>\
    \ a >> b;\n            query[i] = {t, {l, r, a, b}};\n        } else {\n     \
    \       ll x; cin >> x;\n            query[i] = {t, {x, 0, 0, 0}};\n         \
    \   xs.push_back(x);\n        }\n    }\n\n    CHT_offline_get_min<ll> cht(xs);\n\
    \    for(auto [l, r, a, b] : lines) cht.add_segment(a, b, l, r);\n    rep(i,Q)\
    \ {\n        int t = query[i].first;\n        if(t == 0) {\n            auto [l,\
    \ r, a, b] = query[i].second;\n            cht.add_segment(a, b, l, r);\n    \
    \    } else {\n            auto [x, _, __, ___] = query[i].second;\n         \
    \   ll ans = cht.query(x);\n            if(ans == cht.inf) {\n               \
    \ cout << \"INFINITY\" << '\\n';\n            } else {\n                cout <<\
    \ ans << '\\n';\n            }\n        }\n    }\n}\n"
  dependsOn:
  - src/cp-template.hpp
  - src/data_structure/cht_offline_get_min.hpp
  isVerificationFile: true
  path: verify/library_checker/data_structure/segment_add_get_min.test.cpp
  requiredBy: []
  timestamp: '2023-03-26 03:29:33+09:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/library_checker/data_structure/segment_add_get_min.test.cpp
layout: document
redirect_from:
- /verify/verify/library_checker/data_structure/segment_add_get_min.test.cpp
- /verify/verify/library_checker/data_structure/segment_add_get_min.test.cpp.html
title: verify/library_checker/data_structure/segment_add_get_min.test.cpp
---