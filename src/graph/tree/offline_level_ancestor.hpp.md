---
data:
  _extendedDependsOn: []
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':warning:'
  attributes:
    links: []
  bundledCode: "#line 1 \"src/graph/tree/offline_level_ancestor.hpp\"\nnamespace off_la\
    \ {\n\nusing func = function<void(int, vector<int>&, vector<int>&)>;\nvoid sol(const\
    \ vector<vector<int>> &tree, int root, const vector<int> &query, func &ans) {\n\
    \    int N = tree.size();\n    vector<int> path, pos(N, -1);\n    vector<vector<int>>\
    \ q(N);\n    rep(i,int(query.size())) q[query[i]].push_back(i);\n\n    function<void(int,int)>\
    \ dfs = [&](int v, int p) -> void {\n        pos[v] = path.size();\n        path.push_back(v);\n\
    \        for(int i : q[v]) ans(i, path, pos);\n        for(int to : tree[v]) if(to\
    \ != p) dfs(to, v);\n        pos[v] = -1;\n        path.pop_back();\n    }; dfs(root,\
    \ -1);\n}\n\n}\n"
  code: "namespace off_la {\n\nusing func = function<void(int, vector<int>&, vector<int>&)>;\n\
    void sol(const vector<vector<int>> &tree, int root, const vector<int> &query,\
    \ func &ans) {\n    int N = tree.size();\n    vector<int> path, pos(N, -1);\n\
    \    vector<vector<int>> q(N);\n    rep(i,int(query.size())) q[query[i]].push_back(i);\n\
    \n    function<void(int,int)> dfs = [&](int v, int p) -> void {\n        pos[v]\
    \ = path.size();\n        path.push_back(v);\n        for(int i : q[v]) ans(i,\
    \ path, pos);\n        for(int to : tree[v]) if(to != p) dfs(to, v);\n       \
    \ pos[v] = -1;\n        path.pop_back();\n    }; dfs(root, -1);\n}\n\n}\n"
  dependsOn: []
  isVerificationFile: false
  path: src/graph/tree/offline_level_ancestor.hpp
  requiredBy: []
  timestamp: '2023-03-28 01:17:00+09:00'
  verificationStatus: LIBRARY_NO_TESTS
  verifiedWith: []
documentation_of: src/graph/tree/offline_level_ancestor.hpp
layout: document
title: Offline Level Ancestor
---

参考: https://noshi91.hatenablog.com/entry/2019/09/22/114149

## クエリ例
$\mathrm{LA}(v, d) \coloneqq v$ の祖先であって、深さが $d$ である頂点
```cpp
off_la::func la = [&](int i, vector<int> &path, vector<int> &pos) {
    ans[i] = path[d[i]];
};
```

$v$ から根の方向に $k$ 回親を辿った頂点
```cpp
off_la::func up = [&](int i, vector<int> &path, vector<int> &pos) {
    ans[i] = path[int(path.size()) - 1 - k[i]];
};
```

$x$ は $y$ の祖先とし、距離を $d$ とする。 $x$ と $y$ を結ぶパス上の頂点を順に $x = v_0, v_1, \dots, v_d = y$ とするときの $v_1$

https://atcoder.jp/contests/pakencamp-2022-day1/tasks/pakencamp_2022_day1_g
```cpp
off_la::func get_ans = [&](int i, vector<int> &path, vector<int> &pos) {
    ans[i] = path[pos[x[i]] + 1];
};
```