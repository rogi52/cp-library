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
  bundledCode: "Traceback (most recent call last):\n  File \"/opt/hostedtoolcache/Python/3.10.6/x64/lib/python3.10/site-packages/onlinejudge_verify/documentation/build.py\"\
    , line 71, in _render_source_code_stat\n    bundled_code = language.bundle(stat.path,\
    \ basedir=basedir, options={'include_paths': [basedir]}).decode()\n  File \"/opt/hostedtoolcache/Python/3.10.6/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus.py\"\
    , line 187, in bundle\n    bundler.update(path)\n  File \"/opt/hostedtoolcache/Python/3.10.6/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 401, in update\n    self.update(self._resolve(pathlib.Path(included), included_from=path))\n\
    \  File \"/opt/hostedtoolcache/Python/3.10.6/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 260, in _resolve\n    raise BundleErrorAt(path, -1, \"no such header\"\
    )\nonlinejudge_verify.languages.cplusplus_bundle.BundleErrorAt: src/union_find.hpp:\
    \ line -1: no such header\n"
  code: "#include \"src/union_find.hpp\"\n\ntemplate < class abel_monoid > class union_fold\
    \ : union_find {\n  public:\n    union_fold(vector< T > &a) : union_find(a.size()),\
    \ value(a) {}\n    int unite(int x, int y) {\n        if((x = union_find::unite(x,\
    \ y)) != -1)\n            value[x] = abel_monoid::op(value[x], value[y]);\n  \
    \      return x;\n    }\n    T fold(int x) { return value[root(x)]; }\n    void\
    \ add(int x, T v) {\n        x = root(x);\n        value[x] = abel_monoid::op(value[x],\
    \ v);\n    }\n\n  private:\n    using T = typename abel_monoid::set;\n    vector<\
    \ T > value;\n};\n"
  dependsOn: []
  isVerificationFile: false
  path: src/data_structure/union_fold.hpp
  requiredBy: []
  timestamp: '1970-01-01 00:00:00+00:00'
  verificationStatus: LIBRARY_NO_TESTS
  verifiedWith: []
documentation_of: src/data_structure/union_fold.hpp
layout: document
redirect_from:
- /library/src/data_structure/union_fold.hpp
- /library/src/data_structure/union_fold.hpp.html
title: src/data_structure/union_fold.hpp
---