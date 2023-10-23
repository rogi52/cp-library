#pragma once
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using ld = long double;
using uint = unsigned int;
using ull  = unsigned long long;
using i32 = int;
using u32 = unsigned int;
using i64 = long long;
using u64 = unsigned long long;
using i128 = __int128_t;
template < class T > bool chmin(T& a, T b) { if(a > b) { a = b; return true; } return false; }
template < class T > bool chmax(T& a, T b) { if(a < b) { a = b; return true; } return false; }

#include "./utility/rep_itr.hpp"
#include "./utility/io.hpp"
#include "./utility/key_val.hpp"
#include "./utility/vec_op.hpp"
#include "./utility/heap.hpp"

#include "./algorithm/bin_search.hpp"
#include "./algorithm/argsort.hpp"
