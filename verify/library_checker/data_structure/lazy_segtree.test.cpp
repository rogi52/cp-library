#define PROBLEM "https://judge.yosupo.jp/problem/range_affine_range_sum"

#include "src/cp-template.hpp"
#include "src/number/modint.hpp"
#include "src/data_structure/lazy_segtree.hpp"
#include "src/algebra/range_affine_range_sum.hpp"

int main(){
    cin.tie(0);
    ios::sync_with_stdio(0);
    
    int N,Q; cin >> N >> Q;
    using S = range_affine_range_sum<mint>::value_structure::set;
    vector< S > a(N);
    rep(i,N) {
        mint x; cin >> x;
        a[i] = S{x, 1};
    }

    lazy_segtree< range_affine_range_sum<mint> > lzst(a);
    rep(_,Q) {
        int t; cin >> t;
        if(t == 0) {
            int l,r; mint b,c; cin >> l >> r >> b >> c;
            lzst.apply(l, r, {b, c});
        }

        if(t == 1) {
            int l,r; cin >> l >> r;
            cout << lzst.prod(l, r).first << "\n";
        }
    }
}