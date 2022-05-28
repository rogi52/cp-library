template < class AbelGroup > class FenwickTree2D {
  public:
    using T = typename AbelGroup::set;

  private:
    int h,w;
    vector< vector< T > > data[2][2];

  public:
    FenwickTree2D(int h, int w) : h(h), w(w) {
        for(int i = 0; i < 2; i++)
            for(int j = 0; j < 2; j++)
                data[i][j].assign(h + 1, vector< T >(w + 1, AbelGroup::id));
    }

    // [0, x) * [0, y)
    void add(int x, int y, T v) {
        sub_add(0, 0, x, y, AbelGroup::pow(v, x * y));
        sub_add(1, 0, x, y, AbelGroup::pow(AbelGroup::inv(v), x));
        sub_add(0, 1, x, y, AbelGroup::pow(AbelGroup::inv(v), y));
        sub_add(1, 1, x, y, v);
    }

    // [x1, x2) * [y1, y2)
    void add(int x1, int x2, int y1, int y2, T v) {
        add(x1, y1, v);
        add(x1, y2, AbelGroup::inv(v));
        add(x2, y1, AbelGroup::inv(v));
        add(x2, y2, v);
    }

    // [0, x) * [0, y)
    T fold(int x, int y) {
        T s = AbelGroup::id;
        s = AbelGroup::op(sub_fold(0, 0, x, y), s);
        s = AbelGroup::op(AbelGroup::pow(sub_fold(1, 0, x, y), y), s);
        s = AbelGroup::op(AbelGroup::pow(sub_fold(0, 1, x, y), x), s);
        s = AbelGroup::op(AbelGroup::pow(sub_fold(1, 1, x, y), x * y), s);
        return s;
    }

    // [x1, x2) * [y1, y2)
    T fold(int x1, int x2, int y1, int y2) {
        T s = AbelGroup::id;
        s = AbelGroup::op(fold(x2, y2), s);
        s = AbelGroup::op(AbelGroup::inv(fold(x1, y2)), s);
        s = AbelGroup::op(AbelGroup::inv(fold(x2, y1)), s);
        s = AbelGroup::op(fold(x1, y1), s);
        return s;
    }

    T get(int x, int y) {
        return fold(x, x + 1, y, y + 1);
    }

  private:
    void sub_add(int f, int g, int x, int y, T v) {
        for(int i = x + 1; i <= h; i += i & -i)
            for(int j = y + 1; j <= w; j += j & -j)
                data[f][g][i - 1][j - 1] = AbelGroup::op(data[f][g][i - 1][j - 1], v);
    }

    T sub_fold(int f, int g, int x, int y) {
        T s = AbelGroup::id;
        for(int i = x; i > 0; i -= i & -i)
            for(int j = y; j > 0; j -= j & -j)
                s = AbelGroup::op(data[f][g][i - 1][j - 1], s);
        return s;
    }
};
