template < class T > class sum_monoid {
  public:
    using set = T;
    static constexpr T op(const T &l, const T &r) { return l + r; }
    static constexpr T id() { return T(0); }
    static constexpr T inv(const T &x) { return -x; }
    static constexpr T pow(const T &x, const ll n) { return x * n; }
    static constexpr bool comm = true;
};
