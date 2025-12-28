#include <bits/stdc++.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
using namespace std;

// ============================================================
// Wheel105 compressed sieve for pi(N), multi-thread by segments
// - odd index z = (n-1)/2
// - keep only z with gcd(2z+1,105)=1 => 48 residues per 105
// - mark composites using next0/next1 jump tables on residues
// - count primes by popcount(~bitset)
// - OpenMP parallelization: each segment independent, thread-local bitset
// ============================================================

static inline int64_t isqrt_i64(int64_t x) {
    int64_t r = (int64_t)std::sqrt((long double)x);
    while ((r + 1) > 0 && (r + 1) * (r + 1) <= x) r++;
    while (r * r > x) r--;
    return r;
}

struct Wheel105 {
    static constexpr int M = 105;
    static constexpr int K = 48;

    array<uint8_t, M> valid{};
    array<uint16_t, M + 1> prefix_lt{};
    array<array<uint8_t, M>, M> next0{};
    array<array<uint8_t, M>, M> next1{};

    Wheel105() {
        for (int r = 0; r < M; r++) {
            int n = 2 * r + 1;
            if (std::gcd(n, M) == 1) valid[r] = 1;
        }
        int c = 0;
        for (int r = 0; r < M; r++) {
            prefix_lt[r] = (uint16_t)c;
            if (valid[r]) c++;
        }
        prefix_lt[M] = (uint16_t)c;
        if (c != K) throw runtime_error("Wheel105: expected 48 valid residues.");

        for (int a = 0; a < M; a++) {
            for (int r = 0; r < M; r++) {
                int t = 0, rr = r;
                while (!valid[rr] && t < M) { t++; rr += a; rr %= M; }
                next0[a][r] = (uint8_t)((t < M) ? t : 0);

                t = 1; rr = (r + a) % M;
                while (!valid[rr] && t < M) { t++; rr += a; rr %= M; }
                next1[a][r] = (uint8_t)((t <= M) ? t : 1);
            }
        }
    }

    inline int64_t count_valid_before_z(int64_t z) const {
        int64_t block = z / M;
        int rem = (int)(z - block * M);
        return block * K + prefix_lt[rem];
    }
    inline int64_t global_bit_index(int64_t z) const {
        int64_t block = z / M;
        int rem = (int)(z - block * M);
        return block * K + prefix_lt[rem];
    }
};

static inline int popcount_u64(uint64_t x) {
    return __builtin_popcountll(x);
}

static vector<int> base_primes_upto(int limit) {
    vector<int> primes;
    if (limit < 2) return primes;
    primes.push_back(2);
    if (limit < 3) return primes;

    int n_odds = (limit - 1) / 2;
    vector<uint8_t> comp(n_odds, 0);

    int r = (int)std::sqrt((long double)limit);
    int zmax = (r - 1) / 2;

    for (int z = 1; z <= zmax; z++) {
        if (!comp[z]) {
            int p = 2 * z + 1;
            int start = (p * p - 1) / 2;
            for (int i = start; i < n_odds; i += p) comp[i] = 1;
        }
    }
    for (int z = 1; z < n_odds; z++) if (!comp[z]) primes.push_back(2 * z + 1);
    return primes;
}

// Count primes (in compressed domain) in a single segment [zL..zR] (z-space)
// Returns count of primes in this segment (>= 11), i.e. excludes 2,3,5,7.
static inline int64_t count_segment(
    int64_t zL, int64_t zR,
    int64_t low_n, int64_t high_n,
    const Wheel105& W,
    const vector<int>& primes,
    vector<uint64_t>& compw // thread-local reusable buffer
) {
    int64_t base_rank = W.count_valid_before_z(zL);
    int64_t end_rank  = W.count_valid_before_z(zR + 1);
    int64_t seg_bits  = end_rank - base_rank;
    if (seg_bits <= 0) return 0;

    int64_t nwords = (seg_bits + 63) >> 6;
    if ((int64_t)compw.size() < nwords) compw.resize((size_t)nwords);
    std::fill(compw.begin(), compw.begin() + (size_t)nwords, 0ULL);

    // mark composites for p>7
    for (size_t k = 0; k < primes.size(); k++) {
        int p = primes[k];
        if (p <= 7) continue;

        int64_t pp = 1LL * p * p;
        if (pp > high_n) break;

        int64_t start = pp;
        if (start < low_n) start = ((low_n + p - 1) / p) * 1LL * p;
        if ((start & 1LL) == 0) start += p;
        if (start > high_n) continue;

        int64_t z = (start - 1) / 2;
        if (z < zL) z = zL;

        int a = p % Wheel105::M;
        int res = (int)(z % Wheel105::M);

        if (!W.valid[res]) {
            int t0 = W.next0[a][res];
            z += 1LL * p * t0;
            res = (res + a * t0) % Wheel105::M;
        }

        while (z <= zR) {
            int64_t gb  = W.global_bit_index(z);
            int64_t loc = gb - base_rank;
            compw[(size_t)(loc >> 6)] |= (1ULL << (loc & 63));

            int t = W.next1[a][res];
            z += 1LL * p * t;
            res = (res + a * t) % Wheel105::M;
        }
    }

    // count primes = number of 0 bits among seg_bits
    int64_t cnt = 0;
    int64_t full_words = seg_bits >> 6;
    for (int64_t w = 0; w < full_words; w++) cnt += popcount_u64(~compw[(size_t)w]);

    int rem = (int)(seg_bits & 63);
    if (rem) {
        uint64_t last = ~compw[(size_t)full_words];
        last &= ((rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL));
        cnt += popcount_u64(last);
    }
    return cnt;
}

// Full pi(N), wheel105 compressed, OpenMP by segments
int64_t pi_wheel105_omp(int64_t N, int64_t seg_mb, int threads) {
    if (N < 2) return 0;

    Wheel105 W;

    int64_t lim = isqrt_i64(N);
    auto primes = base_primes_upto((int)lim);

    // seg_bits from MB
    int64_t seg_bits_cap = seg_mb * 1024LL * 1024LL * 8LL;
    if (seg_bits_cap < (1 << 20)) seg_bits_cap = (1 << 20);

    // density 48/105 => seg_z_span ≈ seg_bits * (105/48)
    int64_t seg_z_span = (seg_bits_cap * Wheel105::M) / Wheel105::K;
    if (seg_z_span < (1 << 20)) seg_z_span = (1 << 20);

    // count small primes explicitly
    int64_t total = 0;
    if (N >= 2) total++;
    if (N >= 3) total++;
    if (N >= 5) total++;
    if (N >= 7) total++;
    if (N < 11) return total;

    int64_t z_start = (11 - 1) / 2;
    int64_t z_end   = (N  - 1) / 2;

    // build segment boundaries
    vector<int64_t> segL;
    for (int64_t zL = z_start; zL <= z_end; ) {
        segL.push_back(zL);
        int64_t zR = zL + seg_z_span - 1;
        if (zR > z_end) zR = z_end;
        zL = zR + 1;
    }
    const int64_t nseg = (int64_t)segL.size();

    // thread-local buffers
#ifdef _OPENMP
    if (threads <= 0) threads = omp_get_max_threads();
    threads = max(1, threads);
    omp_set_num_threads(threads);
#else
    threads = 1;
#endif
    vector<vector<uint64_t>> buffers((size_t)threads);

    int64_t sum_segments = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(+:sum_segments)
#endif
    for (int64_t i = 0; i < nseg; i++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        int64_t zL = segL[(size_t)i];
        int64_t zR = zL + seg_z_span - 1;
        if (zR > z_end) zR = z_end;

        int64_t low_n  = 2 * zL + 1;
        int64_t high_n = 2 * zR + 1;

        sum_segments += count_segment(zL, zR, low_n, high_n, W, primes, buffers[(size_t)tid]);
    }

    total += sum_segments;
    return total;
}

// simple bench helper (median of runs)
static double median_time_s(function<void(void)> fn, int runs) {
    vector<double> t; t.reserve(runs);
    for (int i=0;i<runs;i++){
        auto a = chrono::high_resolution_clock::now();
        fn();
        auto b = chrono::high_resolution_clock::now();
        t.push_back(chrono::duration<double>(b-a).count());
    }
    nth_element(t.begin(), t.begin()+runs/2, t.end());
    return t[runs/2];
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int64_t N = 100000000;
    int64_t seg_mb = 16;
    int threads = 4;
    int runs = 3;

    if (argc >= 2) N = atoll(argv[1]);
    if (argc >= 3) seg_mb = atoll(argv[2]);
    if (argc >= 4) threads = atoi(argv[3]);
    if (argc >= 5) runs = atoi(argv[4]);

    int64_t ans = 0;
    auto fn = [&](){ ans = pi_wheel105_omp(N, seg_mb, threads); };

    double med = median_time_s(fn, runs);

    cout << "N=" << N << " seg_mb=" << seg_mb << " threads=" << threads << " runs=" << runs << "\n";
    cout << "pi(N)=" << ans << "\n";
    cout << "median_time_s=" << fixed << setprecision(6) << med << "\n";
    return 0;
}
