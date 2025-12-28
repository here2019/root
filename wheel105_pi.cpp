#include <bits/stdc++.h>
using namespace std;

// ============================================================
// Wheel105 compressed sieve for pi(N)
// - represent odd numbers via z = (n-1)/2
// - compress z by keeping only residues r mod 105 such that gcd(2r+1,105)=1
// - mark composites using jump tables next0/next1 on residues
// - count primes with popcount on ~bitset
//
// Counts primes <= N (pi(N)).
// Handles primes 2,3,5,7 explicitly.
// ============================================================

static inline int64_t isqrt_i64(int64_t x) {
    return (int64_t)std::sqrt((long double)x);
}

struct Wheel105 {
    static constexpr int M = 105;
    static constexpr int K = 48;

    array<uint8_t, M> valid{};
    array<uint16_t, M+1> prefix_lt{};
    // next0[a][r] minimal t>=0 to reach valid
    // next1[a][r] minimal t>=1 to reach valid
    array<array<uint8_t, M>, M> next0{};
    array<array<uint8_t, M>, M> next1{};

    Wheel105() {
        // build valid residues
        for (int r=0; r<M; r++){
            int n = 2*r + 1;
            if (std::gcd(n, M) == 1) valid[r] = 1;
        }
        // prefix_lt
        int c=0;
        for(int r=0;r<M;r++){
            prefix_lt[r] = (uint16_t)c;
            if(valid[r]) c++;
        }
        prefix_lt[M] = (uint16_t)c;
        if(c != K) {
            throw runtime_error("Wheel105: expected 48 valid residues.");
        }
        // next tables
        for(int a=0;a<M;a++){
            for(int r=0;r<M;r++){
                // next0
                int t=0;
                int rr=r;
                while(!valid[rr] && t < M){
                    t++;
                    rr += a;
                    rr %= M;
                }
                next0[a][r] = (uint8_t)((t < M) ? t : 0);

                // next1
                t=1;
                rr = (r + a) % M;
                while(!valid[rr] && t < M){
                    t++;
                    rr += a;
                    rr %= M;
                }
                next1[a][r] = (uint8_t)((t <= M) ? t : 1);
            }
        }
    }

    // count valid indices in [0..z-1]
    inline int64_t count_valid_before_z(int64_t z) const {
        int64_t block = z / M;
        int rem = (int)(z - block * M);
        return block * K + prefix_lt[rem];
    }

    // global compressed bit index for valid z
    inline int64_t global_bit_index(int64_t z) const {
        int64_t block = z / M;
        int rem = (int)(z - block * M);
        return block * K + prefix_lt[rem];
    }
};

static vector<int> base_primes_upto(int limit) {
    // odd-only sieve for base primes <= limit
    vector<int> primes;
    if(limit < 2) return primes;
    primes.push_back(2);
    if(limit < 3) return primes;

    int n_odds = (limit - 1) / 2; // z indices for odds up to limit
    vector<uint8_t> comp(n_odds, 0);

    int r = (int)std::sqrt((long double)limit);
    int zmax = (r - 1) / 2;

    for(int z=1; z<=zmax; z++){
        if(!comp[z]){
            int p = 2*z + 1;
            int start = (p*p - 1) / 2;
            for(int i=start; i<n_odds; i += p) comp[i] = 1;
        }
    }
    for(int z=1; z<n_odds; z++){
        if(!comp[z]) primes.push_back(2*z + 1);
    }
    return primes;
}

static inline int popcount_u64(uint64_t x) {
    return __builtin_popcountll(x);
}

int64_t pi_wheel105_compressed(int64_t N, int64_t seg_mb = 16) {
    if(N < 2) return 0;

    Wheel105 W;

    // base primes <= sqrt(N)
    int64_t lim = isqrt_i64(N);
    while((lim+1)*(lim+1) <= N) lim++;
    while(lim*lim > N) lim--;

    auto primes = base_primes_upto((int)lim);

    // seg_bits capacity from MB
    int64_t seg_bits_cap = seg_mb * 1024LL * 1024LL * 8LL;
    if(seg_bits_cap < (1<<20)) seg_bits_cap = (1<<20);

    // density = 48/105 => seg_bits ≈ seg_z_span*(48/105)
    int64_t seg_z_span = (seg_bits_cap * Wheel105::M) / Wheel105::K;
    if(seg_z_span < (1<<20)) seg_z_span = (1<<20);

    // Count small primes explicitly (2,3,5,7)
    int64_t count = 0;
    if(N >= 2) count++;
    if(N >= 3) count++;
    if(N >= 5) count++;
    if(N >= 7) count++;
    if(N < 11) return count;

    int64_t z_start = (11 - 1) / 2;
    int64_t z_end   = (N  - 1) / 2;

    vector<uint64_t> compw; // reused buffer

    for(int64_t zL = z_start; zL <= z_end; ){
        int64_t zR = zL + seg_z_span - 1;
        if(zR > z_end) zR = z_end;

        int64_t low_n  = 2*zL + 1;
        int64_t high_n = 2*zR + 1;

        int64_t base_rank = W.count_valid_before_z(zL);
        int64_t end_rank  = W.count_valid_before_z(zR + 1);
        int64_t seg_bits  = end_rank - base_rank;

        if(seg_bits <= 0){
            zL = zR + 1;
            continue;
        }

        int64_t nwords = (seg_bits + 63) >> 6;
        compw.assign((size_t)nwords, 0ULL);

        // mark composites for p>7
        for(size_t k=0; k<primes.size(); k++){
            int p = primes[k];
            if(p <= 7) continue;

            int64_t pp = 1LL*p*p;
            if(pp > high_n) break;

            int64_t start = pp;
            if(start < low_n) start = ((low_n + p - 1) / p) * 1LL*p;
            if((start & 1LL) == 0) start += p;
            if(start > high_n) continue;

            int64_t z = (start - 1) / 2;
            if(z < zL) z = zL;

            int a = p % Wheel105::M;
            int res = (int)(z % Wheel105::M);

            if(!W.valid[res]){
                int t0 = W.next0[a][res];
                z += 1LL*p*t0;
                res = (res + a * t0) % Wheel105::M;
            }

            while(z <= zR){
                int64_t gb = W.global_bit_index(z);
                int64_t loc = gb - base_rank;
                compw[(size_t)(loc >> 6)] |= (1ULL << (loc & 63));

                int t = W.next1[a][res];
                z += 1LL*p*t;
                res = (res + a * t) % Wheel105::M;
            }
        }

        // count primes in segment = number of 0 bits among seg_bits
        int64_t full_words = seg_bits >> 6;
        for(int64_t w=0; w<full_words; w++){
            count += popcount_u64(~compw[(size_t)w]);
        }
        int rem = (int)(seg_bits & 63);
        if(rem){
            uint64_t last = ~compw[(size_t)full_words];
            last &= ((rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL));
            count += popcount_u64(last);
        }

        zL = zR + 1;
    }

    return count;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int64_t N = 100000000;
    int64_t seg_mb = 16;

    if(argc >= 2) N = atoll(argv[1]);
    if(argc >= 3) seg_mb = atoll(argv[2]);

    auto t0 = chrono::high_resolution_clock::now();
    int64_t ans = pi_wheel105_compressed(N, seg_mb);
    auto t1 = chrono::high_resolution_clock::now();

    double sec = chrono::duration<double>(t1 - t0).count();
    cout << "N=" << N << " seg_mb=" << seg_mb << "\n";
    cout << "pi(N)=" << ans << "\n";
    cout << "time_s=" << fixed << setprecision(6) << sec << "\n";
    return 0;
}
