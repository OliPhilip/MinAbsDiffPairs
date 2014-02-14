/// Algorithms and Data Structures
#ifndef cf_al_ds

#define cf_al_ds

#ifndef cf_stdlib
    #define cf_stdlib
    #include <stdlib.h>
#endif
#ifndef cf_chrono
    #define cf_chrono
    #include <chrono>
#endif
#ifndef cf_iostream
    #define cf_iostream
    #include <iostream>
#endif

#define al_mom_nb 3
#define ds_prg_length 624
#define ds_prg_bitMask_32 0xffffffff
#define ds_prg_bitPow_31 0x80000000

namespace algo {
    static uint32_t bitmask[] = {0x1,0x3,0xf,0xff,0xffff};
    // Order statistic (3 ways) indexed with call back
    template<class big_obj> size_t OrderStat(big_obj *container, size_t *index,
                                             size_t start, size_t stop, size_t order,
                                             bool (*is_inferior)(big_obj, big_obj)) {
        big_obj pivot = container[index[start]];
        size_t lower = start, pos = start, upper = stop - 1, exg;
        while (upper >= pos) {
            if (is_inferior(container[index[pos]], pivot)) { exg = index[lower]; index[lower++] = index[pos]; index[pos++] = exg; }
            else if (is_inferior(pivot, container[index[pos]])) { exg = index[pos]; index[pos] = index[upper]; index[upper--] = exg; }
            else ++pos;
        }
        if (order < lower) return OrderStat(container, index, start, lower, order, is_inferior);
        if (order > upper) return OrderStat(container, index, ++upper, stop, order, is_inferior);
        return order;
    }
    // Median of medians indexed with call back
    template<class big_obj> size_t MedianOfMedians(big_obj *container, size_t *index,
                                                   size_t start, size_t stop, size_t nb_median,
                                                   bool (*is_inferior)(big_obj, big_obj)) {
        size_t ext_med[nb_median], nb_ext_med = (stop - start) / nb_median, nb_ext_empty = 0;
        for (size_t i = 0; i < nb_median; ++i) {
            size_t nb_int_med = nb_ext_med;
            if (i == 0) nb_int_med += (stop - start) % nb_median;
            if (nb_int_med == 0) ++nb_ext_empty;
            else {
                size_t int_med[nb_median], nb_int_empty = 0;
                for (size_t j = 0; j < nb_median; ++j) {
                    size_t nb = nb_int_med / nb_median;
                    if (j == 0) nb += nb_int_med % nb_median;
                    if (nb == 0) ++nb_int_empty;
                    else { int_med[j - nb_int_empty] = index[start + (nb / 2)]; start += nb; }
                }
                nb_int_empty = nb_median - nb_int_empty;
                nb_int_empty = OrderStat(container, int_med, (size_t)0, nb_int_empty, nb_int_empty / 2, is_inferior);
                ext_med[i - nb_ext_empty] = int_med[nb_int_empty];
            }
        }
        nb_median -= nb_ext_empty;
        nb_median = OrderStat(container, ext_med, (size_t)0, nb_median, nb_median / 2, is_inferior);
        return ext_med[nb_median];
    }
    // Recursive Quick sort (3 ways) indexed with call back
    template<class big_obj> void QuickSort(big_obj *container, size_t *index,
                                           size_t start, size_t stop,
                                           bool (*is_inferior)(big_obj, big_obj), bool med = true) {
        // Pivot can be selected randomly or by median of medians method
        big_obj pivot;
        if (med) pivot = container[MedianOfMedians(container, index, start, stop, al_mom_nb, is_inferior)];
        else pivot = container[index[start]]; // not recommended if not shuffled
        size_t lower = start, pos = start, upper = stop - 1, exg;
        while (upper >= pos) {
            if (is_inferior(container[index[pos]], pivot)) {
                exg = index[lower]; index[lower] = index[pos]; index[pos] = exg; ++lower; ++pos;
            } else if (is_inferior(pivot, container[index[pos]])) {
                exg = index[pos]; index[pos] = index[upper]; index[upper] = exg; --upper;
            } else ++pos;
        }
        if ((lower - start) > 1) QuickSort(container, index, start, lower, is_inferior);
        ++upper;
        if ((stop - upper) > 1) QuickSort(container, index, upper, stop, is_inferior);
    }
    // Indexed Radix sort LSD stable
    template<class big_obj> void RadixLSD(big_obj *container, size_t *index,
                                          size_t start, size_t stop, bool sign = false) {
        size_t ans[stop - start], *src, *dst;
        bool way = true;
        size_t i, nb[257], rot, byte = sizeof(big_obj) << 3, ind, b, e, o;
        for (rot = 0; rot < byte; rot += 8) {
            if (way) { src = index; dst = ans; b = start; e = stop; o = 0; }
            else { src = ans; dst = index; b = 0; e = stop - start; o = start; }
            way = !way;
            if (sign && (rot + 8) == byte) {
                for (i = 0; i < 257; ++i) nb[i] = 0;
                for (i = b; i < e; ++i) ++nb[(((container[src[i]] >> rot) & 0xff) ^ 0x80) + 1];
                for (i = 2; i < 256; ++i)  nb[i] += nb[i - 1];
                for (i = b; i < e; ++i) { ind = ((container[src[i]] >> rot) & 0xff) ^ 0x80; dst[nb[ind] + o] = src[i]; ++nb[ind]; }
            } else {
                for (i = 0; i < 257; ++i) nb[i] = 0;
                for (i = b; i < e; ++i) ++nb[((container[src[i]] >> rot) & 0xff) + 1];
                for (i = 2; i < 256; ++i)  nb[i] += nb[i - 1];
                for (i = b; i < e; ++i) { ind = (container[src[i]] >> rot) & 0xff; dst[nb[ind] + o] = src[i]; ++nb[ind]; }
            }
        }
    }
    // Mersenne twister pseudo-random number generator
    class RandomGenerator{
        uint32_t *mt, idx = 0;
        void _gen() {
            for(uint32_t i = 0; i < ds_prg_length; ++i) {
                uint32_t y = (mt[i] & ds_prg_bitPow_31) + (mt[(i + 1) % ds_prg_length] & (ds_prg_bitPow_31 - 1));
                mt[i] = mt[(i + 397) % ds_prg_length] ^ (y >> 1);
                if (y % 2) mt[i] ^= 2567483615;
            }
        }
        public:
        RandomGenerator(uint32_t seed) {
            mt = (uint32_t*)malloc(ds_prg_length * sizeof(uint32_t)); mt[0] = seed;
            for (uint32_t i = 1; i < ds_prg_length; ++i) mt[i] = (1812433253 * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i) & ds_prg_bitMask_32;
        }
        ~RandomGenerator() { free(mt); }
        uint32_t get() {
            if (idx == 0) _gen();
            uint32_t y = mt[idx];
            y ^= y >> 11; y ^= (y << 7) & 2636928640; y ^= (y << 15) & 4022730752; y ^= y >> 18; idx = (idx + 1) % ds_prg_length;
            return y;
        }
    };
    // Multiplication modulo
    inline uint32_t MultiModulo(uint32_t a, uint32_t m, uint32_t n) {
        if (!a || !m) return 0;
        uint32_t ans = 0;
        while (m) {
            if (m & 1) ans = (ans + a) % n;
            m = m >> 1; a = (a << 1) % n;
        }
        return ans;
    }
    // Power modulo
    inline uint32_t PowerModulo(uint32_t a, uint32_t m, uint32_t n) {
        if (!a) return 0;
        uint32_t ans = 1;
        while (m) {
            if (m & 1) ans = MultiModulo(ans, a, n);
            m = m >> 1; a = MultiModulo(a, a, n);
        }
        return ans;
    }
    // Probabilistic primality test
    inline bool MillerRabin(uint32_t n, size_t k = 6) {
        RandomGenerator gen = RandomGenerator(std::chrono::system_clock::now().time_since_epoch().count());
        uint32_t a, s, d, r = n - 1;
        for (size_t i = 0; i < k; ++i) {
            s = 1; d = n; a = 1 + (gen.get() % (n - 2));
            while (((d = (d >> 1)) & 1) == 0) s = s << 1;
            a = PowerModulo(a, d, n);
            if (a != 1 && a != r) {
                while ((s = s >> 1) != 0) { if ((a = PowerModulo(a, 2, n)) == r) break; }
                if (s == 0) return false;
            }
        }
        return true;
    }
    // Hash function
    class hash1 {
        algo::RandomGenerator gen = algo::RandomGenerator(std::chrono::system_clock::now().time_since_epoch().count());
        #define al_hash1_nbp 4
        public:
        hash1() {}
        uint32_t hashing(const uint8_t *arg, uint32_t argsize, uint32_t bucket, uint32_t *param, uint32_t os = 0) {
            uint32_t res = 0, i;
            os *= al_hash1_nbp;
            for (i = 0; i < argsize; ++i) res = (i + res + (uint32_t)(arg[i]) * param[os + (i % al_hash1_nbp)]) % bucket;
            return res;
        }
        void newparams(uint32_t *param, uint32_t bucket, uint32_t numhash = 0) {
            uint32_t mod = (bucket > UINT8_MAX ? bucket : UINT8_MAX) - 1, i;
            for (i = 0; i < al_hash1_nbp; ++i) param[(numhash * al_hash1_nbp) + i] = 1 + (gen.get() % mod);
        }
        uint32_t NbParam() { return al_hash1_nbp; }
    };
    // Hash function - mimdiff
    class hash_mindiff {
        algo::RandomGenerator gen = algo::RandomGenerator(std::chrono::system_clock::now().time_since_epoch().count());
        #define al_hashmimdiff_nbp 0
        public:
        uint32_t shift = 0, offset = 0;
        hash_mindiff() {}
        uint32_t hashing(const uint8_t *arg, uint32_t argsize, uint32_t bucket, uint32_t *param, uint32_t os = 0) {
            uint32_t res = *((uint32_t*)arg);
            return ((res >> shift) + offset) % bucket;
        }
        void newparams(uint32_t *param, uint32_t bucket, uint32_t numhash = 0) {
            uint32_t mod = (bucket > UINT8_MAX ? bucket : UINT8_MAX) - 1, i;
            for (i = 0; i < al_hashmimdiff_nbp; ++i) param[(numhash * al_hashmimdiff_nbp) + i] = 1 + (gen.get() % mod);
        }
        uint32_t NbParam() { return al_hashmimdiff_nbp; }
    };
    // Shuffle array
    template<class big_obj> void PermuteArray(big_obj *arg, uint32_t sz) {
        algo::RandomGenerator gen = algo::RandomGenerator(std::chrono::system_clock::now().time_since_epoch().count());
        uint32_t pos;
        big_obj val;
        while (sz > 1) { pos = gen.get() % sz; val = arg[pos]; arg[pos] = arg[--sz]; arg[sz] = val; }
    }
}
namespace data_struct {
    template<class big_obj> class hash_info {
        public:
        big_obj *data;
        hash_info() {}
        hash_info(big_obj *arg) { data = arg; }
        uint32_t GetSize() { return sizeof(big_obj); }
        uint8_t* GetData() { return (uint8_t*)data; }
        bool IsEqual(hash_info *arg) { return (*data == *(arg->data)); }
    };
    // Hash table - open address - hopscotch 32bits
    template<class big_obj, class hfct> class hash_hopscotch {
        private:
        uint32_t bucket_min, nb_param, nb_obj, *param;
        float threshold;
        bool _insert(big_obj *arg, uint32_t key, big_obj **tab, uint32_t *bittab, uint32_t buck) {
            uint32_t pos, temp, diff, ins = key;
            while (tab[ins % buck]) ++ins;
            if (bittab[key] == ds_prg_bitMask_32) { tab[ins % buck] = arg; return false; }
            while ((ins - key) > 31) {
                pos = 31; diff = 0;
                while (!((temp = bittab[(ins - pos) % buck]) & (ds_prg_bitMask_32 >> (32 - pos)))) {
                    if (!(--pos)) { tab[ins % buck] = arg; return false; }
                }
                if (!(temp & algo::bitmask[4])) { diff += 16; temp = temp >> 16; }
                if (!(temp & algo::bitmask[3])) { diff +=  8; temp = temp >>  8; }
                if (!(temp & algo::bitmask[2])) { diff +=  4; temp = temp >>  4; }
                if (!(temp & algo::bitmask[1])) { diff +=  2; temp = temp >>  2; }
                if (!(temp & algo::bitmask[0])) { diff +=  1; }
                tab[ins % buck] = tab[(ins - pos + diff) % buck]; bittab[(ins - pos) % buck] ^= (0x1 << diff) | (0x1 << pos);
                ins -= pos - diff;
            }
            bittab[key] |= 0x1 << (ins - key); tab[ins % buck] = arg;
            return true;
        }
        big_obj * _lookup(big_obj *arg, uint32_t key, bool del = false) {
            big_obj *val;
            uint32_t mod = key, mask = 0x1;
            while (mask) {
                if ((bittable[key] & mask) && arg->IsEqual(table[mod])) {
                    val = table[mod];
                    if (del) { bittable[key] ^= mask; table[mod] = nullptr; }
                    return val;
                }
                mod = (mod + 1) % bucket; mask = mask << 1;
            }
            return nullptr;
        }
        void _rehash(uint32_t newbucket) {
            big_obj **newtable = (big_obj**)calloc((size_t)newbucket, sizeof(big_obj*));
            uint32_t *newbittable = (uint32_t*)calloc((size_t)newbucket, sizeof(uint32_t)), pos = 0;
            fct.newparams(param, newbucket);
            while (pos < bucket) {
                if (table[pos] && !_insert(table[pos], fct.hashing(table[pos]->GetData(),
                                                                   table[pos]->GetSize(),
                                                                   newbucket,
                                                                   param), newtable, newbittable, newbucket)) {
                    for (uint32_t i = 0; i < newbucket; ++i) { newtable[i] = nullptr; newbittable[i] = 0; }
                    fct.newparams(param, newbucket); pos = 0;
                } else ++pos;
            }
            free(table); free(bittable); table = newtable; bittable = newbittable; bucket = newbucket;
        }
        public:
        hfct fct;
        big_obj **table;
        uint32_t *bittable, bucket;
        hash_hopscotch(uint32_t _bucket, float _threshold = 0.9) {
            _bucket = _bucket | 1; while (!algo::MillerRabin(_bucket)) _bucket += 2;
            bucket_min = _bucket; bucket = _bucket; threshold = _threshold; nb_obj = 0; nb_param = fct.NbParam();
            table = (big_obj**)calloc(bucket, sizeof(big_obj*)); bittable = (uint32_t*)calloc(bucket, sizeof(uint32_t));
            param = (uint32_t*)malloc(nb_param * sizeof(uint32_t)); fct.newparams(param, bucket);
        }
        ~hash_hopscotch() { free(table); free(bittable); free(param); }
        bool Insert(big_obj *arg) {
            uint32_t key;
            return Insert(arg, key);
        }
        bool Insert(big_obj *arg, uint32_t &key) {
            key = fct.hashing(arg->GetData(), arg->GetSize(), bucket, param);
            if (_lookup(arg, key) != nullptr) return false;
            bool reh = false;
            uint32_t newbucket = bucket;
            if (!_insert(arg, key, table, bittable, bucket)) {
                if (!nb_param) { std::cout << "ERR " << nb_obj << std::endl; return false; }
                reh = true;
            }
            if (((float)(++nb_obj) / (float)bucket) >= threshold) {
                reh = true; newbucket = (2 * bucket) | 1;
                while (!algo::MillerRabin(newbucket)) newbucket += 2;
            }
            if (reh) _rehash(newbucket);
            return true;
        }
        big_obj * Lookup(big_obj *arg) { return _lookup(arg, fct.hashing(arg->GetData(), arg->GetSize(), bucket, param)); }
        void Erase() {
            for (uint32_t i = 0; i < bucket; ++i) {
                if (table[i]) delete table[i];
                table[i] = nullptr; bittable[i] = 0;
            }
            nb_obj = 0;
        }
    };
}
namespace algo {
    template<class big_obj> bool IsInferior(big_obj obj1, big_obj obj2) { return (obj1 < obj2); }
    // No duplicate random array
    template<class big_obj> big_obj* RandomArray(size_t nb, big_obj mask) {
        typedef data_struct::hash_info<big_obj> HI;
        big_obj *ans = (big_obj*)malloc(nb * sizeof(big_obj));
        data_struct::hash_hopscotch<HI, hash1> ha(nb * 10 / 9 + 1);
        RandomGenerator gen = RandomGenerator(std::chrono::system_clock::now().time_since_epoch().count());
        HI *arg;
        size_t i = 0;
        while (i < nb) {
            arg = new HI(); arg->data = &ans[i]; ans[i] = (big_obj)(gen.get() & mask);
            while (!ha.Insert(arg)) ans[i] = (big_obj)(gen.get() & mask);
            ++i;
        }
        ha.Erase();
        return ans;
    }
    // Find minimum absolute difference - radix sort
    template<class big_obj> void MinDiffPairs_radix(big_obj *data, size_t nb) {
        size_t index[nb], nbmin = 1;
        for (size_t i = 0; i < nb; ++i) index[i] = i;
        RadixLSD(data, index, 0, nb, false);
        big_obj minval = data[index[1]] - data[index[0]], temp;
        for (size_t i = 2; i < nb; ++i) {
            if ((temp = data[index[i]] - data[index[i - 1]]) < minval) { minval = temp; nbmin = 0; }
            if (temp == minval) ++nbmin;
        }
        std::cout << "Diff: " << minval << " - #: " << nbmin << " ";
    }
    // Find minimum absolute difference - quick sort
    template<class big_obj> void MinDiffPairs_quick(big_obj *data, size_t nb) {
        size_t index[nb], nbmin = 1;
        for (size_t i = 0; i < nb; ++i) index[i] = i;
        QuickSort(data, index, 0, nb, IsInferior);
        big_obj minval = data[index[1]] - data[index[0]], temp;
        for (size_t i = 2; i < nb; ++i) {
            if ((temp = data[index[i]] - data[index[i - 1]]) < minval) { minval = temp; nbmin = 0; }
            if (temp == minval) ++nbmin;
        }
        std::cout << "Diff: " << minval << " - #: " << nbmin << " ";
    }
    // Helper function - finding min diff of a newly added number within a hash bucket
    template<class big_obj> big_obj GetMinDiff(data_struct::hash_info<big_obj> **table,
                                               data_struct::hash_info<big_obj> *arg,
                                               uint32_t bits, uint32_t key, uint32_t bucket) {
        uint32_t os = 0;
        big_obj mindiff = (big_obj)0, temp;
        while (bits) {
            if (bits & 0x1) {
                temp = *(table[(key + os) % bucket]->data);
                if (temp < *(arg->data)) temp = *(arg->data) - temp;
                else temp -= *(arg->data);
                if (temp && (!mindiff || temp < mindiff)) mindiff = temp;
            }
            ++os; bits = bits >> 1;
        }
        return mindiff;
    }
    // Find minimum absolute difference - hash table
    template<class big_obj> void MinDiffPairs_hash(big_obj *data, size_t nb) {
        typedef data_struct::hash_info<big_obj> HI;
        HI index[nb];
        size_t i, j, maxi = sizeof(big_obj) << 3;
        uint32_t key, cpt = 0, badcol = 0;
        big_obj maxdiff, mindiff, temp;
        for (i = 0; i < maxi; ++i) {
            data_struct::hash_hopscotch<HI, hash_mindiff> ha(nb * 4);
            ha.fct.shift = i; maxdiff = (big_obj)((0x1 << (i + 1)) - 1); mindiff = (big_obj)0;
            for (j = 0; j < nb; ++j) {
                if (!i) index[j].data = &data[j];
                ha.fct.offset = 0; if (!ha.Insert(&index[j], key)) return;
                temp = GetMinDiff(ha.table, &index[j], ha.bittable[key], key, ha.bucket);
                if (temp) {
                    if (temp > maxdiff) ++badcol;
                    else if (temp == mindiff) ++cpt;
                    else if (!mindiff || temp < mindiff) { mindiff = temp; cpt = 1; }
                }
                ha.fct.offset = 1; if (!ha.Insert(&index[j], key)) return;
                temp = GetMinDiff(ha.table, &index[j], ha.bittable[key], key, ha.bucket);
                if (temp) {
                    if (temp > maxdiff) ++badcol;
                    else if (temp == mindiff) ++cpt;
                    else if (!mindiff || temp < mindiff) { mindiff = temp; cpt = 1; }
                }
            }
            if (mindiff) { std::cout << "Diff: " << mindiff << " - #: " << cpt << " (" << 100 * badcol / (nb * (i + 1)) << "%) "; return; }
        }
    }
}
#endif // cf_al_ds
