/*
 * A052075 — Primes p such that nextprime(p) is substring of p^3
 * ==============================================================
 * Version 4.0 — Post multi-LLM review (Round 3) — Production release
 *
 * High-performance search using:
 *   - Segmented sieve of Eratosthenes for prime generation
 *   - __uint128_t for p^3 when p < 2^42 (~6.98e12)
 *   - GMP (mpz_t) for p^3 when p >= 2^42
 *   - Trailing-match modular filter with TRUE early exit
 *   - OpenMP macro-block parallelism (8+ cores)
 *
 * Compile (generic):
 *   gcc -O2 -march=native -Wall -Wextra -Wpedantic -fopenmp \
 *       -o a052075 a052075_v4.c -lgmp -lm
 *
 * Compile (Ryzen 9 7940HS / Zen 4, optimized):
 *   gcc -O3 -march=znver4 -mtune=znver4 -flto -fopenmp \
 *       -Wall -Wextra -Wpedantic \
 *       -o a052075 a052075_v4.c -lgmp -lm
 *
 * Usage:
 *   ./a052075                    # search from 2 to 10^13
 *   ./a052075 START END          # search in [START, END]
 *   ./a052075 --verify           # verify oracle set only
 *
 * Changelog v3 -> v4 (review Round 3 adjudication):
 *   [R3-05] Fixed 1-byte memset overflow in sieve_segment (Deepseek)
 *   [R3-01] Made pre_seg allocation failure fatal (ChatGPT, Grok)
 *   [R3-02] Dynamic wave_ends allocation (ChatGPT, Grok)
 *   [R3-04] 8x blocks per wave for dynamic load balancing (Gemini)
 *   [R3-03] Extended last block beyond end for final pair (ChatGPT)
 *
 * Author: Carlo Corti, 2026
 * License: MIT
 */

#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <gmp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ──────────────────────────────────────────────────────────
 *  Configuration
 * ────────────────────────────────────────────────────────── */

#define SIEVE_SEGMENT_BYTES  (256 * 1024)
#define SIEVE_SEGMENT_RANGE  ((uint64_t)SIEVE_SEGMENT_BYTES * 16ULL)

#define MAX_LIMIT  10000000000000000ULL  /* 10^16 */
#define P3_128_LIMIT  6980000000000ULL   /* ~2^(128/3) [O8] */

#define PROGRESS_INTERVAL  30
#define ORACLE_SIZE 15
#define P3_STR_BUF  64
#define NP_STR_BUF  24

/* Size of each macro-block for OpenMP parallelism.
 * Each thread processes one macro-block at a time.
 * 10^9 numbers per block ≈ 50M primes ≈ a few seconds of work. */
#define MACRO_BLOCK_SIZE  1000000000ULL  /* 10^9 */

/* ──────────────────────────────────────────────────────────
 *  Powers of 10 lookup table [O6, B4]
 * ────────────────────────────────────────────────────────── */

static const uint64_t POW10[20] = {
    1ULL, 10ULL, 100ULL, 1000ULL, 10000ULL, 100000ULL,
    1000000ULL, 10000000ULL, 100000000ULL, 1000000000ULL,
    10000000000ULL, 100000000000ULL, 1000000000000ULL,
    10000000000000ULL, 100000000000000ULL, 1000000000000000ULL,
    10000000000000000ULL, 100000000000000000ULL,
    1000000000000000000ULL, 10000000000000000000ULL
};

static inline int count_digits(uint64_t n)
{
    if (n == 0) return 1;
    if (n < POW10[8]) {
        if (n < POW10[4]) {
            if (n < POW10[2]) return (n < POW10[1]) ? 1 : 2;
            return (n < POW10[3]) ? 3 : 4;
        }
        if (n < POW10[6]) return (n < POW10[5]) ? 5 : 6;
        return (n < POW10[7]) ? 7 : 8;
    }
    if (n < POW10[16]) {
        if (n < POW10[12]) {
            if (n < POW10[10]) return (n < POW10[9]) ? 9 : 10;
            return (n < POW10[11]) ? 11 : 12;
        }
        if (n < POW10[14]) return (n < POW10[13]) ? 13 : 14;
        return (n < POW10[15]) ? 15 : 16;
    }
    if (n < POW10[18]) return (n < POW10[17]) ? 17 : 18;
    return (n < POW10[19]) ? 19 : 20;
}

/* ──────────────────────────────────────────────────────────
 *  Integer square root [B3, R2-03]
 * ────────────────────────────────────────────────────────── */

static uint64_t isqrt_u64(uint64_t n)
{
    if (n <= 1) return n;
    uint64_t x = (uint64_t)sqrtl((long double)n);
    while (x > 0 && x > n / x) x--;
    while ((x + 1) <= n / (x + 1)) x++;
    return x;
}

/* ──────────────────────────────────────────────────────────
 *  Oracle set (frozen)
 * ────────────────────────────────────────────────────────── */

typedef struct {
    int index;
    uint64_t p;
    uint64_t nextprime;
    const char *p_cubed_str;
    int match_position;
} oracle_entry_t;

static const oracle_entry_t ORACLE[ORACLE_SIZE] = {
    { 1,              11,              13, "1331",                                        0},
    { 2,             101,             103, "1030301",                                     0},
    { 3,            2239,            2243, "11224377919",                                 2},
    { 4,           34297,           34301, "40343019516073",                              2},
    { 5,           43789,           43793, "83964379378069",                              4},
    { 6,           53549,           53551, "153551511228149",                             1},
    { 7,          535487,          535489, "153548930496746303",                          1},
    { 8,        59897017,        59897039, "214889691497505989703913",                   14},
    { 9,       430784719,       430784737, "79943078473759892945966959",                  3},
    {10,      2549592677ULL,   2549592733ULL, "16573430415736921632549592733",            19},
    {11,      2837138669ULL,   2837138677ULL, "22837138677705447754568672309",             1},
    {12,     97969345967ULL,  97969346063ULL, "940309072235302647342697969346063",         22},
    {13,    100000000019ULL, 100000000057ULL, "1000000000570000000108300000006859",         0},
    {14,    328096840219ULL, 328096840223ULL, "35318816603250425999328096840223459",       20},
    {15,   4110739763869ULL,4110739763909ULL, "69464026243759115461642614110739763909",    25},
};

/* ──────────────────────────────────────────────────────────
 *  Small primes sieve
 * ────────────────────────────────────────────────────────── */

static uint32_t *small_primes = NULL;
static uint64_t  small_primes_count = 0;

static void generate_small_primes(uint64_t limit)
{
    uint64_t sieve_size = limit / 2 + 1;
    uint8_t *sieve = calloc(sieve_size, 1);
    if (!sieve) {
        fprintf(stderr, "ERROR: cannot allocate small prime sieve\n");
        exit(1);
    }

    for (uint64_t i = 1; (2 * i + 1) * (2 * i + 1) <= limit; i++) {
        if (!sieve[i]) {
            uint64_t p = 2 * i + 1;
            for (uint64_t j = 2 * i * (i + 1); j < sieve_size; j += p)
                sieve[j] = 1;
        }
    }

    uint64_t count = 1;
    for (uint64_t i = 1; i < sieve_size; i++)
        if (!sieve[i]) count++;

    small_primes = malloc(count * sizeof(uint32_t));
    if (!small_primes) { fprintf(stderr, "ERROR: malloc small_primes\n"); exit(1); }
    small_primes[0] = 2;
    uint64_t idx = 1;
    for (uint64_t i = 1; i < sieve_size; i++)
        if (!sieve[i])
            small_primes[idx++] = (uint32_t)(2 * i + 1);
    small_primes_count = idx;
    free(sieve);
}

/* ──────────────────────────────────────────────────────────
 *  Segmented sieve (unchanged, thread-safe: no global state)
 * ────────────────────────────────────────────────────────── */

typedef void (*prime_pair_callback_t)(uint64_t p, uint64_t q, void *user_data);

static uint64_t sieve_segment(uint64_t lo, uint64_t hi,
                              uint8_t *segment,
                              uint64_t prev_prime,
                              prime_pair_callback_t callback,
                              void *user_data)
{
    uint64_t range = hi - lo;
    /* [R3-05] Cap memset size at allocation size to prevent overflow.
     * Formula gives (range+15)/16+1 which can exceed SIEVE_SEGMENT_BYTES by 1. */
    uint64_t seg_bytes = (range + 15) / 16 + 1;
    if (seg_bytes > SIEVE_SEGMENT_BYTES) seg_bytes = SIEVE_SEGMENT_BYTES;
    memset(segment, 0, seg_bytes);

    for (uint64_t i = 1; i < small_primes_count; i++) {
        uint64_t p = small_primes[i];
        if ((uint64_t)p * p >= hi) break;
        uint64_t start;
        if ((uint64_t)p * p >= lo) {
            start = (uint64_t)p * p;
        } else {
            start = lo + ((p - (lo % p)) % p);
            if (start % 2 == 0) start += p;
        }
        for (uint64_t j = start; j < hi; j += 2 * p) {
            uint64_t b = (j - lo) / 2;
            segment[b / 8] |= (1 << (b % 8));
        }
    }

    uint64_t last_prime = prev_prime;
    for (uint64_t i = (lo % 2 == 0) ? 1 : 0; i < range; i += 2) {
        uint64_t n = lo + i;
        uint64_t b = i / 2;
        if (!(segment[b / 8] & (1 << (b % 8)))) {
            if (last_prime > 0 && callback)
                callback(last_prime, n, user_data);
            last_prime = n;
        }
    }
    return last_prime;
}

/* ──────────────────────────────────────────────────────────
 *  p^3 computation
 * ────────────────────────────────────────────────────────── */

static int cube_to_str_128(uint64_t p, char *buf)
{
    __uint128_t p128 = (__uint128_t)p;
    __uint128_t cube = p128 * p128 * p128;
    char tmp[48];
    int len = 0;
    if (cube == 0) { buf[0] = '0'; buf[1] = '\0'; return 1; }
    while (cube > 0) { tmp[len++] = '0' + (int)(cube % 10); cube /= 10; }
    for (int i = 0; i < len; i++) buf[i] = tmp[len - 1 - i];
    buf[len] = '\0';
    return len;
}

static int cube_to_str_gmp(uint64_t p, char *buf, int bufsize,
                           mpz_t tmp_p, mpz_t tmp_cube)
{
    mpz_set_ui(tmp_p, p);
    mpz_pow_ui(tmp_cube, tmp_p, 3);
    size_t needed = mpz_sizeinbase(tmp_cube, 10) + 2;
    if ((int)needed > bufsize) {
        fprintf(stderr, "ERROR: p^3 overflow for p=%"PRIu64"\n", p);
        exit(1);
    }
    mpz_get_str(buf, 10, tmp_cube);
    return (int)strlen(buf);
}

/* ──────────────────────────────────────────────────────────
 *  Trailing match modular filter
 * ────────────────────────────────────────────────────────── */

static int trailing_match_filter(uint64_t p, uint64_t g)
{
    uint64_t q = p + g;
    int d = count_digits(q);
    if (d > 19) return 0;
    uint64_t mod = POW10[d];
    __uint128_t p128 = (__uint128_t)p;
    __uint128_t p2_mod = (p128 * p128) % mod;
    __uint128_t p2m1_mod = (p2_mod + mod - 1) % mod;
    __uint128_t result = (p128 % mod) * p2m1_mod % mod;
    return (result == (__uint128_t)g);
}

/* ──────────────────────────────────────────────────────────
 *  Thread-local search state [P1]
 * ────────────────────────────────────────────────────────── */

typedef struct {
    uint64_t p, q, gap;
    int match_pos, digits_p3;
    char match_type[12];
} hit_record_t;

typedef struct {
    /* Counters (thread-local, merged after parallel region) */
    uint64_t primes_tested;
    uint64_t trailing_filter_hits;
    uint64_t full_check_count;
    uint64_t primes_per_decade[20];
    uint64_t hits_per_decade[20];

    /* Hit records (thread-local) */
    hit_record_t *hits;
    int hits_count;
    int hits_capacity;

    /* GMP temporaries (thread-local, NOT shared) */
    mpz_t tmp_p, tmp_cube;

    /* String buffers (thread-local) */
    char p3_buf[P3_STR_BUF];
    char np_buf[NP_STR_BUF];
} thread_state_t;

static thread_state_t *create_thread_state(void)
{
    thread_state_t *ts = calloc(1, sizeof(thread_state_t));
    if (!ts) { perror("calloc"); exit(1); }
    ts->hits_capacity = 64;
    ts->hits = malloc(ts->hits_capacity * sizeof(hit_record_t));
    if (!ts->hits) { perror("malloc"); exit(1); }
    mpz_init(ts->tmp_p);
    mpz_init(ts->tmp_cube);
    return ts;
}

static void destroy_thread_state(thread_state_t *ts)
{
    mpz_clear(ts->tmp_p);
    mpz_clear(ts->tmp_cube);
    free(ts->hits);
    free(ts);
}

static void record_hit_local(thread_state_t *ts, uint64_t p, uint64_t q,
                             int match_pos, int digits_p3, const char *match_type)
{
    if (ts->hits_count >= ts->hits_capacity) {
        ts->hits_capacity *= 2;
        ts->hits = realloc(ts->hits, ts->hits_capacity * sizeof(hit_record_t));
        if (!ts->hits) { perror("realloc"); exit(1); }
    }
    hit_record_t *h = &ts->hits[ts->hits_count++];
    h->p = p; h->q = q; h->gap = q - p;
    h->match_pos = match_pos; h->digits_p3 = digits_p3;
    strncpy(h->match_type, match_type, sizeof(h->match_type) - 1);
    h->match_type[sizeof(h->match_type) - 1] = '\0';
}

/* ──────────────────────────────────────────────────────────
 *  Per-pair check [R2-01 — true early exit]
 * ────────────────────────────────────────────────────────── */

static void check_pair(uint64_t p, uint64_t q, void *user_data)
{
    thread_state_t *ts = (thread_state_t *)user_data;
    ts->primes_tested++;

    int d = count_digits(p);
    if (d < 20) ts->primes_per_decade[d]++;

    uint64_t gap = q - p;

    /* Stage 1: Trailing match filter (true early exit) */
    if (trailing_match_filter(p, gap)) {
        ts->trailing_filter_hits++;
        int digits_p3;
        if (p < P3_128_LIMIT) {
            __uint128_t cube = (__uint128_t)p * p * p;
            digits_p3 = 0;
            __uint128_t tmp = cube;
            do { digits_p3++; tmp /= 10; } while (tmp > 0);
        } else {
            mpz_set_ui(ts->tmp_p, p);
            mpz_pow_ui(ts->tmp_cube, ts->tmp_p, 3);
            digits_p3 = (int)mpz_sizeinbase(ts->tmp_cube, 10);
        }
        int np_digits = count_digits(q);
        if (d < 20) ts->hits_per_decade[d]++;
        record_hit_local(ts, p, q, digits_p3 - np_digits, digits_p3, "TRAILING");
        return;
    }

    /* Stage 2: Full substring check */
    ts->full_check_count++;

    int p3_len;
    if (p < P3_128_LIMIT) {
        p3_len = cube_to_str_128(p, ts->p3_buf);
    } else {
        p3_len = cube_to_str_gmp(p, ts->p3_buf, P3_STR_BUF,
                                 ts->tmp_p, ts->tmp_cube);
    }

    int np_len = snprintf(ts->np_buf, NP_STR_BUF, "%"PRIu64, q);
    const char *pos = strstr(ts->p3_buf, ts->np_buf);

    if (pos) {
        int match_pos = (int)(pos - ts->p3_buf);
        const char *match_type;
        if (match_pos == 0) match_type = "LEADING";
        else if (match_pos == p3_len - np_len) match_type = "TRAILING";
        else match_type = "interior";
        if (d < 20) ts->hits_per_decade[d]++;
        record_hit_local(ts, p, q, match_pos, p3_len, match_type);
    }
}

/* ──────────────────────────────────────────────────────────
 *  Process one macro-block [P1]
 *
 *  Each thread calls this independently. The function:
 *  1. Finds the last prime just before block_lo (overlap sieve)
 *  2. Sieves [block_lo, block_hi) in segments
 *  3. Stores results in thread-local state
 * ────────────────────────────────────────────────────────── */

static void process_block(uint64_t block_lo, uint64_t block_hi,
                          thread_state_t *ts)
{
    if (block_lo >= block_hi) return;
    if (block_lo < 3) block_lo = 3;
    if (block_lo % 2 == 0) block_lo++;

    uint64_t seg_range = SIEVE_SEGMENT_RANGE;

    /* Allocate thread-local sieve segment */
    uint8_t *segment = NULL;
    if (posix_memalign((void **)&segment, 64, SIEVE_SEGMENT_BYTES) != 0) {
        fprintf(stderr, "ERROR: cannot allocate sieve segment (thread)\n");
        exit(1);
    }

    /* Find prev_prime: sieve a small range before block_lo.
     * [R3-01] Allocation failure is fatal — if pre_seg can't be allocated,
     * prev_prime stays 0 and we'd lose the cross-boundary pair. */
    uint64_t prev_prime = 0;
    if (block_lo > 3) {
        uint64_t pre_lo = (block_lo > seg_range) ? block_lo - seg_range : 3;
        if (pre_lo % 2 == 0) pre_lo--;
        uint8_t *pre_seg = NULL;
        if (posix_memalign((void **)&pre_seg, 64, SIEVE_SEGMENT_BYTES) != 0) {
            fprintf(stderr, "ERROR: cannot allocate pre-sieve segment\n");
            free(segment);
            exit(1);
        }
        prev_prime = sieve_segment(pre_lo, block_lo, pre_seg, 0, NULL, NULL);
        free(pre_seg);
    }

    /* Sieve block_lo..block_hi in segments */
    uint64_t seg_lo = block_lo;
    while (seg_lo < block_hi) {
        uint64_t seg_hi = seg_lo + seg_range;
        if (seg_hi > block_hi) seg_hi = block_hi;
        prev_prime = sieve_segment(seg_lo, seg_hi, segment,
                                   prev_prime, check_pair, ts);
        seg_lo = seg_hi;
    }

    free(segment);
}

/* ──────────────────────────────────────────────────────────
 *  Checkpoint support
 * ────────────────────────────────────────────────────────── */

static void write_checkpoint(const char *filename, uint64_t last_block_hi,
                             uint64_t total_primes, uint64_t total_hits,
                             uint64_t trailing_hits, uint64_t full_checks,
                             double elapsed)
{
    FILE *f = fopen(filename, "w");
    if (!f) return;
    fprintf(f, "CHECKPOINT A052075 v4\n");
    fprintf(f, "last_block_hi=%"PRIu64"\n", last_block_hi);
    fprintf(f, "primes_tested=%"PRIu64"\n", total_primes);
    fprintf(f, "hits_found=%"PRIu64"\n", total_hits);
    fprintf(f, "trailing_filter_hits=%"PRIu64"\n", trailing_hits);
    fprintf(f, "full_check_count=%"PRIu64"\n", full_checks);
    fprintf(f, "elapsed=%.1f\n", elapsed);
    fclose(f);
}

static uint64_t read_checkpoint(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if (!f) return 0;
    char line[256];
    uint64_t last_hi = 0;
    while (fgets(line, sizeof(line), f))
        sscanf(line, "last_block_hi=%"SCNu64, &last_hi);
    fclose(f);
    if (last_hi > 0)
        printf("Checkpoint loaded: resuming from %"PRIu64"\n", last_hi);
    return last_hi;
}

/* ──────────────────────────────────────────────────────────
 *  Oracle verification
 * ────────────────────────────────────────────────────────── */

static int verify_oracle(void)
{
    printf("======================================================\n");
    printf("ORACLE SET VERIFICATION — A052075 (v4)\n");
    printf("======================================================\n");

    mpz_t tmp_p, tmp_cube;
    mpz_init(tmp_p);
    mpz_init(tmp_cube);
    int all_ok = 1;

    for (int i = 0; i < ORACLE_SIZE; i++) {
        const oracle_entry_t *o = &ORACLE[i];
        char p3_buf[128];
        int p3_len;
        if (o->p < P3_128_LIMIT)
            p3_len = cube_to_str_128(o->p, p3_buf);
        else
            p3_len = cube_to_str_gmp(o->p, p3_buf, sizeof(p3_buf), tmp_p, tmp_cube);

        if (strcmp(p3_buf, o->p_cubed_str) != 0) {
            printf("  FAIL a(%d): p^3 mismatch\n", o->index);
            all_ok = 0; continue;
        }

        char np_buf[24];
        int np_len = snprintf(np_buf, sizeof(np_buf), "%"PRIu64, o->nextprime);
        const char *match = strstr(p3_buf, np_buf);
        int pos = match ? (int)(match - p3_buf) : -1;

        if (pos != o->match_position) {
            printf("  FAIL a(%d): match at %d, expected %d\n", o->index, pos, o->match_position);
            all_ok = 0; continue;
        }

        uint64_t gap = o->nextprime - o->p;
        int trailing = trailing_match_filter(o->p, gap);
        int is_trailing = (pos == p3_len - np_len);
        if (is_trailing && !trailing) {
            printf("  FAIL a(%d): trailing filter missed!\n", o->index);
            all_ok = 0; continue;
        }

        printf("  OK   a(%2d) = %15"PRIu64"  |  p^3 OK  |  match@%d/%d",
               o->index, o->p, pos, p3_len);
        if (is_trailing) printf("  [TRAILING, filter=OK]");
        else if (trailing) printf("  [filter=false positive]");
        printf("\n");
    }

    printf("------------------------------------------------------\n");
    printf(all_ok ? "ALL %d TERMS VERIFIED SUCCESSFULLY.\n"
                  : "*** VERIFICATION FAILED ***\n", ORACLE_SIZE);
    printf("\n");
    mpz_clear(tmp_p);
    mpz_clear(tmp_cube);
    return all_ok;
}

/* ──────────────────────────────────────────────────────────
 *  Main search driver [P1 — OpenMP parallel]
 * ────────────────────────────────────────────────────────── */

static void run_search(uint64_t start, uint64_t end)
{
    printf("======================================================\n");
    printf("A052075 SEARCH: [%"PRIu64", %"PRIu64"]\n", start, end);
    printf("======================================================\n\n");

    uint64_t sqrt_end = isqrt_u64(end) + 1;
    printf("Generating small primes up to %"PRIu64"...\n", sqrt_end);
    generate_small_primes(sqrt_end);
    printf("Found %"PRIu64" small primes (largest: %u)\n\n",
           small_primes_count, small_primes[small_primes_count - 1]);

    /* Open log file */
    char logname[256];
    snprintf(logname, sizeof(logname),
             "a052075_log_%"PRIu64"_%"PRIu64".csv", start, end);
    FILE *logfile = fopen(logname, "a");
    if (logfile) {
        fprintf(logfile, "# A052075 search log (v3, OpenMP)\n");
        fprintf(logfile, "# Range: [%"PRIu64", %"PRIu64"]\n", start, end);
        fflush(logfile);
    }

    /* Check for checkpoint */
    uint64_t resume_from = read_checkpoint("a052075_checkpoint.txt");
    if (resume_from > start) {
        start = resume_from;
        printf("Resuming from %"PRIu64"\n\n", start);
    }

    if (start < 3) start = 3;
    if (start % 2 == 0) start++;

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
#endif
    printf("Using %d threads (macro-block size: %"PRIu64")\n\n",
           num_threads, (uint64_t)MACRO_BLOCK_SIZE);

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    double start_time = t0.tv_sec + t0.tv_nsec * 1e-9;
    double last_progress = start_time;

    /* Global accumulators */
    uint64_t total_primes = 0;
    uint64_t total_trailing = 0;
    uint64_t total_full = 0;
    uint64_t total_hits = 0;
    uint64_t global_primes_per_decade[20] = {0};
    uint64_t global_hits_per_decade[20] = {0};

    /* Collect ALL hits from all threads, then sort at the end */
    hit_record_t *all_hits = malloc(256 * sizeof(hit_record_t));
    int all_hits_count = 0;
    int all_hits_capacity = 256;

    /* Process in waves of macro-blocks.
     * [R3-04] Dispatch up to 8× num_threads blocks per wave to enable
     * true dynamic load balancing (schedule(dynamic,1)).
     * [R3-03] Extend the final block by one segment beyond `end` to
     * ensure the last prime in range gets its nextprime tested.
     * After each wave, merge results, report progress, checkpoint. */
    int max_blocks_per_wave = num_threads * 8;
    uint64_t wave_lo = start;

    /* [R3-03] Effective search end: extend by one sieve segment so that
     * the last prime before `end` can be paired with its nextprime. */
    uint64_t search_end = end + SIEVE_SEGMENT_RANGE;

    while (wave_lo < search_end) {
        /* [R3-02] Dynamic allocation for wave block boundaries */
        int blocks_in_wave = 0;
        uint64_t *wave_ends = malloc((size_t)max_blocks_per_wave * sizeof(uint64_t));
        if (!wave_ends) { perror("malloc wave_ends"); exit(1); }

        uint64_t block_lo = wave_lo;
        while (block_lo < search_end && blocks_in_wave < max_blocks_per_wave) {
            uint64_t block_hi = block_lo + MACRO_BLOCK_SIZE;
            if (block_hi > search_end) block_hi = search_end;
            wave_ends[blocks_in_wave] = block_hi;
            blocks_in_wave++;
            block_lo = block_hi;
        }
        uint64_t wave_hi = wave_ends[blocks_in_wave - 1];

        /* Allocate per-thread state */
        thread_state_t **thread_states = malloc(blocks_in_wave * sizeof(thread_state_t *));
        for (int i = 0; i < blocks_in_wave; i++)
            thread_states[i] = create_thread_state();

        /* ---- PARALLEL REGION ---- */
        #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (int b = 0; b < blocks_in_wave; b++) {
            uint64_t blo = (b == 0) ? wave_lo : wave_ends[b - 1];
            uint64_t bhi = wave_ends[b];
            process_block(blo, bhi, thread_states[b]);
        }

        /* ---- MERGE RESULTS (sequential) ---- */
        for (int b = 0; b < blocks_in_wave; b++) {
            thread_state_t *ts = thread_states[b];
            total_primes += ts->primes_tested;
            total_trailing += ts->trailing_filter_hits;
            total_full += ts->full_check_count;

            for (int d = 0; d < 20; d++) {
                global_primes_per_decade[d] += ts->primes_per_decade[d];
                global_hits_per_decade[d] += ts->hits_per_decade[d];
            }

            /* Merge hits */
            for (int h = 0; h < ts->hits_count; h++) {
                if (all_hits_count >= all_hits_capacity) {
                    all_hits_capacity *= 2;
                    all_hits = realloc(all_hits, all_hits_capacity * sizeof(hit_record_t));
                    if (!all_hits) { perror("realloc"); exit(1); }
                }
                all_hits[all_hits_count++] = ts->hits[h];
                total_hits++;

                /* Print immediately */
                hit_record_t *r = &all_hits[all_hits_count - 1];
                printf("  HIT #%"PRIu64": p = %"PRIu64", nextprime = %"PRIu64
                       ", gap = %"PRIu64", pos = %d/%d, type = %s\n",
                       total_hits, r->p, r->q, r->gap,
                       r->match_pos, r->digits_p3, r->match_type);
                fflush(stdout);

                if (logfile) {
                    fprintf(logfile, "HIT,%"PRIu64",%"PRIu64",%"PRIu64
                            ",%d,%d,%s\n",
                            r->p, r->q, r->gap,
                            r->match_pos, r->digits_p3, r->match_type);
                    fflush(logfile);
                }
            }

            destroy_thread_state(ts);
        }
        free(thread_states);
        free(wave_ends);  /* [R3-02] */

        /* Progress reporting */
        struct timespec tnow;
        clock_gettime(CLOCK_MONOTONIC, &tnow);
        double now = tnow.tv_sec + tnow.tv_nsec * 1e-9;
        double elapsed = now - start_time;

        if (now - last_progress > PROGRESS_INTERVAL || wave_hi >= search_end) {
            double rate = (double)total_primes / elapsed;
            /* Show progress relative to original end, not search_end */
            double pct_raw = 100.0 * (double)(wave_hi - start) / (double)(end - start);
            double pct = (pct_raw > 100.0) ? 100.0 : pct_raw;
            printf("  ... wave done: [%"PRIu64", %"PRIu64"), "
                   "total %"PRIu64" primes, %.0f/sec, "
                   "hits: %"PRIu64", %.1f%%, %.0fs elapsed\n",
                   wave_lo, wave_hi, total_primes, rate,
                   total_hits, pct, elapsed);
            fflush(stdout);
            last_progress = now;
        }

        /* Checkpoint */
        write_checkpoint("a052075_checkpoint.txt", wave_hi,
                         total_primes, total_hits, total_trailing,
                         total_full, elapsed);

        wave_lo = wave_hi;
    }

    /* Final statistics */
    struct timespec tend;
    clock_gettime(CLOCK_MONOTONIC, &tend);
    double total_time = tend.tv_sec + tend.tv_nsec * 1e-9 - start_time;

    printf("\n======================================================\n");
    printf("SEARCH COMPLETE\n");
    printf("======================================================\n\n");
    printf("Range: [%"PRIu64", %"PRIu64"]\n", start, end);
    printf("Threads: %d\n", num_threads);
    printf("Primes tested: %"PRIu64"\n", total_primes);
    printf("Hits found: %"PRIu64"\n", total_hits);
    printf("Trailing filter hits: %"PRIu64"\n", total_trailing);
    printf("Full checks: %"PRIu64"\n", total_full);
    printf("Time: %.1f seconds (%.2f hours)\n", total_time, total_time / 3600);
    printf("Rate: %.0f primes/sec (%.0f/sec/thread)\n",
           (double)total_primes / total_time,
           (double)total_primes / total_time / num_threads);

    printf("\nPer-decade statistics:\n");
    printf("  Digits  Primes tested    Hits\n");
    for (int d = 2; d < 20; d++)
        if (global_primes_per_decade[d] > 0)
            printf("  %5d   %13"PRIu64"   %4"PRIu64"\n",
                   d, global_primes_per_decade[d], global_hits_per_decade[d]);

    /* Sort hits by p */
    for (int i = 0; i < all_hits_count - 1; i++)
        for (int j = i + 1; j < all_hits_count; j++)
            if (all_hits[j].p < all_hits[i].p) {
                hit_record_t tmp = all_hits[i];
                all_hits[i] = all_hits[j];
                all_hits[j] = tmp;
            }

    printf("\nHits found (sorted by p):\n");
    for (int i = 0; i < all_hits_count; i++) {
        hit_record_t *h = &all_hits[i];
        printf("  p = %"PRIu64", nextprime = %"PRIu64
               ", gap = %"PRIu64", pos = %d/%d, type = %s\n",
               h->p, h->q, h->gap, h->match_pos, h->digits_p3, h->match_type);
    }

    if (logfile) {
        fprintf(logfile, "# COMPLETE: %"PRIu64" primes, %"PRIu64
                " hits, %.1fs, %d threads\n",
                total_primes, total_hits, total_time, num_threads);
        fclose(logfile);
    }

    free(all_hits);
    free(small_primes);
}

/* ──────────────────────────────────────────────────────────
 *  Main
 * ────────────────────────────────────────────────────────── */

int main(int argc, char *argv[])
{
    printf("A052075 — Primes p such that nextprime(p) is substring of p^3\n");
    printf("Version 4.0 — Carlo Corti, 2026\n\n");

    if (argc >= 2 && strcmp(argv[1], "--verify") == 0)
        return verify_oracle() ? 0 : 1;

    uint64_t start = 2;
    uint64_t end = 10000000000000ULL;

    if (argc >= 3) {
        start = strtoull(argv[1], NULL, 10);
        end   = strtoull(argv[2], NULL, 10);
    } else if (argc == 2) {
        end = strtoull(argv[1], NULL, 10);
    }

    if (end > MAX_LIMIT) {
        fprintf(stderr, "ERROR: limit exceeds MAX_LIMIT\n");
        return 1;
    }

    if (!verify_oracle()) {
        fprintf(stderr, "FATAL: oracle verification failed.\n");
        return 1;
    }

    run_search(start, end);
    return 0;
}
