#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <stdio.h>
#include <math.h>
#include <queue>
#include <unordered_map>

#if defined(__GNUC__) || defined(__clang__)
#  define SSAM_MAYBE_UNUSED __attribute__((unused))
#else
#  define SSAM_MAYBE_UNUSED
#endif

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
#endif

/*
 * SIMD backend selection.
 *
 *   x86/x86-64 : AVX-512F, detected at runtime with CPUID.
 *   aarch64    : NEON (always present in the base ISA) plus an optional SVE
 *                variant that is compiled with a function-level target
 *                attribute and selected at runtime through AT_HWCAP, so a
 *                single binary runs on SVE and non-SVE cores alike.
 *   otherwise  : portable scalar fallback.
 *
 * Two escape hatches exist for testing and for troubleshooting a build:
 * -DSSAM_NO_SVE drops the SVE kernel, -DSSAM_NO_SIMD drops all of them and
 * leaves only the scalar kernel.
 */
#if defined(SSAM_NO_SIMD)
   /* no architecture macro defined: scalar only */
#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#  define SSAM_ARCH_X86 1
#elif defined(__aarch64__) || defined(_M_ARM64)
#  define SSAM_ARCH_AARCH64 1
#endif

#if defined(SSAM_ARCH_X86)
#  ifdef _WIN32
#    include <intrin.h>
#  else
#    include <cpuid.h>
#  endif
#  include <immintrin.h> // AVX intrinsics
#  define SSAM_HAVE_AVX512 1
#  if defined(_MSC_VER)
#    define SSAM_AVX512_TARGET /* MSVC exposes the intrinsics unconditionally */
#  else
#    define SSAM_AVX512_TARGET __attribute__((target("avx512f")))
#  endif
#elif defined(SSAM_ARCH_AARCH64)
#  include <arm_neon.h>
#  define SSAM_HAVE_NEON 1
#  if !defined(SSAM_NO_SVE)
#    if defined(__ARM_FEATURE_SVE)
#      include <arm_sve.h>
#      define SSAM_HAVE_SVE 1
#      define SSAM_SVE_ALWAYS 1
#      define SSAM_SVE_TARGET
#    elif defined(__linux__) && \
          ((defined(__clang__) && __clang_major__ >= 16) || \
           (!defined(__clang__) && defined(__GNUC__) && __GNUC__ >= 12))
#      include <arm_sve.h>
#      include <sys/auxv.h>
#      define SSAM_HAVE_SVE 1
#      define SSAM_SVE_RUNTIME 1
#      define SSAM_SVE_TARGET __attribute__((target("+sve")))
#      ifndef HWCAP_SVE
#        define HWCAP_SVE (1 << 22)
#      endif
#    endif
#  endif
#endif

#include <Python.h>
#include "numpy/npy_math.h"
#include "numpy/arrayobject.h"

#define I2D(X, Y, YL) ((X) * (YL) + (Y))
#define I3D(X, Y, Z, YL, ZL) (((X) * (YL) * (ZL)) + ((Y) * (ZL)) + (Z))

struct pos2d {
    long x;
    long y;
};

struct pos3d {
    long x;
    long y;
    long z;
};

namespace std {
    template <>
    struct hash<pos3d> {
        std::size_t operator()(const pos3d &p) const {
            const std::size_t prime1 = 73856093;
            const std::size_t prime2 = 19349663;
            const std::size_t prime3 = 83492791;
            
            std::size_t hx = std::hash<long>()(p.x);
            std::size_t hy = std::hash<long>()(p.y);
            std::size_t hz = std::hash<long>()(p.z);
            
            return hx * prime1 ^ hy * prime2 ^ hz * prime3;
        }
    };
}

bool operator==(const pos3d& a, const pos3d& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

static double gauss_kernel(double x, double y, double z) {
    return exp(-0.5 * (x*x + y*y + z*z));
}

void kde(std::unordered_map<pos3d, double> &arr, double *xx, double *yy, double *zz, int *shape, int npts, double bandwidth, double prune_coeff, int ncores) {
    int maxdist;
    if (prune_coeff > 0) {
        maxdist = static_cast<int>(bandwidth * prune_coeff);
    } else {
        maxdist = -1;
    }

    #pragma omp parallel num_threads(ncores)
    {
        std::unordered_map<pos3d, double> local_map;

        #pragma omp for
        for (int i = 0; i < npts; i++) {
            int x0 = static_cast<int>(xx[i]);
            int y0 = static_cast<int>(yy[i]);
            int z0 = static_cast<int>(zz[i]);

            int xs = (maxdist > 0) ? std::max(0, x0 - maxdist) : 0;
            int xe = (maxdist > 0) ? std::min(shape[0], x0 + maxdist + 1) : shape[0];
            int ys = (maxdist > 0) ? std::max(0, y0 - maxdist) : 0;
            int ye = (maxdist > 0) ? std::min(shape[1], y0 + maxdist + 1) : shape[1];
            int zs = (maxdist > 0) ? std::max(0, z0 - maxdist) : 0;
            int ze = (maxdist > 0) ? std::min(shape[2], z0 + maxdist + 1) : shape[2];

            for (int x = xs; x < xe; x++) {
                for (int y = ys; y < ye; y++) {
                    for (int z = zs; z < ze; z++) {
                        pos3d key{x, y, z};
                        double value = gauss_kernel((x - xx[i]) / bandwidth, (y - yy[i]) / bandwidth, (z - zz[i]) / bandwidth);
                        local_map[key] += value;
                    }
                }
            }
        }

        #pragma omp critical
        for (const auto& pair : local_map) {
            arr[pair.first] += pair.second;
        }
    }
}

/*
 * Pearson's correlation coefficient.
 * Used by the scalar, NEON, SVE kernels for a single-pass "shifted data" formulation.
 */
static inline double corr_finalize(double n, double Sa, double Sb,
                                   double Saa, double Sbb, double Sab) {
    const double ma = Sa / n;
    const double mb = Sb / n;
    const double var_a = Saa / n - ma * ma;
    const double var_b = Sbb / n - mb * mb;
    const double cov = Sab / n - ma * mb;

    if (var_a <= 0.0 || var_b <= 0.0)
        return 0.0;

    return cov / (sqrt(var_a) * sqrt(var_b));
}

SSAM_MAYBE_UNUSED
static double corr_scalar(const double *a, const double *b, int ngene) {
    if (ngene < 2)
        return 0.0;

    const double sa = a[0], sb = b[0];
    double Sa = 0, Sb = 0, Saa = 0, Sbb = 0, Sab = 0;

    for (int i = 0; i < ngene; i++) {
        const double da = a[i] - sa;
        const double db = b[i] - sb;
        Sa += da;
        Sb += db;
        Saa += da * da;
        Sbb += db * db;
        Sab += da * db;
    }

    return corr_finalize((double)ngene, Sa, Sb, Saa, Sbb, Sab);
}

#if defined(SSAM_HAVE_NEON)
/*
 * NEON kernel. A NEON vector holds 2 doubles, so the loop is unrolled 8-wide
 * to keep enough independent accumulator chains in flight to hide the ~4 cycle
 * FMA latency of current Arm cores. The linear sums use two chains and the
 * three quadratic sums four each: 16 vector accumulators, which still leaves
 * room for the eight in-flight loads in the 32 entry register file.
 */
static double corr_neon(const double *a, const double *b, int ngene) {
    if (ngene < 2)
        return 0.0;

    const float64x2_t vsa = vdupq_n_f64(a[0]);
    const float64x2_t vsb = vdupq_n_f64(b[0]);
    const float64x2_t zero = vdupq_n_f64(0.0);

    float64x2_t Sa0 = zero, Sa1 = zero;
    float64x2_t Sb0 = zero, Sb1 = zero;
    float64x2_t Saa0 = zero, Saa1 = zero, Saa2 = zero, Saa3 = zero;
    float64x2_t Sbb0 = zero, Sbb1 = zero, Sbb2 = zero, Sbb3 = zero;
    float64x2_t Sab0 = zero, Sab1 = zero, Sab2 = zero, Sab3 = zero;

    int i = 0;
    for (; i <= ngene - 8; i += 8) {
        const float64x2_t da0 = vsubq_f64(vld1q_f64(a + i + 0), vsa);
        const float64x2_t da1 = vsubq_f64(vld1q_f64(a + i + 2), vsa);
        const float64x2_t da2 = vsubq_f64(vld1q_f64(a + i + 4), vsa);
        const float64x2_t da3 = vsubq_f64(vld1q_f64(a + i + 6), vsa);
        const float64x2_t db0 = vsubq_f64(vld1q_f64(b + i + 0), vsb);
        const float64x2_t db1 = vsubq_f64(vld1q_f64(b + i + 2), vsb);
        const float64x2_t db2 = vsubq_f64(vld1q_f64(b + i + 4), vsb);
        const float64x2_t db3 = vsubq_f64(vld1q_f64(b + i + 6), vsb);

        Sa0 = vaddq_f64(Sa0, vaddq_f64(da0, da1));
        Sa1 = vaddq_f64(Sa1, vaddq_f64(da2, da3));
        Sb0 = vaddq_f64(Sb0, vaddq_f64(db0, db1));
        Sb1 = vaddq_f64(Sb1, vaddq_f64(db2, db3));

        Saa0 = vfmaq_f64(Saa0, da0, da0);
        Saa1 = vfmaq_f64(Saa1, da1, da1);
        Saa2 = vfmaq_f64(Saa2, da2, da2);
        Saa3 = vfmaq_f64(Saa3, da3, da3);

        Sbb0 = vfmaq_f64(Sbb0, db0, db0);
        Sbb1 = vfmaq_f64(Sbb1, db1, db1);
        Sbb2 = vfmaq_f64(Sbb2, db2, db2);
        Sbb3 = vfmaq_f64(Sbb3, db3, db3);

        Sab0 = vfmaq_f64(Sab0, da0, db0);
        Sab1 = vfmaq_f64(Sab1, da1, db1);
        Sab2 = vfmaq_f64(Sab2, da2, db2);
        Sab3 = vfmaq_f64(Sab3, da3, db3);
    }

    for (; i <= ngene - 2; i += 2) {
        const float64x2_t da = vsubq_f64(vld1q_f64(a + i), vsa);
        const float64x2_t db = vsubq_f64(vld1q_f64(b + i), vsb);
        Sa0 = vaddq_f64(Sa0, da);
        Sb0 = vaddq_f64(Sb0, db);
        Saa0 = vfmaq_f64(Saa0, da, da);
        Sbb0 = vfmaq_f64(Sbb0, db, db);
        Sab0 = vfmaq_f64(Sab0, da, db);
    }

    double Sa = vaddvq_f64(vaddq_f64(Sa0, Sa1));
    double Sb = vaddvq_f64(vaddq_f64(Sb0, Sb1));
    double Saa = vaddvq_f64(vaddq_f64(vaddq_f64(Saa0, Saa1), vaddq_f64(Saa2, Saa3)));
    double Sbb = vaddvq_f64(vaddq_f64(vaddq_f64(Sbb0, Sbb1), vaddq_f64(Sbb2, Sbb3)));
    double Sab = vaddvq_f64(vaddq_f64(vaddq_f64(Sab0, Sab1), vaddq_f64(Sab2, Sab3)));

    /* At most one element left over. */
    for (; i < ngene; i++) {
        const double da = a[i] - a[0];
        const double db = b[i] - b[0];
        Sa += da;
        Sb += db;
        Saa += da * da;
        Sbb += db * db;
        Sab += da * db;
    }

    return corr_finalize((double)ngene, Sa, Sb, Saa, Sbb, Sab);
}
#endif /* SSAM_HAVE_NEON */

#if defined(SSAM_HAVE_SVE)
/*
 * SVE kernel. Vector-length agnostic: the same code runs on 128-bit
 * implementations (Neoverse V2, Cortex-X/A7xx) and on wider ones (Neoverse V1
 * at 256-bit, A64FX at 512-bit) without recompilation. Predication makes the
 * tail free, so short gene vectors do not fall back to scalar code.
 */
SSAM_SVE_TARGET
static double corr_sve(const double *a, const double *b, int ngene) {
    if (ngene < 2)
        return 0.0;

    const uint64_t vl = svcntd();
    const svfloat64_t vsa = svdup_f64(a[0]);
    const svfloat64_t vsb = svdup_f64(b[0]);
    const svbool_t ptrue = svptrue_b64();

    const svfloat64_t zero = svdup_f64(0.0);
    svfloat64_t Sa0 = zero, Sa1 = zero;
    svfloat64_t Sb0 = zero, Sb1 = zero;
    svfloat64_t Saa0 = zero, Saa1 = zero, Saa2 = zero, Saa3 = zero;
    svfloat64_t Sbb0 = zero, Sbb1 = zero, Sbb2 = zero, Sbb3 = zero;
    svfloat64_t Sab0 = zero, Sab1 = zero, Sab2 = zero, Sab3 = zero;

    uint64_t i = 0;
    const uint64_t n = (uint64_t)ngene;

    /* Four vectors per iteration. On a 128-bit implementation this matches
     * the NEON kernel instruction for instruction; on wider ones it covers the
     * same number of elements in proportionally fewer iterations. */
    if (n >= 4 * vl) {
        for (; i <= n - 4 * vl; i += 4 * vl) {
            const svfloat64_t da0 = svsub_f64_x(ptrue, svld1_f64(ptrue, a + i          ), vsa);
            const svfloat64_t da1 = svsub_f64_x(ptrue, svld1_f64(ptrue, a + i +      vl), vsa);
            const svfloat64_t da2 = svsub_f64_x(ptrue, svld1_f64(ptrue, a + i + 2 *  vl), vsa);
            const svfloat64_t da3 = svsub_f64_x(ptrue, svld1_f64(ptrue, a + i + 3 *  vl), vsa);
            const svfloat64_t db0 = svsub_f64_x(ptrue, svld1_f64(ptrue, b + i          ), vsb);
            const svfloat64_t db1 = svsub_f64_x(ptrue, svld1_f64(ptrue, b + i +      vl), vsb);
            const svfloat64_t db2 = svsub_f64_x(ptrue, svld1_f64(ptrue, b + i + 2 *  vl), vsb);
            const svfloat64_t db3 = svsub_f64_x(ptrue, svld1_f64(ptrue, b + i + 3 *  vl), vsb);

            Sa0 = svadd_f64_x(ptrue, Sa0, svadd_f64_x(ptrue, da0, da1));
            Sa1 = svadd_f64_x(ptrue, Sa1, svadd_f64_x(ptrue, da2, da3));
            Sb0 = svadd_f64_x(ptrue, Sb0, svadd_f64_x(ptrue, db0, db1));
            Sb1 = svadd_f64_x(ptrue, Sb1, svadd_f64_x(ptrue, db2, db3));

            Saa0 = svmla_f64_x(ptrue, Saa0, da0, da0);
            Saa1 = svmla_f64_x(ptrue, Saa1, da1, da1);
            Saa2 = svmla_f64_x(ptrue, Saa2, da2, da2);
            Saa3 = svmla_f64_x(ptrue, Saa3, da3, da3);

            Sbb0 = svmla_f64_x(ptrue, Sbb0, db0, db0);
            Sbb1 = svmla_f64_x(ptrue, Sbb1, db1, db1);
            Sbb2 = svmla_f64_x(ptrue, Sbb2, db2, db2);
            Sbb3 = svmla_f64_x(ptrue, Sbb3, db3, db3);

            Sab0 = svmla_f64_x(ptrue, Sab0, da0, db0);
            Sab1 = svmla_f64_x(ptrue, Sab1, da1, db1);
            Sab2 = svmla_f64_x(ptrue, Sab2, da2, db2);
            Sab3 = svmla_f64_x(ptrue, Sab3, da3, db3);
        }
    }

    /* Predicated remainder: no scalar epilogue. */
    for (; i < n; i += vl) {
        const svbool_t pg = svwhilelt_b64(i, n);
        const svfloat64_t da = svsub_f64_z(pg, svld1_f64(pg, a + i), vsa);
        const svfloat64_t db = svsub_f64_z(pg, svld1_f64(pg, b + i), vsb);

        Sa0 = svadd_f64_m(pg, Sa0, da);
        Sb0 = svadd_f64_m(pg, Sb0, db);
        Saa0 = svmla_f64_m(pg, Saa0, da, da);
        Sbb0 = svmla_f64_m(pg, Sbb0, db, db);
        Sab0 = svmla_f64_m(pg, Sab0, da, db);
    }

    Saa0 = svadd_f64_x(ptrue, svadd_f64_x(ptrue, Saa0, Saa1), svadd_f64_x(ptrue, Saa2, Saa3));
    Sbb0 = svadd_f64_x(ptrue, svadd_f64_x(ptrue, Sbb0, Sbb1), svadd_f64_x(ptrue, Sbb2, Sbb3));
    Sab0 = svadd_f64_x(ptrue, svadd_f64_x(ptrue, Sab0, Sab1), svadd_f64_x(ptrue, Sab2, Sab3));

    return corr_finalize((double)ngene,
                         svaddv_f64(ptrue, svadd_f64_x(ptrue, Sa0, Sa1)),
                         svaddv_f64(ptrue, svadd_f64_x(ptrue, Sb0, Sb1)),
                         svaddv_f64(ptrue, Saa0),
                         svaddv_f64(ptrue, Sbb0),
                         svaddv_f64(ptrue, Sab0));
}
#endif /* SSAM_HAVE_SVE */

#if defined(SSAM_SVE_RUNTIME)
/* 0 = NEON, 1 = SVE, -1 = not probed yet. */
static int ssam_use_sve = -1;

static int ssam_probe_sve(void) {
    return (getauxval(AT_HWCAP) & HWCAP_SVE) != 0;
}
#endif

#if defined(SSAM_HAVE_AVX512)
SSAM_AVX512_TARGET
static double corr_avx512(double *a, double *b, int ngene) {
    __m512d sum_a = _mm512_setzero_pd();
    __m512d sum_b = _mm512_setzero_pd();
    __m512d sum_aa = _mm512_setzero_pd();
    __m512d sum_bb = _mm512_setzero_pd();
    __m512d sum_ab_diff = _mm512_setzero_pd();

    int i;

    // Handle bulk of data with AVX-512
    for (i = 0; i <= ngene - 8; i += 8) {
        __m512d va = _mm512_loadu_pd(&a[i]);
        __m512d vb = _mm512_loadu_pd(&b[i]);

        sum_a = _mm512_add_pd(sum_a, va);
        sum_b = _mm512_add_pd(sum_b, vb);
        
        sum_aa = _mm512_add_pd(sum_aa, _mm512_mul_pd(va, va));
        sum_bb = _mm512_add_pd(sum_bb, _mm512_mul_pd(vb, vb));
    }

    double arr_a[8], arr_b[8], arr_aa[8], arr_bb[8];
    _mm512_storeu_pd(arr_a, sum_a);
    _mm512_storeu_pd(arr_b, sum_b);
    _mm512_storeu_pd(arr_aa, sum_aa);
    _mm512_storeu_pd(arr_bb, sum_bb);

    double a_mean = 0, b_mean = 0, aa_mean = 0, bb_mean = 0;
    for(int j = 0; j < 8; j++) {
        a_mean += arr_a[j];
        b_mean += arr_b[j];
        aa_mean += arr_aa[j];
        bb_mean += arr_bb[j];
    }
    a_mean /= ngene;
    b_mean /= ngene;
    aa_mean /= ngene;
    bb_mean /= ngene;

    // Handle remainder
    for (; i < ngene; i++) {
        a_mean += a[i] / ngene;
        b_mean += b[i] / ngene;
        aa_mean += a[i] * a[i] / ngene;
        bb_mean += b[i] * b[i] / ngene;
    }

    double a_std = sqrt(aa_mean - a_mean * a_mean);
    double b_std = sqrt(bb_mean - b_mean * b_mean);
    double rtn = 0;

    if (a_std != 0 && b_std != 0) {
        for (i = 0; i <= ngene - 8; i += 8) {
            __m512d va = _mm512_loadu_pd(&a[i]);
            __m512d vb = _mm512_loadu_pd(&b[i]);

            va = _mm512_sub_pd(va, _mm512_set1_pd(a_mean));
            vb = _mm512_sub_pd(vb, _mm512_set1_pd(b_mean));

            sum_ab_diff = _mm512_add_pd(sum_ab_diff, _mm512_mul_pd(va, vb));
        }

        double arr_ab_diff[8];
        _mm512_storeu_pd(arr_ab_diff, sum_ab_diff);
        for(int j = 0; j < 8; j++) {
            rtn += arr_ab_diff[j];
        }

        for (; i < ngene; i++) {
            rtn += (a[i] - a_mean) * (b[i] - b_mean);
        }

        rtn /= (a_std * b_std * ngene);
    }

    return rtn;
}

static int ssam_probe_avx512f(void) {
    int info[4];
    
    #ifdef _WIN32
    __cpuidex(info, 7, 0);
    #else
    __cpuid_count(7, 0, info[0], info[1], info[2], info[3]);
    #endif

    return (info[1] & (1 << 16)) != 0;  // Check if the 16th bit of EBX is set
}
#endif /* SSAM_HAVE_AVX512 */

/*
 * Dispatch to the widest kernel this machine actually supports.
 */
static double __corr__(double *a, double *b, int ngene) {
#if defined(SSAM_HAVE_AVX512)
    static int use_avx512 = -1;
    if (use_avx512 < 0)
        use_avx512 = ssam_probe_avx512f();
    if (use_avx512)
        return corr_avx512(a, b, ngene);
    return corr_scalar(a, b, ngene);
#elif defined(SSAM_SVE_ALWAYS)
    return corr_sve(a, b, ngene);
#elif defined(SSAM_SVE_RUNTIME)
    if (ssam_use_sve < 0)
        ssam_use_sve = ssam_probe_sve();
    if (ssam_use_sve)
        return corr_sve(a, b, ngene);
    return corr_neon(a, b, ngene);
#elif defined(SSAM_HAVE_NEON)
    return corr_neon(a, b, ngene);
#else
    return corr_scalar(a, b, ngene);
#endif
}

/* Name of the kernel __corr__ will actually run on this machine. */
static const char *ssam_simd_backend(void) {
#if defined(SSAM_HAVE_AVX512)
    return ssam_probe_avx512f() ? "avx512f" : "scalar";
#elif defined(SSAM_SVE_ALWAYS)
    return "sve";
#elif defined(SSAM_SVE_RUNTIME)
    return ssam_probe_sve() ? "sve" : "neon";
#elif defined(SSAM_HAVE_NEON)
    return "neon";
#else
    return "scalar";
#endif
}

static PyObject *simd_backend(PyObject *self, PyObject *args) {
    return PyUnicode_FromString(ssam_simd_backend());
}

static PyObject *calc_kde(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    PyObject *arg4 = NULL;
    PyArrayObject *arr1 = NULL;
    PyArrayObject *arr2 = NULL;
    PyArrayObject *arr3 = NULL;
    PyArrayObject *arr4 = NULL;
    int ncores = omp_get_max_threads();
    PyObject *rtn, *poslist, *xlist, *ylist, *zlist, *vlist;
    double *x, *y, *z;
    int *shape;
    double h, prune_coeff;
    int kernel = 0;
    unsigned int npts;
    int cnt;
    std::unordered_map<pos3d, double> oarr_map;
    std::unordered_map<pos3d, double>::iterator it;

    static const char *kwlist[] = { "h", "x", "y", "z", "shape", "prune_coeff", "kernel", "ncores", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dOOOOd|ii", const_cast<char **>(kwlist), &h, &arg1, &arg2, &arg3, &arg4, &prune_coeff, &kernel, &ncores)) return NULL;
    if ((arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) return NULL;
    if ((arr2 = (PyArrayObject*)PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) goto fail;
    if ((arr3 = (PyArrayObject*)PyArray_FROM_OTF(arg3, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) goto fail;
    if ((arr4 = (PyArrayObject*)PyArray_FROM_OTF(arg4, NPY_INT, NPY_ARRAY_IN_ARRAY)) == NULL) goto fail;
    
    if (PyArray_NDIM(arr1) != 1 || PyArray_NDIM(arr2) != 1 || PyArray_NDIM(arr3) != 1 || PyArray_NDIM(arr4) != 1)
    {
        goto fail;
    }

    npts = PyArray_DIMS(arr1)[0];

    x = (double *)PyArray_DATA(arr1);
    y = (double *)PyArray_DATA(arr2);
    z = (double *)PyArray_DATA(arr3);
    shape = (int *)PyArray_DATA(arr4);

    kde(oarr_map, x, y, z, shape, npts, h, prune_coeff, ncores);
    rtn = (PyObject *)PyTuple_New(2);
    poslist = (PyObject *)PyList_New(3);
    xlist = (PyObject *)PyList_New(oarr_map.size());
    ylist = (PyObject *)PyList_New(oarr_map.size());
    zlist = (PyObject *)PyList_New(oarr_map.size());
    vlist = (PyObject *)PyList_New(oarr_map.size());
    
    cnt = 0;
    for (const auto& pair : oarr_map) {
        PyList_SetItem(xlist, cnt, PyLong_FromLong((long)(pair.first.x)));
        PyList_SetItem(ylist, cnt, PyLong_FromLong((long)(pair.first.y)));
        PyList_SetItem(zlist, cnt, PyLong_FromLong((long)(pair.first.z)));
        PyList_SetItem(vlist, cnt, PyFloat_FromDouble(pair.second));
        cnt++;
    }
    PyList_SetItem(poslist, 0, xlist);
    PyList_SetItem(poslist, 1, ylist);
    PyList_SetItem(poslist, 2, zlist);
    PyTuple_SetItem(rtn, 0, poslist);
    PyTuple_SetItem(rtn, 1, vlist);
    
    Py_DECREF(arr1);
    Py_DECREF(arr2);
    Py_DECREF(arr3);
    Py_DECREF(arr4);
    
    return (PyObject *) rtn;
    
fail:
    Py_XDECREF(arr1);
    Py_XDECREF(arr2);
    Py_XDECREF(arr3);
    Py_XDECREF(arr4);
    return NULL;
}

static PyObject *flood_fill(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject* filled_poslist = NULL;
    PyArrayObject *arr1 = NULL;
    PyArrayObject *arr2 = NULL;
    long nvec, nd, ngene = 0;
    long *pos, x, y, z, cnt;
    double r = 0.6, *vf;
    npy_intp *dimsp;
    int min_pixels = 10, max_pixels=2000;
    int i;
    bool *mask;

    static const char *kwlist[] = { "pos", "vf", "r", "min_pixels", "max_pixels", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|dii", const_cast<char **>(kwlist), &arg1, &arg2, &r, &min_pixels, &max_pixels)) return NULL;
    if ((arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_LONG, NPY_ARRAY_IN_ARRAY)) == NULL) return NULL;
    if ((arr2 = (PyArrayObject*)PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) goto fail;
    if (PyArray_NDIM(arr1) != 1) goto fail;
    nd = PyArray_NDIM(arr2);
    dimsp = PyArray_DIMS(arr2);
    nvec = 1;
    for (i=0; i<nd-1; i++)
        nvec *= dimsp[i];
    ngene = dimsp[nd-1];
    filled_poslist = PyList_New(0);
    pos = (long *)PyArray_DATA(arr1);
    vf = (double *)PyArray_DATA(arr2);
    mask = (bool *) calloc(nvec, sizeof(bool));
    cnt = 0;
    if (nd == 3) {
        // 2D
        std::queue<pos2d> queue2d;
        queue2d.push(pos2d());
        queue2d.back().x = pos[0];
        queue2d.back().y = pos[1];
        while (queue2d.size() > 0) {
            x = queue2d.front().x;
            y = queue2d.front().y;
            PyObject *t = PyTuple_New(2);
            PyTuple_SetItem(t, 0, PyLong_FromLong(x));
            PyTuple_SetItem(t, 1, PyLong_FromLong(y));
            cnt += 1;
            if (cnt > max_pixels) 
                break;
            PyList_Append(filled_poslist, t);
            queue2d.pop();
            if (x < dimsp[0] - 1 && mask[I2D(x + 1, y, dimsp[1])] == false &&
                    __corr__(vf + (I2D(pos[0], pos[1], dimsp[1]) * ngene),
                             vf + (I2D(x + 1, y, dimsp[1]) * ngene), ngene) > r) {
                mask[I2D(x + 1, y, dimsp[1])] = true;
                queue2d.push(pos2d());
                queue2d.back().x = x + 1;
                queue2d.back().y = y;
            }
            if (x > 1 && mask[I2D(x - 1, y, dimsp[1])] == false &&
                    __corr__(vf + (I2D(pos[0], pos[1], dimsp[1]) * ngene),
                             vf + (I2D(x - 1, y, dimsp[1]) * ngene), ngene) > r) {
                mask[I2D(x - 1, y, dimsp[1])] = true;
                queue2d.push(pos2d());
                queue2d.back().x = x - 1;
                queue2d.back().y = y;
            }
            if (y < dimsp[1] - 1 && mask[I2D(x, y + 1, dimsp[1])] == false &&
                    __corr__(vf + (I2D(pos[0], pos[1], dimsp[1]) * ngene),
                             vf + (I2D(x, y + 1, dimsp[1]) * ngene), ngene) > r) {
                mask[I2D(x, y + 1, dimsp[1])] = true;
                queue2d.push(pos2d());
                queue2d.back().x = x;
                queue2d.back().y = y + 1;
            }
            if (y > 1 && mask[I2D(x, y - 1, dimsp[1])] == false &&
                    __corr__(vf + (I2D(pos[0], pos[1], dimsp[1]) * ngene),
                             vf + (I2D(x, y - 1, dimsp[1]) * ngene), ngene) > r) {
                mask[I2D(x, y - 1, dimsp[1])] = true;
                queue2d.push(pos2d());
                queue2d.back().x = x;
                queue2d.back().y = y - 1;
            }
        }
    } else if (nd == 4) {
        // 3D
        std::queue<pos3d> queue3d;
        queue3d.push(pos3d());
        queue3d.back().x = pos[0];
        queue3d.back().y = pos[1];
        queue3d.back().z = pos[2];
        while (queue3d.size() > 0) {
            x = queue3d.front().x;
            y = queue3d.front().y;
            z = queue3d.front().z;
            PyObject *t = PyTuple_New(3);
            PyTuple_SetItem(t, 0, PyLong_FromLong(x));
            PyTuple_SetItem(t, 1, PyLong_FromLong(y));
            PyTuple_SetItem(t, 2, PyLong_FromLong(z));
            PyList_Append(filled_poslist, t);
            cnt += 1;
            if (cnt > max_pixels) 
                break;
            queue3d.pop();
            if (x < dimsp[0] - 1 && mask[I3D(x + 1, y, z, dimsp[1], dimsp[2])] == false &&
                    __corr__(vf + I3D(pos[0], pos[1], pos[2], dimsp[1], dimsp[2]) * ngene,
                             vf + I3D(x + 1, y, z, dimsp[1], dimsp[2]) * ngene, ngene) > r) {
                mask[I3D(x + 1, y, z, dimsp[1], dimsp[2])] = true;
                queue3d.push(pos3d());
                queue3d.back().x = x + 1;
                queue3d.back().y = y;
                queue3d.back().z = z;
            }
            if (x > 1 && mask[I3D(x - 1, y, z, dimsp[1], dimsp[2])] == false &&
                    __corr__(vf + I3D(pos[0], pos[1], pos[2], dimsp[1], dimsp[2]) * ngene,
                             vf + I3D(x - 1, y, z, dimsp[1], dimsp[2]) * ngene, ngene) > r) {
                mask[I3D(x - 1, y, z, dimsp[1], dimsp[2])] = true;
                queue3d.push(pos3d());
                queue3d.back().x = x - 1;
                queue3d.back().y = y;
                queue3d.back().z = z;
            }
            if (y < dimsp[1] - 1 && mask[I3D(x, y + 1, z, dimsp[1], dimsp[2])] == false &&
                    __corr__(vf + I3D(pos[0], pos[1], pos[2], dimsp[1], dimsp[2]) * ngene,
                             vf + I3D(x, y + 1, z, dimsp[1], dimsp[2]) * ngene, ngene) > r) {
                mask[I3D(x, y + 1, z, dimsp[1], dimsp[2])] = true;
                queue3d.push(pos3d());
                queue3d.back().x = x;
                queue3d.back().y = y + 1;
                queue3d.back().z = z;
            }
            if (y > 1 && mask[I3D(x, y - 1, z, dimsp[1], dimsp[2])] == false &&
                    __corr__(vf + I3D(pos[0], pos[1], pos[2], dimsp[1], dimsp[2]) * ngene,
                             vf + I3D(x, y - 1, z, dimsp[1], dimsp[2]) * ngene, ngene) > r) {
                mask[I3D(x, y - 1, z, dimsp[1], dimsp[2])] = true;
                queue3d.push(pos3d());
                queue3d.back().x = x;
                queue3d.back().y = y - 1;
                queue3d.back().z = z;
            }
            if (z < dimsp[2] - 1 && mask[I3D(x, y, z + 1, dimsp[1], dimsp[2])] == false &&
                    __corr__(vf + I3D(pos[0], pos[1], pos[2], dimsp[1], dimsp[2]) * ngene,
                             vf + I3D(x, y, z + 1, dimsp[1], dimsp[2]) * ngene, ngene) > r) {
                mask[I3D(x, y, z, dimsp[1], dimsp[2])] = true;
                queue3d.push(pos3d());
                queue3d.back().x = x;
                queue3d.back().y = y;
                queue3d.back().z = z + 1;
            }
            if (z > 1 && mask[I3D(x, y, z - 1, dimsp[1], dimsp[2])] == false &&
                    __corr__(vf + I3D(pos[0], pos[1], pos[2], dimsp[1], dimsp[2]) * ngene,
                             vf + I3D(x, y, z - 1, dimsp[1], dimsp[2]) * ngene, ngene) > r) {
                mask[I3D(x, y, z - 1, dimsp[1], dimsp[2])] = true;
                queue3d.push(pos3d());
                queue3d.back().x = x;
                queue3d.back().y = y;
                queue3d.back().z = z - 1;
            }
        }
    }
    free((void*)mask);
    Py_DECREF(arr1);
    Py_DECREF(arr2);
    if (cnt > max_pixels || cnt < min_pixels) 
        PyList_SetSlice(filled_poslist, 0, PyList_Size(filled_poslist), NULL);
    return (PyObject *) filled_poslist;
 
fail:
    Py_XDECREF(arr1);
    Py_XDECREF(arr2);
    return NULL;
}

static PyObject *calc_corrmap(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *arg1 = NULL;
    PyArrayObject *arr1 = NULL;
    PyArrayObject *oarr = NULL;
    long i, x, y, z, dx, dy, dz;
    long nvec, nd, ngene = 0;
    double *vecs, *corrmap;
    npy_intp *dimsp;
    int ncores = omp_get_max_threads();
    int csize = 1;
    double *tmpvec;

    static const char *kwlist[] = { "vf", "ncores", "size", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ii", const_cast<char **>(kwlist), &arg1, &ncores, &csize)) return NULL;
    if ((arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) return NULL;   
    nd = PyArray_NDIM(arr1);
    if (nd != 3 && nd != 4) goto fail; // only 2D or 3D array is expected
    dimsp = PyArray_DIMS(arr1);
    oarr = (PyArrayObject*)PyArray_ZEROS(nd - 1, dimsp, NPY_DOUBLE, NPY_CORDER);
    ngene = dimsp[nd-1];
    corrmap = (double *)PyArray_DATA(oarr);
    vecs = (double *)PyArray_DATA(arr1);
    nvec = 1;
    for (i=0; i<nd-1; i++)
        nvec *= dimsp[i];
    
    // initialize corrmap with NANs
    #pragma omp parallel for num_threads(ncores)
    for (i=0; i<nvec; i++)
        corrmap[i] = NPY_NAN;

    if (nd == 3) {
        // 2D
        #pragma omp parallel num_threads(ncores) private(tmpvec)
        {
            tmpvec = (double *)calloc(ngene, sizeof(double)); // zero initialized
            #pragma omp for collapse(2)
            for (x=csize; x<dimsp[0]-csize; x++) {
                for (y=csize; y<dimsp[1]-csize; y++) {
                    for (i=0; i<ngene; i++)
                        tmpvec[i] = 0;
                    for (dx=-csize; dx<csize+1; dx++) {
                        for (dy=-csize; dy<csize+1; dy++) {
                            if (dx == 0 && dy == 0) continue;
                            for (i=0; i<ngene; i++)
                                tmpvec[i] += (vecs + I2D(x+dx, y+dy, dimsp[1])*ngene)[i];
                        }
                    }
                    // tmpvec[i] /= (csize * 2 + 1) * (csize * 2 + 1) - 1;
                    corrmap[I2D(x, y, dimsp[1])] = __corr__(vecs + I2D(x, y, dimsp[1])*ngene, tmpvec, ngene);
                }
            }
            free((void*)tmpvec);
        }
    } else {
        // 3D (nd == 4)
        #pragma omp parallel num_threads(ncores) private(tmpvec)
        {
            tmpvec = (double *)calloc(ngene, sizeof(double));
            #pragma omp for collapse(3)
            for (x=csize; x<dimsp[0]-csize; x++) {
                for (y=csize; y<dimsp[1]-csize; y++) {
                    for (z=csize; z<dimsp[2]-csize; z++) {
                        for (i=0; i<ngene; i++)
                            tmpvec[i] = 0;
                        for (dx=-csize; dx<csize+1; dx++) {
                            for (dy=-csize; dy<csize+1; dy++) {
                                for (dz=-csize; dz<csize+1; dz++) {
                                    if (dx == 0 && dy == 0 && dz == 0) continue;
                                    for (i=0; i<ngene; i++)
                                        tmpvec[i] += (vecs + I3D(x+dx, y+dy, z+dz, dimsp[1], dimsp[2])*ngene)[i];
                                }
                            }
                        }
                        //for (i=0; i<ngene; i++)
                        //    tmpvec[i] /= (csize * 2 + 1) * (csize * 2 + 1) * (csize * 2 + 1) - 1;
                        corrmap[I3D(x, y, z, dimsp[1], dimsp[2])] =
                            __corr__(vecs + I3D(x, y, z, dimsp[1], dimsp[2])*ngene, tmpvec, ngene);
                    }
                }
            }
            free((void*)tmpvec);
        }
    }
    Py_DECREF(arr1);

    return (PyObject *) oarr;
 fail:
    Py_XDECREF(arr1);
    return NULL;
}

static PyObject *calc_corrmap_2(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *arg1 = NULL;
    PyArrayObject *arr1 = NULL;
    PyArrayObject *oarr = NULL;
    long i, k, x, y, z, dx, dy, dz;
    long nvec, nd, ngene = 0;
    double *vecs, *corrmap;
    npy_intp *dimsp;
    npy_intp dimsp2[4];
    int ncores = omp_get_max_threads();
    int csize = 1;

    static const char *kwlist[] = { "vf", "ncores", "size", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ii", const_cast<char **>(kwlist), &arg1, &ncores, &csize)) return NULL;
    if ((arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) return NULL;   
    nd = PyArray_NDIM(arr1);
    if (nd != 3 && nd != 4) goto fail; // only 2D or 3D array is expected
    dimsp = PyArray_DIMS(arr1);
    for (i=0; i<nd-1; i++)
        dimsp2[i] = dimsp[i];
    dimsp2[nd-1] = (csize * 2 + 1) * (csize * 2 + 1) - 1;
    oarr = (PyArrayObject*)PyArray_ZEROS(nd, dimsp2, NPY_DOUBLE, NPY_CORDER);
    ngene = dimsp[nd-1];
    corrmap = (double *)PyArray_DATA(oarr);
    vecs = (double *)PyArray_DATA(arr1);
    nvec = 1;
    for (i=0; i<nd; i++)
        nvec *= dimsp2[i];

    // initialize corrmap with NANs
    #pragma omp parallel for num_threads(ncores)
    for (i=0; i<nvec; i++)
        corrmap[i] = NPY_NAN;

    if (nd == 3) {
        // 2D
        #pragma omp parallel for collapse(2) private(k, dx, dy)
        for (x=csize; x<dimsp[0]-csize; x++) {
            for (y=csize; y<dimsp[1]-csize; y++) {
                k = 0;
                for (dx=-csize; dx<csize+1; dx++) {
                    for (dy=-csize; dy<csize+1; dy++) {
                        if (dx == 0 && dy == 0) continue;
                        corrmap[I2D(x, y, dimsp2[1])*dimsp2[2] + (k++)] = 
                            __corr__(vecs + I2D(x, y, dimsp[1])*ngene,
                                     vecs + I2D(x+dx, y+dy, dimsp[1])*ngene, ngene);
                    }
                }
            }
        }
    } else {
        // 3D
        #pragma omp parallel for collapse(3) private(k, dx, dy, dz)
        for (x=csize; x<dimsp[0]-csize; x++) {
            for (y=csize; y<dimsp[1]-csize; y++) {
                for (z=csize; z<dimsp[2]-csize; z++) {
                    k = 0;
                    for (dx=-csize; dx<csize+1; dx++) {
                        for (dy=-csize; dy<csize+1; dy++) {
                            for (dz=-csize; dz<csize+1; dz++) {
                                if (dx == 0 && dy == 0 && dz == 0) continue;
                                corrmap[I3D(x, y, z, dimsp2[1], dimsp2[2])*dimsp2[3] + (k++)] =
                                    __corr__(vecs + I3D(x, y, z, dimsp[1], dimsp[2])*ngene, vecs + I3D(x+dx, y+dy, z+dz, dimsp[1], dimsp[2])*ngene, ngene);
                            }
                        }
                    }
                }
            }
        }
    }
    Py_DECREF(arr1);

    return (PyObject *) oarr;
 fail:
    Py_XDECREF(arr1);
    return NULL;
}

static PyObject *calc_ctmap(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyArrayObject *arr1 = NULL;
    PyArrayObject *arr2 = NULL;
    PyArrayObject *oarr = NULL;
    long nvec, nd, ngene = 0;
    double *cent, *vecs, *scores;
    npy_intp *dimsp;
    int ncores = omp_get_max_threads();
    int i;

    static const char *kwlist[] = { "vec", "vf", "ncores", NULL };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|i", const_cast<char **>(kwlist), &arg1, &arg2, &ncores)) return NULL;
    if ((arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) return NULL;
    if ((arr2 = (PyArrayObject*)PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) goto fail;
    if (PyArray_NDIM(arr1) != 1) goto fail;
    nd = PyArray_NDIM(arr2);
    if((ngene = *PyArray_DIMS(arr1)) != PyArray_DIMS(arr2)[nd-1]) goto fail;

    dimsp = PyArray_DIMS(arr2);
    oarr = (PyArrayObject*)PyArray_ZEROS(nd - 1, dimsp, NPY_DOUBLE, NPY_CORDER);

    nvec = 1;
    for (i=0; i<nd-1; i++)
        nvec *= dimsp[i];

    scores = (double *)PyArray_DATA(oarr);
    cent = (double *)PyArray_DATA(arr1);
    vecs = (double *)PyArray_DATA(arr2);

    #pragma omp parallel for num_threads(ncores)
    for (i=0; i<nvec; i++) {
        scores[i] = __corr__(cent, vecs + (i*ngene), ngene);
    }

    Py_DECREF(arr1);
    Py_DECREF(arr2);

    return (PyObject *) oarr;
 fail:
    Py_XDECREF(arr1);
    Py_XDECREF(arr2);
    return NULL;
}

static PyObject *corr(PyObject *self, PyObject *args) {
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyArrayObject *arr1 = NULL;
    PyArrayObject *arr2 = NULL;
    long ngene = 0;
    double *a;
    double *b;
    double rtn = 0;

    if (!PyArg_ParseTuple(args, "OO", &arg1, &arg2)) return NULL;
    if ((arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) return NULL;
    if ((arr2 = (PyArrayObject*)PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) goto fail;
    if (PyArray_NDIM(arr1) != 1) goto fail;
    if (PyArray_NDIM(arr2) != 1) goto fail;
    if((ngene = *PyArray_DIMS(arr1)) != *PyArray_DIMS(arr2)) goto fail;

    a = (double *)PyArray_DATA(arr1);
    b = (double *)PyArray_DATA(arr2);

    rtn = __corr__(a, b, ngene);

    Py_DECREF(arr1);
    Py_DECREF(arr2);

    return PyFloat_FromDouble(rtn);
 fail:
    Py_XDECREF(arr1);
    Py_XDECREF(arr2);
    return NULL;
}

static struct PyMethodDef module_methods[] = {
    {"corr", (PyCFunction)corr, METH_VARARGS, "Calculates Pearson's correlation coefficient."},
    {"calc_ctmap", (PyCFunction)calc_ctmap, METH_VARARGS | METH_KEYWORDS, "Creates a cell type map."},
    {"calc_corrmap", (PyCFunction)calc_corrmap, METH_VARARGS | METH_KEYWORDS, "Creates a correlation map."},
    {"calc_corrmap_2", (PyCFunction)calc_corrmap_2, METH_VARARGS | METH_KEYWORDS, "Creates a correlation map."},
    {"calc_kde", (PyCFunction)calc_kde, METH_VARARGS | METH_KEYWORDS, "Run kernel density estimation."},
    {"flood_fill", (PyCFunction)flood_fill, METH_VARARGS | METH_KEYWORDS, "Performs 3d flood fill based on correlation."},
    {"simd_backend", (PyCFunction)simd_backend, METH_NOARGS, "Name of the SIMD kernel in use ('avx512f', 'sve', 'neon' or 'scalar')."},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "analysis_utils",
        NULL,
        -1,
        module_methods
};
#endif

PyMODINIT_FUNC
PyInit_utils(void)
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    Py_InitModule("utils", module_methods);
#endif
    import_array();
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}