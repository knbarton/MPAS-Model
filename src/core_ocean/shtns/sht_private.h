/*
 * Copyright (c) 2010-2020 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@univ-grenoble-alpes.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

/********************************************************************
 * SHTns : Spherical Harmonic Transform for numerical simulations.  *
 *    written by Nathanael Schaeffer / CNRS                         *
 ********************************************************************/

/// \internal \file sht_private.h private data and options.

#include <stdlib.h>
#include <string.h>
#ifndef __cplusplus
  // with C, include complex.h before fftw3.h to make "fftw_complex" match "complex double".
  #include <complex.h>
#endif
#include <math.h>
// FFTW la derivee d/dx = ik	(pas de moins !)
#include "fftw3/fftw3.h"

// config file generated by ./configure
#include "sht_config.h"

#define SHTNS_PRIVATE
#include "shtns.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

#ifdef HAVE_LIBCUFFT
#include <cufft.h>
#include "shtns_cuda.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/// private gpu functions:
int cushtns_init_gpu(shtns_cfg);
void cushtns_release_gpu(shtns_cfg);
int cushtns_use_gpu(int);
#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif

/* BEGIN COMPILE-TIME SETTINGS */

/// defines the maximum amount of memory in megabytes that SHTns should use.
#define SHTNS_MAX_MEMORY 2048

/// Minimum performance improve for DCT in \ref sht_auto mode. If not atained, we may switch back to gauss.
#define MIN_PERF_IMPROVE_DCT 1.05

/// Try to enforce at least this accuracy for DCT in sht_auto mode.
#define MIN_ACCURACY_DCT 1.e-8

/// The default \ref opt_polar threshold (0 disabled, 1.e-6 is aggressive, 1.e-10 is safe, 1.e-14 is VERY safe)
#define SHT_DEFAULT_POLAR_OPT 1.e-10

/// The default \ref norm used by shtns_init
#define SHT_DEFAULT_NORM ( sht_orthonormal )
//#define SHT_DEFAULT_NORM ( sht_schmidt | SHT_NO_CS_PHASE )

/// The maximum order of non-linear terms to be resolved by SH transform by default.
/// 1 : no non-linear terms. 2 : quadratic non-linear terms (default), 3 : triadic, ...
/// must be larger or equal to 1.
#define SHT_DEFAULT_NL_ORDER 1

/// minimum NLAT to consider the use of DCT acceleration.
#define SHT_MIN_NLAT_DCT 64

/// time-limit for timing individual transforms (in seconds)
#define SHT_TIME_LIMIT 0.2

/* END COMPILE-TIME SETTINGS */

// sht variants (std, ltr)
enum sht_variants { SHT_STD, SHT_M, SHT_NVAR };
// sht types (scal synth, scal analys, vect synth, ...)
enum sht_types { SHT_TYP_SSY, SHT_TYP_SAN, SHT_TYP_VSY, SHT_TYP_VAN,
	SHT_TYP_GSP, SHT_TYP_GTO, SHT_TYP_3SY, SHT_TYP_3AN, SHT_NTYP };

// sht grids
enum sht_grids { GRID_NONE, GRID_GAUSS, GRID_REGULAR, GRID_POLES };

// pointer to various function types
typedef void (*pf2l)(shtns_cfg, void*, void*, long int);
typedef void (*pf3l)(shtns_cfg, void*, void*, void*, long int);
typedef void (*pf4l)(shtns_cfg, void*, void*, void*, void*, long int);
typedef void (*pf6l)(shtns_cfg, void*, void*, void*, void*, void*, void*, long int);
typedef void (*pf2ml)(shtns_cfg, int, void*, void*, long int);
typedef void (*pf3ml)(shtns_cfg, int, void*, void*, void*, long int);
typedef void (*pf4ml)(shtns_cfg, int, void*, void*, void*, void*, long int);
typedef void (*pf6ml)(shtns_cfg, int, void*, void*, void*, void*, void*, void*, long int);

/// structure containing useful information about the SHT.
struct shtns_info {		// MUST start with "int nlm;"
/* PUBLIC PART (if modified, shtns.h should be modified acordingly) */
	unsigned int nlm;			///< total number of (l,m) spherical harmonics components.
	unsigned short lmax;		///< maximum degree (lmax) of spherical harmonics.
	unsigned short mmax;		///< maximum order (mmax*mres) of spherical harmonics.
	unsigned short mres;		///< the periodicity along the phi axis.
	unsigned short nlat_2;		///< ...and half of it (using (shtns.nlat+1)/2 allows odd shtns.nlat.)
	unsigned int nlat;			///< number of spatial points in Theta direction (latitude) ...
	unsigned int nphi;			///< number of spatial points in Phi direction (longitude)
	unsigned int nspat;			///< number of real numbers that must be allocated in a spatial field.
	unsigned short *li;			///< degree l for given mode index (size nlm) : li[lm]
	unsigned short *mi;			///< order m for given mode index (size nlm) : mi[lm]
	double *ct, *st;			///< cos(theta) and sin(theta) arrays (size nlat)
	unsigned int nlat_padded;	///< number of spatial points in Theta direction, including padding.
	unsigned int nlm_cplx;		///< number of complex coefficients to represent a complex-valued spatial field.
/* END OF PUBLIC PART */

	short fftc_mode;			///< how to perform the complex fft : -1 = no fft; 0 = interleaved/native; 1 = split/transpose.
	unsigned short nthreads;	///< number of threads (openmp).
	unsigned short *tm;			///< start theta value for SH (polar optimization : near the poles the legendre polynomials go to zero for high m's)
	short robert_form;			///< flag for Robert formulation: if true, the vector synthesis are multiplied by sin(theta) and the analysis are divided by sin(theta).
	int k_stride_a;				///< stride in theta direction
	int m_stride_a;				///< stride in phi direction in intermediate spectral space (m)
	double *wg;					///< Gauss weights for Gauss-Legendre quadrature.
	double *st_1;				///< 1/sin(theta);
	double mpos_scale_analys;	///< scale factor for analysis, handles real-norm (0.5 or 1.0);

	fftw_plan ifftc, fftc;
	fftw_plan ifft_cplx, fft_cplx;		// for complex-valued spatial fields.
	fftw_plan ifftc_block, fftc_block;

	/* Legendre function generation arrays */
	double *alm;	// coefficient list for Legendre function recurrence (size 2*NLM)
	double *blm;	// coefficient list for modified Legendre function recurrence for analysis (size 2*NLM)
	double *l_2;	// array of size (LMAX+1) containing 1./l(l+1) for increasing integer l.
	/* matrices for vector transform (to convert to scalar transforms) */
	double *mx_stdt;	// sparse matrix for  sin(theta).d/dtheta,  couples l-1 and l+1
	double *mx_van;		// sparse matrix for  sin(theta).d/dtheta + 2*cos(theta),  couples l-1 and l+1
	#ifdef SHTNS_ISHIOKA
	/* for the new recurrence of Ishioka */
	double *clm;	// a_lm, b_lm
	double *xlm;	// epsilon_lm * alpha_lm
	double *x2lm;	// epsilon_lm * alpha_lm, scaled for analysis (different from xlm only with Schmidt semi-normalization)
	#endif

	void* ftable[SHT_NVAR][SHT_NTYP];		// pointers to transform functions.

	/* rotation stuff (pseudo-spectral) */
	unsigned npts_rot;		// number of physical points needed
	fftw_plan fft_rot;		// Fourier transform for rotations
	double* ct_rot;			// cos(theta) array
	double* st_rot;			// sin(theta) array

	/* _to_lat stuff */
	double* ylm_lat;
	double ct_lat;
	fftw_plan ifft_lat;		///< fftw plan for SHqst_to_lat
	int nphi_lat;			///< nphi of previous SHqst_to_lat

	#ifdef HAVE_LIBCUFFT
	/* cuda stuff */
	short cu_flags;
	short cu_fft_mode;
	double* d_alm;
	double* d_ct;
	double* d_mx_stdt;
	double* d_mx_van;
	double* gpu_mem;
	double* xfft;
	double* xfft_cpu;
	size_t nlm_stride, spat_stride;
	cudaStream_t xfer_stream, comp_stream;		// the cuda streams
	cufftHandle cufft_plan;						// the cufft Handle
	#endif

	/* other misc informations */
	unsigned char nlorder;	// order of non-linear terms to be resolved by SH transform.
	unsigned char grid;		// store grid type.
	short norm;				// store the normalization of the Spherical Harmonics (enum \ref shtns_norm + \ref SHT_NO_CS_PHASE flag)
	unsigned fftw_plan_mode;
	unsigned layout;		// requested data layout
	double Y00_1, Y10_ct, Y11_st;
	shtns_cfg next;		// pointer to next sht_setup or NULL (records a chained list of SHT setup).
	// the end should be aligned on the size of int, to allow the storage of small arrays.
};

// define shortcuts to sizes.
#define NLM shtns->nlm
#define LMAX shtns->lmax
#define NLAT shtns->nlat
#define NLAT_2 shtns->nlat_2
#define NPHI shtns->nphi
#define MMAX shtns->mmax
#define MRES shtns->mres
#define SHT_NL_ORDER shtns->nlorder

// define index in alm/blm matrices
#define ALM_IDX(shtns, im) ( (im)*(2*(shtns->lmax+1) - ((im)-1)*shtns->mres) )

// SHT_NORM without CS_PHASE
#define SHT_NORM (shtns->norm & 0x0FF)

#ifndef M_PI
# define M_PI 3.1415926535897932384626433832795
#endif
#ifndef M_PIl
# define M_PIl 3.1415926535897932384626433832795L
#endif

// value for on-the-fly transforms is lower because it allows to optimize some more (don't compute l which are not significant).
#define SHT_L_RESCALE_FLY 1000
// set to a value close to the machine accuracy, it allows to speed-up on-the-fly SHTs with very large l (lmax > SHT_L_RESCALE_FLY).
#define SHT_ACCURACY 1.0e-20
// scale factor for extended range numbers (used in on-the-fly transforms to compute recurrence)
#define SHT_SCALE_FACTOR 2.9073548971824275622e+135
//#define SHT_SCALE_FACTOR 2.0370359763344860863e+90

#ifdef __NVCC__
		// disable vector extensions when compiling cuda code.
        #undef _GCC_VEC_
#endif

#include "shtns_simd.h"

#if VSIZE2 >= 8
	// these values must be adjusted for the larger AVX512 vectors
	#undef SHT_L_RESCALE_FLY
	#undef SHT_ACCURACY
	#define SHT_L_RESCALE_FLY 1800
	#define SHT_ACCURACY 1.0e-40
#endif

struct DtDp {		// theta and phi derivatives stored together.
	double t, p;
};

#define GLUE2(a,b) a##b
#define GLUE3(a,b,c) a##b##c

// verbose printing
#if SHT_VERBOSE > 1
  #define PRINT_VERB(msg) printf(msg)
#else
  #define PRINT_VERB(msg) (0)
#endif

/// Convert from vector 2 scalar SH to vector SH
/// Slm = - (I*m*Wlm + MX*Vlm) / (l*(l+1))		=> why does this work ??? (aliasing of 1/sin(theta) ???)
/// Tlm = - (I*m*Vlm - MX*Wlm) / (l*(l+1))
/// m = signed m (for complex SH transform).
/// double* mx = shtns->mx_van + 2*LM(shtns,m,m);	//(im*(2*(LMAX+1)-(m+MRES))) + 2*m;
static void SH_2scal_to_vect_reduce(const double *mx, const double* l_2, int llim, int m, rnd* vw, v2d* Sl, v2d* Tl)
{
	double em = m;
	m = abs(m);
	v2d vl = v2d_reduce(vw[0], vw[1]);
	v2d wl = v2d_reduce(vw[2], vw[3]);
	v2d sl1 = vdup( 0.0 );
	v2d tl1 = vdup( 0.0 );
	for (int l=0; l<=llim-m; l++) {
		s2d mxu = vdup( mx[2*l] );
		s2d mxl = vdup( mx[2*l+1] );		// mxl for next iteration
		v2d sl = sl1 + IxKxZ(em, wl);		// sl1 + I*em*wl
		v2d tl = tl1 + IxKxZ(em, vl);		// tl1 + I*em*vl
		sl1 =  mxl*vl;			// vs for next iter
		tl1 = -mxl*wl;			// wt for next iter
		vl = v2d_reduce(vw[4*l+4], vw[4*l+5]);		// kept for next iteration
		wl = v2d_reduce(vw[4*l+6], vw[4*l+7]);
		sl += mxu*vl;
		tl -= mxu*wl;
		Sl[l] = -sl * vdup(l_2[l+m]);
		Tl[l] = -tl * vdup(l_2[l+m]);
	}
}

/// Convert from vector 2 scalar SH to vector SH
/// Slm = - (I*m*Wlm + MX*Vlm) / (l*(l+1))		=> why does this work ??? (aliasing of 1/sin(theta) ???)
/// Tlm = - (I*m*Vlm - MX*Wlm) / (l*(l+1))
/// m = signed m (for complex SH transform).
static void SH_2scal_to_vect(const double *mx, const double* l_2, int llim, int m, v2d* vw, v2d* Sl, v2d* Tl)
{
#if !defined( _GCC_VEC_) || !defined( __AVX__ )
	double em = m;
	m = abs(m);
	v2d vl = vw[0];
	v2d wl = vw[1];
	v2d sl1 = vdup( 0.0 );
	v2d tl1 = vdup( 0.0 );
	for (int l=0; l<=llim-m; l++) {
		s2d mxu = vdup( mx[2*l] );
		s2d mxl = vdup( mx[2*l+1] );		// mxl for next iteration
		v2d sl = sl1 + IxKxZ(em, wl);	// sl1 + I*em*wl;
		v2d tl = tl1 - IxKxZ(em, vl);	// sl1 - I*em*vl;
		sl1 =  mxl*vl;			// vs for next iter
		tl1 =  mxl*wl;			// wt for next iter
		vl = vw[2*l+2];		// kept for next iteration
		wl = vw[2*l+3];
		sl += mxu*vl;
		tl += mxu*wl;
		Sl[l] = -sl * vdup(l_2[l+m]);
		Tl[l] = tl * vdup(l_2[l+m]);
	}
  #else
	v4d em = _mm256_setr_pd(-m,m, m,-m);
	m = abs(m);
	v4d vwl = vread4(vw, 0);
	v4d stl = em * vreverse4(vwl);
	for (int l=0; l<=llim-m; l++) {
			// SH_vect_to_2scal :: 2 full permutes, 1 mul, 2 fma, 2 128-bit stores, 2 64-bit broadcasts, 1 256-bit load
			// here :: 1 full permutes, 2 mul, 2 fma, 2 128bit-stores, 3 64-bit broadcasts, 1 256-bit load
		v4d vwu = vread4(vw+2*l+2, 0);		// kept for next iteration
		v4d mxu = vall4( mx[2*l] );
		stl += mxu * vwu;
		stl *= vall4(l_2[l+m]);
		Sl[l] = - (v2d) _mm256_castpd256_pd128(stl);
		Tl[l] = _mm256_extractf128_pd(stl, 1);
		stl = em * vreverse4(vwu) + vwl * vall4( mx[2*l+1] );
		vwl = vwu;
	}
  #endif
}

/// post-processing for recurrence relation of Ishioka
/// xlm = shtns->xlm + 3*im*(2*(LMAX+4) -m+MRES)/4;
/// llim_m = llim-m
/// qq[l-m]: input data obtained with ishioka's relation
/// Ql[l-m]: output data, spherical harmonic coefficients of degree l (for fixed m).
/// can operate in-place (Ql = qq)
static void ishioka_to_SH(const double* xlm, const v2d* qq, const int llim_m, v2d* Ql)
{
	long l=0;	long ll=0;
  #if !defined( _GCC_VEC_) || !defined( __AVX__ )
	v2d u0 = vdup(0.0);
	while (l<llim_m) {
		v2d uu = qq[l];
		Ql[l] = uu * vdup(xlm[ll]) + u0;
		Ql[l+1] = qq[l+1] * vdup(xlm[ll+2]);
		u0 = uu * vdup(xlm[ll+1]);
		l+=2;	ll+=3;
	}
	if (l==llim_m) {
		Ql[l] = qq[l] * vdup(xlm[ll]) + u0;
	}
  #else
	v4d z = vall4(0.0);
//	#pragma GCC unroll 2
	while (l<llim_m) {
		v4d x = vread4(xlm+ll, 0);	x = vdup_even4(x);
		v4d uu = vread4(qq+l, 0);	// [qq[l], qq[l+1]]
		vstor4(Ql+l, 0, uu*x + z);
		z = _mm256_castpd128_pd256( (v2d)_mm256_castpd256_pd128(uu) * vdup(xlm[ll+1]) );		// upper part of z is zeroed.
		l+=2;	ll+=3;
	}
	if (l==llim_m) {
		Ql[l] = qq[l] * vdup(xlm[ll]) + (v2d)_mm256_castpd256_pd128(z);
	}
  #endif
}

/*
qnew[l]   = q[l]*a + q[l-2]*b;
qnew[l+1] = q[l+1] * c;

AVX: q1 =  q[l, l+1]
	 q0 =  q[l-2, l-1]
qnew = q1*[a,a,c,c] + q0*[b,b,0,0]   [1x, 1fma, 2 shuffles] => not worth it

AVX512: q1 = q[l, l+1, l+2, l+3]
		q0 = q[l-2, l-1, l, l+1]  => [additional read!]
		read: [a,b,c,a',b'c',x,x] => 2 shuffles (7 cycles)
qnew = q1*[a,a,c,c,a',a',c',c'] + q0*[b,b,0,0,b',b',0,0]   => [1x,1fma,2shuffle+1read] for 4 cplx.
		use _mm512_maskz_permutexvar_pd() to produce second coeff including zeros.
*/

/// Same as \ref ishioka_to_SH, but for two interlaced coefficient lists
/// set llim_m = llim-m+1 for vector transforms that include llim+1 (before post-processing)
/// can operate in-place (vw = VWl)
static void ishioka_to_SH2(const double* xlm, const v2d* vw, const int llim_m, v2d* VWl)
{
	long l=0;	long ll=0;
  #if !defined( _GCC_VEC_) || !defined( __AVX__ )
	v2d v0 = vdup(0.0);
	v2d w0 = vdup(0.0);
	while (l<llim_m) {
		v2d vv = vw[2*l];
		v2d ww = vw[2*l+1];
		VWl[2*l]   = vv * vdup(xlm[ll]) + v0;
		VWl[2*l+1] = ww * vdup(xlm[ll]) + w0;
		VWl[2*l+2] = vdup(xlm[ll+2]) * vw[2*l+2];
		VWl[2*l+3] = vdup(xlm[ll+2]) * vw[2*l+3];
		v0 = vv * vdup(xlm[ll+1]);
		w0 = ww * vdup(xlm[ll+1]);
		l+=2;	ll+=3;
	}
	if (l==llim_m) {
		v2d vv = vw[2*l];
		v2d ww = vw[2*l+1];
		VWl[2*l]   = vv * vdup(xlm[ll]) + v0;
		VWl[2*l+1] = ww * vdup(xlm[ll]) + w0;
	}
  #else
	v4d vw0 = vall4(0.0);
	while (l<llim_m) {
		v4d vwl = vread4( vw + 2*l, 0);
		vw0 += vwl * vall4(xlm[ll]);
		vstor4(VWl + 2*l, 0, vw0);
		v4d y = vread4(vw + 2*l, 1) * vall4(xlm[ll+2]);
		vstor4(VWl + 2*l, 1, y);
		vw0 = vwl * vall4(xlm[ll+1]);
		l+=2;	ll+=3;
	}
	if (l==llim_m) {
		v4d vwl = vread4( vw + 2*l, 0);
		vw0 += vwl * vall4(xlm[ll]);
		vstor4(VWl + 2*l, 0, vw0);
	}
  #endif
}

/// pre-processing for recurrence relation of Ishioka
/// xlm = shtns->xlm + 3*im*(2*(LMAX+4) -m+MRES)/4;
/// llim_m = llim-m
/// Ql[l-m]: intput data, spherical harmonic coefficients of degree l (for fixed m).
/// ql[l-m]: output data, ready for ishioka's recurrence
/// can operate in-place (Ql = qq)
static void SH_to_ishioka(const double* xlm, const v2d* Ql, const int llim_m, v2d* ql)
{
	long l=0;	long ll=0;
  #if !defined( _GCC_VEC_) || !defined( __AVX__ )
	v2d qq = Ql[l] * vdup(xlm[0]);
	while (l<llim_m-1) {
		v2d qq2 = Ql[l+2];
		ql[l]   = (qq  +  qq2 * vdup(xlm[ll+1]));
		ql[l+1] = Ql[l+1] * vdup(xlm[ll+2]);
		ll+=3;	l+=2;
		qq = qq2 * vdup(xlm[ll]);
	}
	ql[l]   = qq;
    qq = vdup(0.0);
    if (l<llim_m)  qq = Ql[l+1] * vdup(xlm[ll+2]);
    ql[l+1] = qq;
  #else
	v4d y = vread4(Ql+l, 0);
/*	while (l<llim_m-4) {
		v4d x = vread4(xlm+ll, 0);	x = vdup_even4(x);
		v4d y2 = vread4(Ql+l+2, 0);
		v4d z =_mm256_castpd128_pd256( vdup(xlm[ll+1]) * _mm256_castpd256_pd128( y2 ) );		// upper part of z is zeroed.
		vstor4(ql+l, 0, x*y + z);
		x = vread4(xlm+ll+3, 0);	x = vdup_even4(x);
		y = vread4(Ql+l+4, 0);
		z =_mm256_castpd128_pd256( vdup(xlm[ll+4]) * _mm256_castpd256_pd128( y ) );		// upper part of z is zeroed.
		vstor4(ql+l+2, 0, x*y2 + z);
		ll+=6;	l+=4;
	}	*/
	while (l<llim_m-1) {
		v4d x = vread4(xlm+ll, 0);	x = vdup_even4(x);
		v4d y = vread4(Ql+l, 0);
		v4d z =_mm256_castpd128_pd256( vdup(xlm[ll+1]) * Ql[l+2] );		// upper part of z is zeroed with AVX
		vstor4(ql+l, 0, x*y + z);
		ll+=3;	l+=2;
	}
	ql[l]   = Ql[l] * vdup(xlm[ll]);
    v2d qq = vdup(0.0);
    if (l<llim_m)  qq = Ql[l+1] * vdup(xlm[ll+2]);
    ql[l+1] = qq;
  #endif
}

/// same as \ref SH_to_ishioka, but handles two interleaved arrays + operates in-place.
/// use llim_m = llim-m+1 for vector datat that goes up to llim+1
static void SH2_to_ishioka(const double* xlm, v2d* VWl, const int llim_m)
{
	long l=0;	long ll=0;
  #if !defined( _GCC_VEC_) || !defined( __AVX__ )
	v2d vv = VWl[2*l]   * vdup(xlm[0]);
	v2d ww = VWl[2*l+1] * vdup(xlm[0]);
	while (l<llim_m-1) {
		v2d vv2 = VWl[2*(l+2)];
		v2d ww2 = VWl[2*(l+2)+1];
		VWl[2*l]   = (vv  +  vv2 * vdup(xlm[ll+1]));
		VWl[2*l+1] = (ww  +  ww2 * vdup(xlm[ll+1]));
		VWl[2*l+2] *= vdup(xlm[ll+2]);
		VWl[2*l+3] *= vdup(xlm[ll+2]);
		ll+=3;	l+=2;
		vv = vv2 * vdup(xlm[ll]);
		ww = ww2 * vdup(xlm[ll]);
	}
	VWl[2*l]   = vv;
	VWl[2*l+1] = ww;
	if (l<=llim_m-1) {
		VWl[2*l+2] *= vdup(xlm[ll+2]);
		VWl[2*l+3] *= vdup(xlm[ll+2]);
	}
  #else
	v4d vwl = vread4(VWl +2*l, 0) * vall4(xlm[0]);
	while (l<llim_m-1) {
		v4d vw2 = vread4(VWl +2*l, 2);
		vwl += vw2 * vall4(xlm[ll+1]);
		vstor4(VWl +2*l, 0, vwl);
		v4d y = vread4(VWl + 2*l, 1) * vall4(xlm[ll+2]);
		vstor4(VWl + 2*l, 1, y);
		ll+=3;	l+=2;
		vwl = vw2 * vall4(xlm[ll]);
	}
	vstor4(VWl + 2*l, 0, vwl);
	if (l<=llim_m-1) {
		v4d y = vread4(VWl + 2*l, 1) * vall4(xlm[ll+2]);
		vstor4(VWl + 2*l, 1, y);
	}
  #endif
}

/// Convert from vector SH to 2 scalar SH
/// Vlm =  st*d(Slm)/dtheta + I*m*Tlm
/// Wlm = -st*d(Tlm)/dtheta + I*m*Slm
/// store interleaved: VWlm(2*l) = Vlm(l);	VWlm(2*l+1) = Wlm(l);
/// m = signed m (for complex SH transform).
static void SH_vect_to_2scal(const double *mx, int llim, int m, cplx* Sl, cplx* Tl, cplx* VWl)
{
	long l;
  #if !defined(_GCC_VEC_) || !defined( __AVX__ )
	double em = m;
	v2d sl = ((v2d*)Sl)[m];
	v2d tl = ((v2d*)Tl)[m];
	v2d vs = IxKxZ(em, tl);
	v2d wt = IxKxZ(em, sl);
	for (l=m; l<llim; l++) {
		v2d sl1 = ((v2d*)Sl)[l+1];		// kept for next iteration
		v2d tl1 = ((v2d*)Tl)[l+1];
		s2d mxu = vdup(mx[2*l]);
		s2d mxl = vdup(mx[2*l+1]);	// mxl for next iteration
		((v2d*)VWl)[2*l]   = vs + mxu*sl1;
		((v2d*)VWl)[2*l+1] = wt - mxu*tl1;
		vs = IxKxZ(em, tl1) + mxl*sl;		// vs += I*em*tl;
		wt = IxKxZ(em, sl1) - mxl*tl;		// wt += I*em*sl;
		sl = sl1;
		tl = tl1;
	}
	//if (l==llim)		// Because m<=llim, this is always true.
	{
		s2d mxl = vdup(mx[2*l+1]);		// mxl for next iteration
		((v2d*)VWl)[2*l]   = vs;
		((v2d*)VWl)[2*l+1] = wt;
		((v2d*)VWl)[2*llim+2] =  mxl*sl;
		((v2d*)VWl)[2*llim+3] = -mxl*tl;
	}
  #else
	v4d em = (v4d) _mm256_setr_pd(-m,m, m,-m);
	v4d stl = v2d_x2_to_v4d( -((v2d*)Sl)[m], ((v2d*)Tl)[m]);
	v4d vswt = em*vreverse4(stl);
	for (l=m; l<llim; l++) {
		// 2 full permutes, 1 mul, 2 fma, 2 128bit-loads, 2 64-bit broadcasts, 1 256-bit store
		v4d stlu = v2d_x2_to_v4d( -((v2d*)Sl)[l+1], ((v2d*)Tl)[l+1]);
		v4d mxu = vall4(mx[2*l]);
		v4d mxl = vall4(mx[2*l+1]);		// mxl for next iteration
		vswt -= mxu*stlu;
		vstor4(VWl +2*l, 0, vswt);
		vswt = em*vreverse4(stlu) - mxl*stl;
		stl = stlu;			// kept for next iter
	}
	//if (l==llim)		// Because m<=llim, this is always true.
	{
		v4d mxl = vall4(mx[2*l+1]);		// mxl for next iteration
		vstor4(VWl +2*l, 0, vswt);
		vstor4(VWl + 2*l+2, 0, -mxl*stl);
	}
  #endif
}

static void SH_vect_to_2scal_alt(const double *mx, int llim, int m, const cplx* Sl, const cplx* Tl, cplx* VWl)
{
	double em = m;
	#ifdef _GCC_VEC_
	const rnd emx = vneg_even_precalc( vall(em) );
	//const rnd emx = _mm256_setr_pd(-em, em, -em, em);
	#endif
	long l=m;
	{
		s2d mxu = vdup(mx[2*l]);
		#ifndef _GCC_VEC_
		v2d s = I*em*Tl[l];
		v2d t = I*em*Sl[l];
		#else
		v2d s = v2d_lo(emx) * vxchg(((v2d*)Tl)[l]);
		v2d t = v2d_lo(emx) * vxchg(((v2d*)Sl)[l]);
		#endif
		if (l<llim) {
			s += mxu*((v2d*)Sl)[l+1];
			t -= mxu*((v2d*)Tl)[l+1];
		}
		((v2d*)VWl)[2*l]   = s;
		((v2d*)VWl)[2*l+1] = t;
		l++;
	}
	#if VSIZE2 >= 4
	#if VSIZE2 == 4
		// AVX: there can be some data forwarding.
		rnd Sll = vread(Sl+l-1, 0);
		rnd Tll = vread(Tl+l-1, 0);
	#endif
	for (; l<=llim-VSIZE2/2; l+=VSIZE2/2) {		// general case 		V[2*l] = mx[2*l-1]*S[l-1]
		// AVX512: 4 in-lane permutes, 2 full permutes, 2 mul, 4 fma, 7 512-bit loads, 2 512-bit stores
		// AVX: 4 in-lane permutes, 2 full permutes, 2 mul, 4 fma, 5 256-bit loads, 2 256-bit stores
		rnd s = emx * vxchg_even_odd( vread(Tl+l,0) );
		rnd t = emx * vxchg_even_odd( vread(Sl+l,0) );
		rnd mxx = vread(mx+2*l-1, 0);
		rnd mxl = vdup_even(mxx);
		rnd mxu = vdup_odd(mxx);
		#if VSIZE2 == 4
			rnd Slu = vread(Sl+l+1, 0);
			rnd Tlu = vread(Tl+l+1, 0);
			s += mxl * Sll + mxu * Slu;
			t -= mxl * Tll + mxu * Tlu;
			Sll = Slu;		Tll = Tlu;		// kept for next iteration
		#else
			s += mxl * vread(Sl+l-1, 0) + mxu * vread(Sl+l+1, 0);
			t -= mxl * vread(Tl+l-1, 0) + mxu * vread(Tl+l+1, 0);
		#endif
		#ifdef __AVX512F__
			vstor(VWl+2*l, 0, _mm512_permutex2var_pd(s,_mm512_setr_epi64(0,1,8,9,2,3,10,11), t) );
			vstor(VWl+2*l, 1, _mm512_permutex2var_pd(s,_mm512_setr_epi64(4,5,12,13,6,7,14,15), t) );
		#elif defined( __AVX__ )
			vstor(VWl+2*l, 0, _mm256_permute2f128_pd(s,t, 0x20) );
			vstor(VWl+2*l, 1, _mm256_permute2f128_pd(s,t, 0x31) );
		#else
			#error "unsupported simd vectors"
		#endif
	}
	#endif
	for (; l<=llim; l++) {		// general case, reminder 		V[2*l] = mx[2*l-1]*S[l-1]
		s2d mxl = vdup(mx[2*l-1]);
		s2d mxu = vdup(mx[2*l]);
		#ifndef _GCC_VEC_
		v2d imt = I*em*Tl[l];
		v2d ims = I*em*Sl[l];
		#else
		v2d imt = v2d_lo(emx) * vxchg(((v2d*)Tl)[l]);
		v2d ims = v2d_lo(emx) * vxchg(((v2d*)Sl)[l]);
		#endif
		if (l<llim) {
			imt += mxu*((v2d*)Sl)[l+1];
			ims -= mxu*((v2d*)Tl)[l+1];
		}
		((v2d*)VWl)[2*l]   = imt + mxl*((v2d*)Sl)[l-1];
		((v2d*)VWl)[2*l+1] = ims - mxl*((v2d*)Tl)[l-1];
	}
	{	//l=llim+1
		s2d mxl = vdup(mx[2*l-1]);
		((v2d*)VWl)[2*l] = mxl * ((v2d*)Sl)[l-1];
		((v2d*)VWl)[2*l+1] = -mxl * ((v2d*)Tl)[l-1];		
	}
}


static void SHsph_to_2scal(const double *mx, int llim, int m, cplx* Sl, cplx* VWl)
{
	double em = m;
	v2d sl = ((v2d*)Sl)[m];
	v2d vs = vdup(0.0);
	v2d wt = IxKxZ(em, sl);
	long l;
	for (l=m; l<llim; l++) {
		v2d sl1 = ((v2d*)Sl)[l+1];
		s2d mxu = vdup(mx[2*l]);
		s2d mxl = vdup(mx[2*l+1]);		// mxl for next iteration
		((v2d*)VWl)[2*l]   = vs + mxu*sl1;
		((v2d*)VWl)[2*l+1] = wt;	// IxKxZ(em, sl);
		vs = mxl*sl;			// vs for next iter
		wt = IxKxZ(em, sl1);
		sl = sl1;		// kept for next iteration
	}
	//if (l==llim)		// Because m<=llim, this is always true.
	{
		s2d mxl = vdup(mx[2*l+1]);		// mxl for next iteration
		((v2d*)VWl)[2*l]   = vs;
		((v2d*)VWl)[2*l+1] = wt;
		((v2d*)VWl)[2*l+2] = mxl*sl;
		((v2d*)VWl)[2*l+3] = vdup(0.0);
	}
}

static void SHtor_to_2scal(const double *mx, int llim, int m, cplx* Tl, cplx* VWl)
{
	double em = -m;
	v2d tl = - ((v2d*)Tl)[m];
	v2d vs = IxKxZ(em, tl);
	v2d wt = vdup(0.0);
	long l;
	for (l=m; l<llim; l++) {
		v2d tl1 = - ((v2d*)Tl)[l+1];
		s2d mxu = vdup(mx[2*l]);
		s2d mxl = vdup(mx[2*l+1]);		// mxl for next iteration
		((v2d*)VWl)[2*l]   = vs;
		((v2d*)VWl)[2*l+1] = wt + mxu*tl1;
		wt = mxl*tl;			// wt for next iter
		vs = IxKxZ(em, tl1);
		tl = tl1;
	}
	//if (l==llim)		// Because m<=llim, this is always true.
	{
		s2d mxl = vdup(mx[2*l+1]);		// mxl for next iteration
		((v2d*)VWl)[2*l]   = vs;
		((v2d*)VWl)[2*l+1] = wt;
		((v2d*)VWl)[2*l+2] = vdup(0.0);
		((v2d*)VWl)[2*l+3] = mxl*tl;
	}
}

static void zero_poles4_vect(v2d* F0, long ofsm, long ofs1, long n) {
	#pragma omp simd
	for (long i=0; i<n*VSIZE2; i++) {
		((double*)F0)[i] = 0.0;
		((double*)(F0+ofs1))[i] = 0.0;
		((double*)(F0+ofsm))[i] = 0.0;
		((double*)(F0+ofsm+ofs1))[i] = 0.0;
	}
}

static void zero_poles2_vect(v2d* F0, long ofsm, long n) {
	#pragma omp simd
	for (long i=0; i<n*VSIZE2; i++) {
		((double*)F0)[i] = 0.0;
		((double*)(F0+ofsm))[i] = 0.0;
	}
}

static void zero_mem(v2d* F0, long n) {
	#pragma omp simd
	for (long i=0; i<n*VSIZE2; i++) {
		((double*)F0)[i] = 0.0;
	}
}
