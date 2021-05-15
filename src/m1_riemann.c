/*
 * Copyright (c) 2017-2020 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * M1-riemann-solver is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#include <math.h>

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

#include "m1_riemann.h"

#define M1_R_ABS_MAX_VALUE          0.99

#define GSL_DERIV_STEP              1E-6
#define GSL_BRENT_MAX_ITER          100
#define GSL_BRENT_ABS_TOL           1E-5
#define GSL_QAGS_MAX_ITERVALS       1000
#define GSL_QAGS_EPS_REL            1E-4
#define GSL_QAGS_EPS_ABS            1E-4

/* Packed arguments structs for GSL function wrapper */
struct args_pack2 {
	int i;
	double xi;
};

struct args_pack4 {
	double rhoL;
	double rhoR;
	double rL;
	double rR;
};

static double get_chi_gsl(const double r, void *args);
static double get_chi_deriv(const double r, const int type);
static void get_eigen_values(const double r, double *eigen1, double *eigen2);
static double get_eigen_value_i(const int i, const double r);
static double get_shock(const int i12, const double rm, const double rp);
static double get_li(const int i, const double rm, const double rp);
static double get_rarefaction(const int i12, const double rm, const double rp);
static double get_func2min_rM_brent(const double x, void *args);
static double get_rM_brent(const double uL[2], const double uR[2]);
static double get_func2min_r_brent(const double x, void *args);
static double get_r_brent(const int i, const double xi);
static double get_func2int_invariant_qags(const double r, void *args);
static double get_invariant_qags(int i, const double rm, const double rp);

/*
 * GSL wrapper (with arguments packed in *args) that computes the value of 
 * the Eddington factor.
 */
double get_chi_gsl(const double r, void *args)
{
	return get_chi(r);
}

/*
 * Compute the derivative of the Eddington factor X'(r) using the GSL lib.
 */
static double get_chi_deriv(const double r, const int type)
{
	double res = 0;
	double abserr;
	gsl_function F;

	F.function = &get_chi_gsl;
	F.params = 0;

	switch (type) {
	case -1:
		gsl_deriv_backward(&F, r, GSL_DERIV_STEP, &res, &abserr);
		break;

	case 0:
		gsl_deriv_central(&F, r, GSL_DERIV_STEP, &res, &abserr);
		break;

	case 1:
		gsl_deriv_forward(&F, r, GSL_DERIV_STEP, &res, &abserr);
		break;
	}
	return res;
}

/*
 * Compute both eigen values (vp1 and vp2) in function of r.
 */
static void get_eigen_values(const double r, double *eigen1, double *eigen2)
{
	const double chi = get_chi(r);
	const double dchi = get_chi_deriv(r, 0);
	const double sqdelta_2 = 0.5 * sqrt(dchi * dchi - 4 * r * dchi + 4 * chi);

	*eigen1 = 0.5 * dchi - sqdelta_2;
	*eigen2 = 0.5 * dchi + sqdelta_2;
}

/*
 * Return the i-th (i=1,2) eigen value in function of r.
 */
static double get_eigen_value_i(const int i, const double r)
{
	const double chi = get_chi(r);
	const double dchi = get_chi_deriv(r, 0);
	const double sqdelta_2 = 0.5 * sqrt(dchi * dchi - 4 * r * dchi + 4 * chi);

	if (i == 1) {
		return 0.5 * dchi - sqdelta_2;
	} else {
		return 0.5 * dchi + sqdelta_2;
	}
}

/*
 * Compute the Li quantity (see article quoted at the top of this file).
 */
static double get_shock(const int i12, const double rm, const double rp)
{
	const double chip = get_chi(rp);
	const double chim = get_chi(rm);
	const double a = chip - rp * rp;
	const double b = -chip - chim + 2 * rp * rm;
	const double c = chim - rm * rm;
	const double d = b * b - 4 * a * c;

	if (i12 == 1) {
		return (-b + sqrt(d)) / (2 * a);
	} else {
		return (-b - sqrt(d)) / (2 * a);
	}
}

/*
 * Compute the L_i quantity (see article quoted at the top of this file).
 */
static double get_li(const int i, const double rm, const double rp)
{
	if (rp < rm) {
		return get_shock(i, rm, rp);
	} else {
		return get_rarefaction(i, rm, rp);
	}
}

static double get_rarefaction(const int i12, const double rm, const double rp)
{
	return exp(get_invariant_qags(i12, rm, rp));
}

/*
 * GSL wrapper (with arguments packed in *args) that computes the value of 
 * the function to minimize for finding the value rM associated to the 
 * middle state.
 */
static double get_func2min_rM_brent(const double x, void *args)
{
	/* Unpacking arguments */
	struct args_pack4 *p = (struct args_pack4 *)args;
	return (p->rhoR / p->rhoL) - get_li(1, p->rL, x) * get_li(2, x, p->rR);
}

/*
 * Returns the value of rM (middle state) for given left/right states uL/uR.
 * We use a Brent method from the GSL librairie for finding rM value. 
 * This algorithm ensures that all intermediates values of rM always stay 
 * in [-1, 1].
 */
static double get_rM_brent(const double uL[2], const double uR[2])
{
	const gsl_root_fsolver_type *t;
	gsl_root_fsolver *s;
	gsl_function f;

	/* Packing arguments rhoL, rhoR, rL, rR*/
	struct args_pack4 args = {
		.rhoL = uL[0], .rhoR = uR[0], .rL = uL[1] / uL[0], .rR = uR[1] / uR[0]
	};

	f.function = &get_func2min_rM_brent;
	f.params = &args;
	t = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(t);

	/* Initial searching interval */
	double x_lo = -M1_R_ABS_MAX_VALUE;
	double x_hi = M1_R_ABS_MAX_VALUE;

	gsl_root_fsolver_set(s, &f, x_lo, x_hi);

	int status = GSL_CONTINUE;
	int iter = 0;
	double res = 0;

	while (status == GSL_CONTINUE && iter < GSL_BRENT_MAX_ITER) {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		res = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi, 0, GSL_BRENT_ABS_TOL);
	}

	gsl_root_fsolver_free(s);
	return res;
}

/*
 * GSL wrapper (with arguments packed in *args) that computes the value of the 
 * function to minimize for finding the value r.
 */
static double get_func2min_r_brent(const double x, void *args)
{
	/* Unpacking arguments */
	struct args_pack2 *p = (struct args_pack2 *)args;

	return p->xi - get_eigen_value_i(p->i, x);
}

/*
 * Returns the value of r for given xi = x/t and i (i=1,2). We use a Brent 
 * method from the GSL librairie for finding r value. This algorithm ensures 
 * that all intermediates values of rM always stay in [-1, 1].
 */
static double get_r_brent(const int i, const double xi)
{
	const gsl_root_fsolver_type *t;
	gsl_root_fsolver *s;
	gsl_function f;

	f.function = &get_func2min_r_brent;
	/* Packing arguments */
	struct args_pack2 args = { .i = i, .xi = xi };
	f.params = &args;

	t = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(t);

	/* Initial searching interval */
	double x_lo = -M1_R_ABS_MAX_VALUE;
	double x_hi = M1_R_ABS_MAX_VALUE;

	gsl_root_fsolver_set(s, &f, x_lo, x_hi);

	int status = GSL_CONTINUE;
	int iter = 0;
	double res = 0;

	while (status == GSL_CONTINUE && iter < GSL_BRENT_MAX_ITER) {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		res = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi, 0, GSL_BRENT_ABS_TOL);
	}

	gsl_root_fsolver_free(s);
	return res;
}

/*
 * GSL wrapper (with arguments packed in *args) that computes the value of the 
 * function to integrate for computing the Riemann invariant.
 */
static double get_func2int_invariant_qags(const double r, void *args)
{
	/* Unpacking arguments */
	const int i = *(int *)args;
	return 1.0 / (get_eigen_value_i(i, r) - r);
}

/*
 * Returns the value of the Riemann invariant for a given i, rm and rp values.
 * To perform the numerical integration we use a QAGS method from the GSL 
 * librairie for finding r value.
 */
static double get_invariant_qags(int i, const double rm, const double rp)
{
	double res;
	double err;

	gsl_integration_workspace *w =
		gsl_integration_workspace_alloc(GSL_QAGS_MAX_ITERVALS);

	gsl_function f;

	f.function = &get_func2int_invariant_qags;

	/* Packing arguments */
	f.params = &i;

	gsl_integration_qags(&f, rm, rp, GSL_QAGS_EPS_ABS, GSL_QAGS_EPS_REL,
						 GSL_QAGS_MAX_ITERVALS, w, &res, &err);

	gsl_integration_workspace_free(w);
	return res;
}

/*
 * Solve the exact one-dimensional riemann problem for the M1 model. 
 * The solution u is computed, given :
 *  - xi : x / t
 *  - uL : left state
 *  - uR : right state
 */
void solve_riemann(const double xi, const double uL[2], const double uR[2],
				   double u[2])
{
	double r;
	double vp1L, vp1R, vp2L, vp2R;
	double vp1Mid, vp2Mid;
	double vp1m, vp2m, vp1p, vp2p;
	double x1, x2;

	double rL = uL[1] / uL[0];
	double rR = uR[1] / uR[0];

	/* Compute middle state */
	double rM = get_rM_brent(uL, uR);
	double rhoM = get_li(1, rL, rM) * uL[0];

	/* Compute characteristics associated to each state  */
	get_eigen_values(rL, &vp1L, &vp2L);
	get_eigen_values(rM, &vp1Mid, &vp2Mid);
	get_eigen_values(rR, &vp1R, &vp2R);

	/* 1-Shock */
	if (rM < rL) {
		x1 = get_shock(1, rL, rM);
		vp1m = (x1 * rM - rL) / (x1 - 1);
		vp1p = vp1m;
		/* 1-Rarefaction */
	} else {
		vp1m = vp1L;
		vp1p = vp1Mid;
	}

	/* 2-Shock */
	if (rM > rR) {
		x2 = get_li(2, rM, rR);
		vp2m = (x2 * rR - rM) / (x2 - 1);
		vp2p = vp2m;
		/* 2-Rarefaction */
	} else {
		vp2m = vp2Mid;
		vp2p = vp2R;
	}

	/* Compute the solution u according to xi from left to right */
	if (xi < vp1m) {
		u[0] = uL[0];
		u[1] = uL[1];
	} else if (xi >= vp1m && xi <= vp1p) {
		r = get_r_brent(1, xi);
		u[0] = uL[0] * get_rarefaction(1, rL, r);
		u[1] = r * u[0];
	} else if (xi >= vp1p && xi < vp2m) {
		u[0] = rhoM;
		u[1] = rhoM * rM;
	} else if (xi >= vp2m && xi <= vp2p) {
		r = get_r_brent(2, xi);
		u[0] = uR[0] * get_rarefaction(2, rR, r);
		u[1] = r * u[0];
	} else if (xi > vp2p) {
		u[0] = uR[0];
		u[1] = uR[1];
	}
}
