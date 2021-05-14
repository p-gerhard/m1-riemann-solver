/*
 * Copyright (c) 2017-2020 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * M1-riemann-solver is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "m1_riemann.h"

double (*get_chi)(const double);

static double m1_s2_get_chi_interp(const double r)
{
	double r_abs = fabs(r);

	if (r_abs > 1 - 1e-4)
		return 1 - 1e-4;

	const double p[7] = { 0.325269105027,	 -1.216531569434, 2.116995591837,
												-2.514343153655, 2.271654209692,	-1.312060956603,
												0.333333397667 };

	const double q[7] = { -0.327086152002, 1.253825802327, -0.780269389781,
												-2.820874826259, 5.614892441680, -3.936171251864,
												1.000000000000 };

	double num = p[0];
	double den = q[0];

	for (unsigned int i = 1; i < 7; i++) {
		num = num * r_abs + p[i];
		den = den * r_abs + q[i];
	}

	double chi = fmin(num / den, 1. - 1e-4);

	return chi;
}

static int dump_solution(const char filename[1024], const double *x,
												 const double *u, const int N)
{
	FILE *f = fopen(filename, "w+");

	if (f != NULL) {
		for (int i = 0; i < N; i++) {
			fprintf(f, "%.12f %.12f %.12f\n", x[i], u[2 * i + 0], u[2 * i + 1]);
		}
		fclose(f);
		return 0;
	} else {
		return -1;
	}
}

int main(int argc, char *argv[])
{
	/* Riemann problem parameters */
	const double xmin = -1;
	const double xmax = 1;
	const double tmax = 1;
	const int N = 500;
	const double uL[2] = { 0.7, 0.69 };
	const double uR[2] = { 0.8, -0.1 };

	/* Select the propper Eddington factor function */
	get_chi = &m1_s2_get_chi_interp;

	/* Generate the mesh */
	double *mesh = (double *)malloc(N * sizeof(double));
	const double dx = fabs(xmax - xmin) / (N - 1);

	for (int i = 0; i < N; i++) {
		mesh[i] = xmin + i * dx;
	}

	/* Solution buffers */
	double u[2] = { 0, 0 };
	double *u_ex = (double *)malloc(2 * N * sizeof(double));

	clock_t start = clock();
	for (int k = 0; k < N; k++) {
		solve_riemann(mesh[k] / tmax, uL, uR, u);
		// printf("u(%f,%f)\n ", u[0], u[1]);
		u_ex[2 * k + 0] = u[0];
		u_ex[2 * k + 1] = u[1];
	}
	
	clock_t stop = clock();
	double elapsed = (double)(stop - start) * 1000.0 / CLOCKS_PER_SEC;	
	printf("Elapsed time %f (s)\n", elapsed / 1000);
	
	/* Dump solution */
	const char filename[1024] = "u_s2_riemann_exact.out";
	dump_solution(filename, mesh, u_ex, N);

	free(mesh);
	free(u_ex);
	// SolutionTrace(1);
	return EXIT_SUCCESS;
}
