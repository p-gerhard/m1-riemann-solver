/*
 * Copyright (c) 2017-2020 Pierre Gerhard <pierre.gerhard@gmail.com>
 *
 * M1-riemann-solver is free software; you can redistribute it and/or modify
 * it under the terms of the GPLv3 license. See LICENSE for details.
 */
#ifndef M1_RIEMANN
#define M1_RIEMANN

/* Used to switch the function computing the Eddington factor Chi(r)*/
extern double (*get_chi)(const double);

void solve_riemann(const double xi, const double uL[2], const double uR[2],
									 double u[2]);

#endif