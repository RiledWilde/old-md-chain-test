#pragma once

#ifndef GLOB_H
#define GLOB_H

/* Energy of the system */
typedef struct {
  double kin;
  double FENE;
  double harmonic;
} energy_struct;
energy_struct energy;

/* FENE potential parameters */
typedef struct {
  double k;
  double r;
  double r2;
} FENE;
FENE fene;

/* Particle data */
typedef struct {
  double r[3];
  double v[3];
  double f[3];
  int nb;
  int *bl;
} Particle;
Particle *part;

/* Polymer data */
typedef struct {
  int N;
  int MPC;
} polymer_struct;
polymer_struct poly;

/* Molecule data */
typedef struct {
  double cm[3];
} polymer_whole_struct;
polymer_whole_struct *mole;

/* Umbrella sampling */
typedef struct {
  unsigned int size;
  int *binn;
  double *xVal;
  double range[2];
  double k;
  double K;
  double d_eq;

  int count;
  double xiBar;
  double xi2Bar;
} UmPotential_struct;
UmPotential_struct um;

/*----------------------------------------------------------------*/
/* Most close integer of double x, give a double number back. */
double dround(double x) { return floor(x + 0.5); }

/*----------------------------------------------------------------*/
/* Compute the square of double x, give double number back. */
double SQR(double x) { return x*x; }

/*----------------------------------------------------------------*/
/* Compute the cube of double x, give double number back. */
double CUBE(double x) { return x*x*x; }

#endif