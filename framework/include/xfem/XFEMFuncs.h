/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef XFEMFUNCS_H
#define XFEMFUNCS_H

#include <vector>
#include "libmesh/point.h"

using namespace libMesh;

void dunavant_rule2(const Real* wts, const Real* a, const Real* b, const unsigned int* permutation_ids,
                    unsigned int n_wts, std::vector<Point> &points, std::vector<Real> &weights);

void stdQuadr2D(unsigned int nen, unsigned int iord, std::vector<std::vector<Real> > &sg2);

void wissmannPoints(unsigned int nqp, std::vector<std::vector<Real> > &wss);

void shapeFunc2D(unsigned int nen, std::vector<Real> &ss, std::vector<Point> &xl,
                 std::vector<std::vector<Real> > &shp, Real &xsj, bool natl_flg);

double r8vec_norm(int n, double a[]);

void r8vec_copy(int n, double a1[], double a2[]);

bool r8vec_eq (int n, double a1[], double a2[]);

double r8vec_dot_product(int n, double a1[], double a2[]);

bool line_exp_is_degenerate_nd(int dim_num, double p1[], double p2[]);

int plane_normal_line_exp_int_3d(double pp[3], double normal[3], double p1[3], double p2[3], double pint[3]);

void normalizePoint(Point & p);

#endif
