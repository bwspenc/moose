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

#ifndef XFEM_CIRCLE_CUT_H
#define XFEM_CIRCLE_CUT_H

#include "XFEMGeometricCut3D.h"

class XFEMCircleCut : public XFEMGeometricCut3D
{
public:

  XFEMCircleCut(std::vector<Real> square_nodes);
  ~XFEMCircleCut();

private:

  std::vector<Point> _vertices;
  Point _center;
  Point _normal;
  Real _radius;
  Real _angle;

  virtual bool intersectWithEdge(Point p1, Point p2, Point &pint);

  virtual bool isInsideCutPlane(Point p);
};

#endif
