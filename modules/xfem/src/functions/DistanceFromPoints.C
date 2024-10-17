//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DistanceFromPoints.h"
#include <limits>

registerMooseObject("XFEMApp", DistanceFromPoints);

InputParameters
DistanceFromPoints::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Function defining the minimum distance from a set of points");
  params.addParam<std::vector<RealVectorValue>>(
      "point_coordinates", "Coordinates of the points from which to compute minimum distance");
  return params;
}

DistanceFromPoints::DistanceFromPoints(const InputParameters & parameters)
  : Function(parameters),
    _point_coordinates(getParam<std::vector<RealVectorValue>>("point_coordinates"))
{
}

Real
DistanceFromPoints::value(Real /*t*/, const Point & p) const
{
  Real min_dist = std::numeric_limits<Real>::max();
  for (unsigned int i = 0; i < _point_coordinates.size(); ++i)
    min_dist = std::min(distFromPoint(_point_coordinates[i], p), min_dist);
  return (min_dist);
}

Real
DistanceFromPoints::distFromPoint(const RealVectorValue & point, const RealVectorValue p) const
{
  return ((p-point).norm());
}
