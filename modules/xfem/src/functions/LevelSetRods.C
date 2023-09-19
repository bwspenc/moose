//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetRods.h"
#include <limits>

registerMooseObject("XFEMApp", LevelSetRods);

InputParameters
LevelSetRods::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Level set function defining a set of rods");
  params.addParam<std::vector<RealVectorValue>>(
      "line_endpoints", "Endpoints of the centerlines of the rods");
  params.addParam<Real>("radius", "The radius of the rods.");
  return params;
}

LevelSetRods::LevelSetRods(const InputParameters & parameters)
  : Function(parameters),
    _line_endpoints(getParam<std::vector<RealVectorValue>>("line_endpoints")),
    _radius(getParam<Real>("radius"))
{
  if (_line_endpoints.size() % 2 != 0)
    paramError("line_endpoints", "'line_endpoints' must be prescribed as sets of 2 groups of 3 real values");
}

Real
LevelSetRods::value(Real /*t*/, const Point & p) const
{
  Real min_dist = std::numeric_limits<Real>::max();
  for (unsigned int i = 0; i < _line_endpoints.size() / 2; ++i)
    min_dist = std::min(distFromLine(_line_endpoints[i * 2], _line_endpoints[i * 2 + 1], p), min_dist);
  return (min_dist - _radius);
}

Real
LevelSetRods::distFromLine(const RealVectorValue & lp1, const RealVectorValue & lp2, const RealVectorValue p) const
{
  const RealVectorValue line_vec = (lp2 - lp1).unit();
  const Real proj_point_dist = line_vec * (p - lp1);
  const RealVectorValue proj_point = lp1 + proj_point_dist * line_vec;
  const RealVectorValue line_to_point = p - proj_point;
  return (line_to_point.norm());
}
