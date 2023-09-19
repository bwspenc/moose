//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "Function.h"

/**
 * Defines level sets for a set of rods
 */
class LevelSetRods : public Function
{
public:
  static InputParameters validParams();

  LevelSetRods(const InputParameters & parameters);

  using Function::value;
  virtual Real value(Real /*t*/, const Point & p) const override;

protected:
  /// End coordinates of the rods
  const std::vector<RealVectorValue> & _line_endpoints;

  /// The radius of the rods
  const Real & _radius;

  /*
   * Get the shortest distance to project the point to the line
   * @param lp1 First point on the line
   * @param lp1 Second point on the line
   * @param p Point to project to the line
   * @return Projected distance
   */
  Real distFromLine(const RealVectorValue & lp1, const RealVectorValue & lp2, const RealVectorValue p) const;
};
