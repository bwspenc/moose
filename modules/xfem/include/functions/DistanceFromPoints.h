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
class DistanceFromPoints : public Function
{
public:
  static InputParameters validParams();

  DistanceFromPoints(const InputParameters & parameters);

  using Function::value;
  virtual Real value(Real /*t*/, const Point & p) const override;

protected:
  /// Point coordinates
  const std::vector<RealVectorValue> & _point_coordinates;

  /*
   * Get the distance between two points
   * @param point Point to be checked against
   * @param p Point to project to the line
   * @return Distance
   */
  Real distFromPoint(const RealVectorValue & point, const RealVectorValue p) const;
};
