//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeometricCut2DUserObject.h"

// Forward declarations

class LineSegmentCutSetUserObject : public GeometricCut2DUserObject
{
public:
  /**
   * Build parameters configuring the line-segment cut set user object.
   *
   * @return Input parameters describing the set of line segment cuts.
   */
  static InputParameters validParams();

  LineSegmentCutSetUserObject(const InputParameters & parameters);

  /**
   * Retrieve the ordered crack front points associated with the line segments.
   *
   * @param num_crack_front_points Number of crack front points expected by the caller.
   * @return Collection of crack front points sampled from the line segments.
   */
  virtual const std::vector<Point>
  getCrackFrontPoints(unsigned int num_crack_front_points) const override;

  /**
   * Return normals associated with the crack front points.
   *
   * @param num_crack_front_points Number of normals requested.
   * @return Crack front normals corresponding to each requested point.
   */
  virtual const std::vector<RealVectorValue>
  getCrackPlaneNormals(unsigned int num_crack_front_points) const override;

  /**
   * Get the cut location information stored for the line segments.
   *
   * @return Vector describing fractional cut locations within host elements.
   */
  virtual std::vector<Real> getCutData() const { return _cut_data; };

protected:
  std::vector<Real> _cut_data;
};
