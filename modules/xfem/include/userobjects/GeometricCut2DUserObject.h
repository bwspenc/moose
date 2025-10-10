//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeometricCutUserObject.h"

// Forward declarations

class GeometricCut2DUserObject : public GeometricCutUserObject
{
public:
  /**
   * Build parameters describing the two-dimensional geometric cut user object.
   *
   * @return Input parameters defining the 2D geometric cut configuration.
   */
  static InputParameters validParams();

  GeometricCut2DUserObject(const InputParameters & parameters);

  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutEdge> & cut_edges,
                                    std::vector<Xfem::CutNode> & cut_nodes) const override;
  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutFace> & cut_faces) const override;

  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_edges,
                                     std::vector<Xfem::CutEdge> & cut_edges) const override;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_faces,
                                     std::vector<Xfem::CutFace> & cut_faces) const override;

protected:
  std::vector<std::pair<Point, Point>> _cut_line_endpoints;

  /**
   * Find the fractional distance along a specified cut line at the current time.
   *
   * This is used for time-based propagation along a line segment.
   *
   * @param cut_num Index of the cut being queried.
   * @return Current fractional distance along the cut line.
   */
  virtual Real cutFraction(unsigned int cut_num) const;

  /// Vector of start/end times for each cut segment
  std::vector<std::pair<Real, Real>> _cut_time_ranges;
};
