/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef RADIALDISPLACEMENTAUX_H
#define RADIALDISPLACEMENTAUX_H

#include "AuxKernel.h"

class RadialDisplacementAux;

template<>
InputParameters validParams<RadialDisplacementAux>();

class RadialDisplacementAux : public AuxKernel
{
public:

  /**
   * Calcualtes the radial displacement for cylindrical or spherical
   * geometries. Works for 3D, 2D axisymmetric, 2D planar, and 1D geometries
   */
  RadialDisplacementAux(const InputParameters & parameters);

  virtual ~RadialDisplacementAux() {}

protected:
  /// Compute the value of the radial displacement
  virtual Real computeValue();

  /// Type of coordinate system
  Moose::CoordinateSystemType _coord_system;

  /// Number of displacment components.
  unsigned int _ndisp;
  /// Coupled variable values of the displacement components.
  std::vector<const VariableValue *> _disp_vals;

  /// First point used to define an axis
  RealVectorValue _axis_point_1;
  /// Second point used to define an axis
  RealVectorValue _axis_point_2;

  /// Point used to define an origin for a cylindrical (2D Cartesian) or spherical
  /// (3D Cartesian) system.
  RealVectorValue _origin;

  /// Has the axis been defined by the user?
  bool _have_axis;
  /// Has the origin been defined by the user?
  bool _have_origin;
};

#endif //RADIALDISPLACEMENTAUX_H
