//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

// Forward Declarations

class SurfaceTractions : public Material
{
public:
  static InputParameters validParams();

  SurfaceTractions(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Base name of the material system
  const std::string _base_name;

  /// Material property storing the stresses
  const MaterialProperty<RankTwoTensor> & _stress;
  /// Material property storing the normals
  MaterialProperty<RealVectorValue> & _normal;
  /// Material property storing the tractions (in Cartesian coordinates)
  MaterialProperty<RealVectorValue> & _tractions;
  /// Material property storing the normal traction
  MaterialProperty<Real> & _normal_traction;
  /// Material property storing the max shear traction
  MaterialProperty<Real> & _max_shear_traction;
  /// Material property storing the min shear traction
  MaterialProperty<Real> & _min_shear_traction;

  ///Normal vectors at the quadrature points
  const MooseArray<Point> & _normals;
};
