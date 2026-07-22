//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "RankTwoTensorForward.h"
#include "RankFourTensorForward.h"

class Function;
class MooseVariableScalar;

class GeneralizedPlaneStrain : public Kernel
{
public:
  static InputParameters validParams();

  GeneralizedPlaneStrain(const InputParameters & parameters);

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeResidualAndJacobian() override;

protected:
  virtual Real computeQpResidual() override { return 0; }

  /// Compute the scalar out-of-plane strain residual contribution
  void computeScalarResidual();

  /// Compute the scalar out-of-plane strain diagonal Jacobian contribution
  void computeScalarJacobian();

  /// Base name of the material system
  const std::string _base_name;

  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  /// The stress tensor
  const MaterialProperty<RankTwoTensor> & _stress;

  /// Function defining applied out-of-plane pressure
  const Function * _out_of_plane_pressure_function;

  /// Material property defining applied out-of-plane pressure
  const MaterialProperty<Real> & _out_of_plane_pressure_material;

  /// Factor applied to out-of-plane pressure applied by function and material
  const Real _pressure_factor;

  /// Scalar out-of-plane strain variable
  const MooseVariableScalar & _scalar_out_of_plane_strain_var;

  /// The direction of the out-of-plane strain scalar variable
  unsigned int _scalar_out_of_plane_strain_direction;
};
