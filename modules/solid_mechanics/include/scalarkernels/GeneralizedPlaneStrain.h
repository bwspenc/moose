//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GenericKernelScalar.h"
#include "RankTwoTensorForward.h"
#include "RankFourTensorForward.h"

class Function;

/**
 * Assembles the generalized plane strain equilibrium equation into the scalar out-of-plane strain
 * variable.
 */
template <bool is_ad>
class GeneralizedPlaneStrainTempl : public GenericKernelScalar<is_ad>
{
public:
  static InputParameters validParams();

  GeneralizedPlaneStrainTempl(const InputParameters & parameters);

  void computeResidualAndJacobian() override;

protected:
  GenericReal<is_ad> computeQpResidual() override { return 0; }
  GenericReal<is_ad> computeScalarQpResidual() override;
  Real computeScalarQpJacobian() override;

  /// Base name of the material system
  const std::string _base_name;

  /// The material Jacobian used by the non-AD specialization
  const MaterialProperty<RankFourTensor> * const _Jacobian_mult;

  /// The stress tensor
  const GenericMaterialProperty<RankTwoTensor, is_ad> & _stress;

  /// Function defining applied out-of-plane pressure
  const Function * const _out_of_plane_pressure_function;

  /// Material property defining applied out-of-plane pressure
  const GenericMaterialProperty<Real, is_ad> & _out_of_plane_pressure_material;

  /// Factor applied to out-of-plane pressure applied by function and material
  const Real _pressure_factor;

  /// The direction of the out-of-plane strain scalar variable
  unsigned int _scalar_out_of_plane_strain_direction;

  usingGenericKernelScalarMembers;
};

using GeneralizedPlaneStrain = GeneralizedPlaneStrainTempl<false>;
using ADGeneralizedPlaneStrain = GeneralizedPlaneStrainTempl<true>;
