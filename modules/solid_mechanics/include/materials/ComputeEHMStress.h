//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeStressBase.h"
#include "GuaranteeConsumer.h"

/**
 * ComputeEHMStress computes the stress following elasticity
 * theory for finite strains
 */
class ComputeEHMStress : public ComputeStressBase, public GuaranteeConsumer
{
public:
  static InputParameters validParams();

  ComputeEHMStress(const InputParameters & parameters);

  void initialSetup() override;
  virtual void initQpStatefulProperties() override;

protected:
  virtual void computeQpStress() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  /// Rotation up to current step "n" to compute anisotropic elasticity tensor
  MaterialProperty<RankTwoTensor> & _rotation_total;
  /// Rotation up to "n - 1" (previous) step to compute anisotropic elasticity tensor
  const MaterialProperty<RankTwoTensor> & _rotation_total_old;
  /// Strain increment material property
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  /// Rotation increment material property
  const MaterialProperty<RankTwoTensor> & _rotation_increment;
  /// Old state of the stress tensor material property
  const MaterialProperty<RankTwoTensor> & _stress_old;

  /**
   * The old elastic strain is used to calculate the old stress in the case
   * of variable elasticity tensors
   */
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;


  std::vector<RankFourTensor> _coefficient_data;

  /// Number of components
  const unsigned int _n_parts;
  /// Number of slip systems
  const unsigned int _n_slip_sys;
  /// Vector of initial slip resistances provided by the user
  const std::vector<Real> _initial_slip_resistance;
  /// Material property for current slip resistances
  MaterialProperty<std::vector<std::vector>> & _slip_resistance;
  /// Old state of material property for current slip resistances
  const MaterialProperty<std::vector<std::vector>> & _slip_resistance_old;
  /// VectorPostprocessor containing material constants
  const VectorPostprocessorValue & _matl_consts_vpp_value;
};
