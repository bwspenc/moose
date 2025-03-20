//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeEHMStress.h"

registerMooseObject("SolidMechanicsApp", ComputeEHMStress);

InputParameters
ComputeEHMStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addRequiredParam<std::vector<Real>>("initial_slip_resistance",
                                             "Vector of initial resistance of the components of the slip system");
  params.addRequiredParam<unsigned int>("num_parts", "Number of parts");
  params.addRequiredParam<unsigned int>("num_slip_systems", "Number of slip systems");
  params.addRequiredParam<VectorPostprocessorName>("material_constants", "Name of vectorpostprocessor that provides material constants");
  params.addClassDescription("Compute stress using EHM"); //TODO expand on this
  return params;
}

ComputeEHMStress::ComputeEHMStress(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
    GuaranteeConsumer(this),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _rotation_total(declareProperty<RankTwoTensor>(_base_name + "rotation_total")),
    _rotation_total_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "rotation_total")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _rotation_increment(
        getMaterialPropertyByName<RankTwoTensor>(_base_name + "rotation_increment")),
    _stress_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "stress")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain")),

    _n_parts(getParam<unsigned int>("num_parts")),
    _n_slip_sys(getParam<unsigned int>("num_slip_systems")),
    _initial_slip_reisistance(getParam<std:vector<Real>>("slip_resistance")),
    _slip_resistance(declareProperty<std::vector<std::vector>>(_base_name + "slip_resistance")),
    _slip_resistance_old(getMaterialPropertyOldByName<std::vector<std::vector>>(_base_name + "slip_resistance")),
    _slip_resistance_old(getMaterialPropertyOldByName<std::vector<std::vector>>(_base_name + "slip_resistance")),
    _matl_consts_vpp_value(getVectorPostprocessorValue("material_constants", "column_0")) //TODO maybe the column name needs to be changed
{

  //TODO error check to make sure _initial_slip_resistance is the right size

  //TODO load data file
  for (unsigned int i = 0; i< num_parts; ++i)
  {
    _coefficient_data[i]= some_rank_four_tensor;
    RankFourTensor temp_tensor;
    std::vector<Real> temp_vec;
    for (unsigned int j = 0; j< 36; ++j)
      temp_vec[j] = _matl_consts_vpp_value[j+i*36]; //TODO need to get this right
    temp_tensor.fillFromInputVector(temp_vec, RankFourTensor::symmetric_isotropic_E_nu); //TODO replace with correct fill method

    _coefficient_data[i] = temp_tensor;
  }


}

void
ComputeEHMStress::initialSetup()
{
}

void
ComputeEHMStress::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  RankTwoTensor identity_rotation(RankTwoTensor::initIdentity);

  _rotation_total[_qp] = identity_rotation;

  //TODO Note that _initial_slip_resistance is a vector of vectors, so need to reshape
  _slip_resistance[_qp].resize(_n_parts*12);
  _slip_resistance[_qp] = _initial_slip_resistance;
}

void
ComputeEHMStress::computeQpStress()
{
  // Calculate the stress in the intermediate configuration
  RankTwoTensor intermediate_stress;

  if (hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
  {
    //TODO: Initialize the current slip resistance. You could copy it in like this:
    _slip_resistance[_qp] = _slip_resistance_old[_qp];
    //or resize and initialize it:
    _slip_resistance[_qp].resize(_n_parts*12);


    intermediate_stress =
        _elasticity_tensor[_qp] * (_elastic_strain_old[_qp] + _strain_increment[_qp]);




    // Compute dstress_dstrain
    _Jacobian_mult[_qp] = _elasticity_tensor[_qp]; // This is NOT the exact jacobian
  }
  else
  {
    // Rotate elasticity tensor to the intermediate configuration
    // That is, elasticity tensor is defined in the previous time step
    // This is consistent with the definition of strain increment
    // The stress is projected onto the current configuration a few lines below
    RankFourTensor elasticity_tensor_rotated = _elasticity_tensor[_qp];
    elasticity_tensor_rotated.rotate(_rotation_total_old[_qp]);

    intermediate_stress =
        elasticity_tensor_rotated * (_elastic_strain_old[_qp] + _strain_increment[_qp]);

    // Update current total rotation matrix to be used in next step
    _rotation_total[_qp] = _rotation_increment[_qp] * _rotation_total_old[_qp];

    // Compute dstress_dstrain
    _Jacobian_mult[_qp] = elasticity_tensor_rotated; // This is NOT the exact jacobian
  }

  // Rotate the stress state to the current configuration
  _stress[_qp] =
      _rotation_increment[_qp] * intermediate_stress * _rotation_increment[_qp].transpose();

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = _mechanical_strain[_qp];
}
