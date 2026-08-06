//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GeneralizedPlaneStrain.h"

#include "Assembly.h"
#include "Function.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("SolidMechanicsApp", GeneralizedPlaneStrain);
registerMooseObject("SolidMechanicsApp", ADGeneralizedPlaneStrain);

template <bool is_ad>
InputParameters
GeneralizedPlaneStrainTempl<is_ad>::validParams()
{
  InputParameters params = GenericKernelScalar<is_ad>::validParams();
  params.addClassDescription(
      "Assembles the generalized plane strain residual into the scalar out-of-plane strain "
      "variable.");
  params.renameCoupledVar("scalar_variable",
                          "scalar_out_of_plane_strain",
                          "Scalar variable for generalized plane strain");
  params.set<bool>("compute_field_residuals") = false;
  params.addParam<FunctionName>("out_of_plane_pressure_function",
                                "Function used to prescribe pressure (applied toward the body) in "
                                "the out-of-plane direction");
  params.addDeprecatedParam<FunctionName>(
      "out_of_plane_pressure",
      "Function used to prescribe pressure (applied toward the body) in the out-of-plane direction "
      "(y for 1D Axisymmetric or z for 2D Cartesian problems)",
      "This has been replaced by 'out_of_plane_pressure_function'");
  params.addParam<MaterialPropertyName>("out_of_plane_pressure_material",
                                        "0",
                                        "Material used to prescribe pressure (applied toward the "
                                        "body) in the out-of-plane direction");
  MooseEnum out_of_plane_direction("x y z", "z");
  params.addParam<MooseEnum>("out_of_plane_direction",
                             out_of_plane_direction,
                             "The direction of the out-of-plane strain.");
  params.addDeprecatedParam<Real>(
      "factor",
      "Scale factor applied to prescribed out-of-plane pressure (both material and function)",
      "This has been replaced by 'pressure_factor'");
  params.addParam<Real>(
      "pressure_factor",
      "Scale factor applied to prescribed out-of-plane pressure (both material and function)");
  params.addParam<std::string>("base_name", "Material properties base name");

  return params;
}

template <bool is_ad>
GeneralizedPlaneStrainTempl<is_ad>::GeneralizedPlaneStrainTempl(const InputParameters & parameters)
  : GenericKernelScalar<is_ad>(parameters),
    _base_name(this->isParamValid("base_name")
                   ? this->template getParam<std::string>("base_name") + "_"
                   : ""),
    _Jacobian_mult(
        is_ad ? nullptr
              : &this->template getMaterialProperty<RankFourTensor>(_base_name + "Jacobian_mult")),
    _stress(this->template getGenericMaterialProperty<RankTwoTensor, is_ad>(_base_name + "stress")),
    _out_of_plane_pressure_function(parameters.isParamSetByUser("out_of_plane_pressure_function")
                                        ? &this->getFunction("out_of_plane_pressure_function")
                                    : parameters.isParamSetByUser("out_of_plane_pressure")
                                        ? &this->getFunction("out_of_plane_pressure")
                                        : nullptr),
    _out_of_plane_pressure_material(
        this->template getGenericMaterialProperty<Real, is_ad>("out_of_plane_pressure_material")),
    _pressure_factor(parameters.isParamSetByUser("pressure_factor")
                         ? this->template getParam<Real>("pressure_factor")
                     : parameters.isParamSetByUser("factor")
                         ? this->template getParam<Real>("factor")
                         : 1.0),
    _scalar_out_of_plane_strain_direction(
        this->template getParam<MooseEnum>("out_of_plane_direction"))
{
  if (parameters.isParamSetByUser("out_of_plane_pressure_function") &&
      parameters.isParamSetByUser("out_of_plane_pressure"))
    this->paramError("out_of_plane_pressure_function",
                     "Cannot specify both 'out_of_plane_pressure_function' and "
                     "'out_of_plane_pressure'.");
  if (parameters.isParamSetByUser("pressure_factor") && parameters.isParamSetByUser("factor"))
    this->paramError("pressure_factor", "Cannot specify both 'pressure_factor' and 'factor'.");
}

template <bool is_ad>
void
GeneralizedPlaneStrainTempl<is_ad>::computeResidualAndJacobian()
{
  if constexpr (is_ad)
    GenericKernelScalar<is_ad>::computeResidualAndJacobian();
  else
  {
    GenericKernelScalar<is_ad>::computeResidual();
    if (this->_is_implicit)
      GenericKernelScalar<is_ad>::computeJacobian();
  }
}

template <bool is_ad>
GenericReal<is_ad>
GeneralizedPlaneStrainTempl<is_ad>::computeScalarQpResidual()
{
  if (this->_assembly.coordSystem() == Moose::COORD_RZ)
    _scalar_out_of_plane_strain_direction = 1;
  else if (this->_assembly.coordSystem() != Moose::COORD_XYZ)
    mooseError("Unsupported coordinate system for generalized plane strain formulation");

  const Real function_pressure =
      _out_of_plane_pressure_function
          ? _out_of_plane_pressure_function->value(this->_t, _q_point[_qp])
          : 0.0;
  const auto pressure =
      (function_pressure + _out_of_plane_pressure_material[_qp]) * _pressure_factor;

  return _stress[_qp](_scalar_out_of_plane_strain_direction,
                      _scalar_out_of_plane_strain_direction) +
         pressure;
}

template <bool is_ad>
Real
GeneralizedPlaneStrainTempl<is_ad>::computeScalarQpJacobian()
{
  mooseAssert(!is_ad, "computeScalarQpJacobian should not be called by ADGeneralizedPlaneStrain");
  if (this->_assembly.coordSystem() == Moose::COORD_RZ)
    _scalar_out_of_plane_strain_direction = 1;
  else if (this->_assembly.coordSystem() != Moose::COORD_XYZ)
    mooseError("Unsupported coordinate system for generalized plane strain formulation");

  return (*_Jacobian_mult)[_qp](_scalar_out_of_plane_strain_direction,
                                _scalar_out_of_plane_strain_direction,
                                _scalar_out_of_plane_strain_direction,
                                _scalar_out_of_plane_strain_direction);
}

template class GeneralizedPlaneStrainTempl<false>;
template class GeneralizedPlaneStrainTempl<true>;
