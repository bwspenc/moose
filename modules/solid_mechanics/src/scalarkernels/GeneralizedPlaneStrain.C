//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GeneralizedPlaneStrain.h"

// MOOSE includes
#include "Assembly.h"
#include "Function.h"
#include "InputParameters.h"
#include "MooseVariableScalar.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

#include "libmesh/quadrature.h"

#include <cmath>

registerMooseObject("SolidMechanicsApp", GeneralizedPlaneStrain);

InputParameters
GeneralizedPlaneStrain::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Generalized Plane Strain scalar residual Kernel");
  params.addRequiredCoupledVar("scalar_out_of_plane_strain",
                               "Scalar variable for generalized plane strain");
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
  MooseEnum outOfPlaneDirection("x y z", "z");
  params.addParam<MooseEnum>(
      "out_of_plane_direction", outOfPlaneDirection, "The direction of the out-of-plane strain.");
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

GeneralizedPlaneStrain::GeneralizedPlaneStrain(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _Jacobian_mult(getMaterialProperty<RankFourTensor>(_base_name + "Jacobian_mult")),
    _stress(getMaterialProperty<RankTwoTensor>(_base_name + "stress")),
    _out_of_plane_pressure_function(parameters.isParamSetByUser("out_of_plane_pressure_function")
                                        ? &getFunction("out_of_plane_pressure_function")
                                    : parameters.isParamSetByUser("out_of_plane_pressure")
                                        ? &getFunction("out_of_plane_pressure")
                                        : nullptr),
    _out_of_plane_pressure_material(getMaterialProperty<Real>("out_of_plane_pressure_material")),
    _pressure_factor(parameters.isParamSetByUser("pressure_factor")
                         ? getParam<Real>("pressure_factor")
                     : parameters.isParamSetByUser("factor") ? getParam<Real>("factor")
                                                             : 1.0),
    _scalar_out_of_plane_strain_var(*getScalarVar("scalar_out_of_plane_strain", 0)),
    _scalar_out_of_plane_strain_direction(getParam<MooseEnum>("out_of_plane_direction"))
{
  if (parameters.isParamSetByUser("out_of_plane_pressure_function") &&
      parameters.isParamSetByUser("out_of_plane_pressure"))
    mooseError("Cannot specify both 'out_of_plane_pressure_function' and 'out_of_plane_pressure'");
  if (parameters.isParamSetByUser("pressure_factor") && parameters.isParamSetByUser("factor"))
    mooseError("Cannot specify both 'pressure_factor' and 'factor'");
}

void
GeneralizedPlaneStrain::computeResidual()
{
  Kernel::computeResidual();
  computeScalarResidual();
}

void
GeneralizedPlaneStrain::computeJacobian()
{
  Kernel::computeJacobian();
  computeScalarJacobian();
}

void
GeneralizedPlaneStrain::computeResidualAndJacobian()
{
  computeResidual();

  if (_is_implicit)
  {
    prepareShapes(_var.number());
    Kernel::computeJacobian();
    computeScalarJacobian();
  }
}

void
GeneralizedPlaneStrain::computeScalarResidual()
{
  if (_assembly.coordSystem() == Moose::COORD_XYZ)
    _scalar_out_of_plane_strain_direction = getParam<MooseEnum>("out_of_plane_direction");
  else if (_assembly.coordSystem() == Moose::COORD_RZ)
    _scalar_out_of_plane_strain_direction = 1;
  else
    mooseError("Unsupported coordinate system for generalized plane strain formulation");

  std::vector<Real> residual(_scalar_out_of_plane_strain_var.dofIndices().size(), 0.0);
  std::vector<Real> reference_residual(_scalar_out_of_plane_strain_var.dofIndices().size(), 0.0);

  for (const auto qp : make_range(_qrule->n_points()))
  {
    _qp = qp;
    const Real out_of_plane_pressure =
        ((_out_of_plane_pressure_function
              ? _out_of_plane_pressure_function->value(_t, _q_point[_qp])
              : 0.0) +
         _out_of_plane_pressure_material[_qp]) *
        _pressure_factor;

    const Real stress =
        _stress[_qp](_scalar_out_of_plane_strain_direction, _scalar_out_of_plane_strain_direction);

    for (const auto i : index_range(residual))
    {
      _i = i;
      residual[_i] += _JxW[_qp] * _coord[_qp] * (stress + out_of_plane_pressure);
      reference_residual[_i] += std::abs(_JxW[_qp] * _coord[_qp] * stress);
    }
  }

  prepareVectorTag(_assembly,
                   _scalar_out_of_plane_strain_var.number(),
                   TaggingInterface::ResidualTagType::NonReference);
  addResiduals(_assembly,
               residual,
               _scalar_out_of_plane_strain_var.dofIndices(),
               _scalar_out_of_plane_strain_var.scalingFactor());

  prepareVectorTag(_assembly,
                   _scalar_out_of_plane_strain_var.number(),
                   TaggingInterface::ResidualTagType::Reference);
  addResiduals(_assembly,
               reference_residual,
               _scalar_out_of_plane_strain_var.dofIndices(),
               _scalar_out_of_plane_strain_var.scalingFactor());
}

void
GeneralizedPlaneStrain::computeScalarJacobian()
{
  if (_assembly.coordSystem() == Moose::COORD_XYZ)
    _scalar_out_of_plane_strain_direction = getParam<MooseEnum>("out_of_plane_direction");
  else if (_assembly.coordSystem() == Moose::COORD_RZ)
    _scalar_out_of_plane_strain_direction = 1;
  else
    mooseError("Unsupported coordinate system for generalized plane strain formulation");

  DenseMatrix<Real> jacobian(_scalar_out_of_plane_strain_var.dofIndices().size(),
                             _scalar_out_of_plane_strain_var.dofIndices().size());

  for (const auto qp : make_range(_qrule->n_points()))
  {
    _qp = qp;
    for (const auto i : make_range(jacobian.m()))
      jacobian(i, i) += _JxW[_qp] * _coord[_qp] *
                        _Jacobian_mult[_qp](_scalar_out_of_plane_strain_direction,
                                            _scalar_out_of_plane_strain_direction,
                                            _scalar_out_of_plane_strain_direction,
                                            _scalar_out_of_plane_strain_direction);
  }

  prepareMatrixTag(_assembly,
                   _scalar_out_of_plane_strain_var.number(),
                   _scalar_out_of_plane_strain_var.number());
  addJacobian(_assembly,
              jacobian,
              _scalar_out_of_plane_strain_var.dofIndices(),
              _scalar_out_of_plane_strain_var.dofIndices(),
              _scalar_out_of_plane_strain_var.scalingFactor());
}
