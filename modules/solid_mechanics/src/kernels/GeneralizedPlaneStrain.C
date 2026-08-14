//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GeneralizedPlaneStrain.h"

#include "Function.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "MooseVariable.h"
#include "MooseVariableScalar.h"
#include "Assembly.h"
#include "SubblockIndexProvider.h"
#include "UserObject.h"

#include "libmesh/quadrature.h"

registerMooseObject("SolidMechanicsApp", GeneralizedPlaneStrain);
registerMooseObject("SolidMechanicsApp", ADGeneralizedPlaneStrain);

template <bool is_ad>
InputParameters
GeneralizedPlaneStrainTempl<is_ad>::validParams()
{
  InputParameters params = GenericKernelScalar<is_ad>::validParams();
  params.addClassDescription("Assembles the contributions of the area integral performed over "
                             "the elements to the residual for the out-of-plane strain scalar "
                             "variable, for both the automatic differentiation (AD) and non-AD "
                             "cases.");
  params.renameCoupledVar("scalar_variable",
                          "scalar_out_of_plane_strain",
                          "Scalar variable for generalized plane strain");
  params.makeParamRequired<std::vector<VariableName>>("scalar_out_of_plane_strain");
  params.addParam<FunctionName>("out_of_plane_pressure_function",
                                "Function used to prescribe pressure (applied toward the body) in "
                                "the out-of-plane direction");
  params.addParam<MaterialPropertyName>("out_of_plane_pressure_material",
                                        "0",
                                        "Material used to prescribe pressure (applied toward the "
                                        "body) in the out-of-plane direction");
  MooseEnum out_of_plane_direction("x y z", "z");
  params.addParam<MooseEnum>(
      "out_of_plane_direction", out_of_plane_direction, "The direction of the out-of-plane strain");
  params.addParam<Real>(
      "pressure_factor",
      1.0,
      "Scale factor applied to prescribed out-of-plane pressure (both material and function)");
  params.addParam<std::string>("base_name", "Material property base name");

  // Non-AD only: the displacement and temperature variables, and the eigenstrain derivatives
  // needed to fill the off-diagonal Jacobian coupling with those variables (AD instead captures
  // this coupling automatically through the stress material property's AD chain). Declared
  // unconditionally so GeneralizedPlaneStrainAction can forward its own params to either kernel
  // type via applyParameters() without special-casing is_ad.
  params.addParam<std::vector<VariableName>>("displacements", "The displacement variables");
  params.addCoupledVar("temperature", "The temperature variable");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", {}, "List of eigenstrains to be applied in this strain calculation");

  // Multi-region support: selects which of a set of per-region scalar_out_of_plane_strain
  // variables (and which elements) this instance is responsible for.
  params.addParam<UserObjectName>("subblock_index_provider",
                                  "SubblockIndexProvider user object name");
  params.addParam<unsigned int>("scalar_out_of_plane_strain_index",
                                0,
                                "The index number of scalar_out_of_plane_strain this kernel acts "
                                "on");

  // Non-AD only: whether this instance owns the out-of-plane strain scalar variable's own
  // residual and diagonal Jacobian contributions. Tautological for AD, which only ever has one
  // instance.
  params.addParam<bool>("primary",
                        true,
                        "Whether this instance owns the out-of-plane strain scalar variable's own "
                        "residual and diagonal Jacobian contributions");

  if constexpr (is_ad)
  {
    // The field row (this instance's attached displacement/temperature variable) is always zero
    // for AD: AD instead gets its d(field)/d(kappa) Jacobian entries for free through the
    // material property AD chain of the ordinary AD stress-divergence/heat-conduction kernels.
    params.set<bool>("compute_field_residuals") = false;
    params.suppressParameter<bool>("compute_field_residuals");
  }
  else
    params.addParam<bool>(
        "reference_residual_excludes_pressure",
        false,
        "Whether to exclude the applied out-of-plane pressure from the 'Reference' residual tag "
        "contribution, matching the behavior of the removed non-AD GeneralizedPlaneStrain "
        "ScalarKernel");

  return params;
}

template <bool is_ad>
GeneralizedPlaneStrainTempl<is_ad>::GeneralizedPlaneStrainTempl(const InputParameters & parameters)
  : DerivativeMaterialInterface<GenericKernelScalar<is_ad>>(parameters),
    _base_name(this->isParamValid("base_name") ? this->template getParam<std::string>("base_name") + "_"
                                               : ""),
    _stress(this->template getGenericMaterialProperty<RankTwoTensor, is_ad>(_base_name + "stress")),
    _out_of_plane_pressure_function(parameters.isParamSetByUser("out_of_plane_pressure_function")
                                        ? &this->getFunction("out_of_plane_pressure_function")
                                        : nullptr),
    _out_of_plane_pressure_material(
        this->template getMaterialProperty<Real>("out_of_plane_pressure_material")),
    _pressure_factor(this->template getParam<Real>("pressure_factor")),
    _out_of_plane_direction(this->template getParam<MooseEnum>("out_of_plane_direction")),
    _primary(this->template getParam<bool>("primary")),
    _subblock_id_provider(this->isParamValid("subblock_index_provider")
                              ? &this->template getUserObject<SubblockIndexProvider>(
                                    "subblock_index_provider")
                              : nullptr),
    _scalar_var_id(this->template getParam<unsigned int>("scalar_out_of_plane_strain_index")),
    _reference_residual_excludes_pressure(
        is_ad ? false : this->template getParam<bool>("reference_residual_excludes_pressure")),
    _Jacobian_mult(is_ad ? nullptr
                        : &this->template getMaterialProperty<RankFourTensor>(_base_name +
                                                                              "Jacobian_mult")),
    _temp_var(this->isCoupled("temperature") ? this->getVar("temperature", 0) : nullptr),
    _in_plane_component(-1)
{
  // The remaining members are only meaningful for the non-AD case: AD needs none of this because
  // it captures the coupling between the scalar variable and the displacement/temperature field
  // variables automatically through the stress material property's AD chain, whereas the non-AD
  // kernel must fill that off-diagonal Jacobian coupling explicitly, one field variable at a time.
  if constexpr (!is_ad)
  {
    if (this->isParamValid("displacements"))
    {
      const auto & disp_names = this->template getParam<std::vector<VariableName>>("displacements");
      for (const auto & disp_name : disp_names)
        _disp_var.push_back(&this->_subproblem.getStandardVariable(this->_tid, disp_name));

      if (_disp_var.size() >= 1 && _var.number() == _disp_var[0]->number())
        _in_plane_component = 0;
      else if (_disp_var.size() >= 2 && _var.number() == _disp_var[1]->number())
        _in_plane_component = 1;
    }

    const auto & eigenstrain_names =
        this->template getParam<std::vector<MaterialPropertyName>>("eigenstrain_names");
    if (_temp_var && _var.number() == _temp_var->number())
      for (const auto & eigenstrain_name : eigenstrain_names)
        _deigenstrain_dT.push_back(&this->template getMaterialPropertyDerivative<RankTwoTensor>(
            _base_name + eigenstrain_name, _temp_var->name()));
  }
}

template <bool is_ad>
void
GeneralizedPlaneStrainTempl<is_ad>::initialSetup()
{
  if (this->getBlockCoordSystem() == Moose::COORD_RZ)
    _out_of_plane_direction = 1;
  else if (this->getBlockCoordSystem() != Moose::COORD_XYZ)
    this->paramError("out_of_plane_direction",
                     "Generalized plane strain supports only Cartesian and axisymmetric "
                     "coordinate systems");
}

template <bool is_ad>
bool
GeneralizedPlaneStrainTempl<is_ad>::onOwnSubblock() const
{
  return !_subblock_id_provider ||
        _subblock_id_provider->getSubblockIndex(*_current_elem) == _scalar_var_id;
}

template <bool is_ad>
GenericReal<is_ad>
GeneralizedPlaneStrainTempl<is_ad>::computeQpResidual()
{
  return 0;
}

template <bool is_ad>
GenericReal<is_ad>
GeneralizedPlaneStrainTempl<is_ad>::computeScalarQpResidual()
{
  if (!_primary || !onOwnSubblock())
    return 0;

  const Real out_of_plane_pressure =
      ((_out_of_plane_pressure_function ? _out_of_plane_pressure_function->value(this->_t,
                                                                                 _q_point[_qp])
                                        : 0.0) +
       _out_of_plane_pressure_material[_qp]) *
      _pressure_factor;

  return _stress[_qp](_out_of_plane_direction, _out_of_plane_direction) + out_of_plane_pressure;
}

template <bool is_ad>
Real
GeneralizedPlaneStrainTempl<is_ad>::computeScalarQpJacobian()
{
  // This function will never be called for the AD version: AD differentiates
  // computeScalarQpResidual() directly instead. But because C++ does not support an optional
  // function declaration based on a template parameter, we must keep this template for all cases.
  mooseAssert(!is_ad,
             "In ADGeneralizedPlaneStrain, computeScalarQpJacobian should not be called. Check "
             "computeJacobian implementation.");

  if (!_primary || !onOwnSubblock())
    return 0;

  return (*_Jacobian_mult)[_qp](
      _out_of_plane_direction, _out_of_plane_direction, _out_of_plane_direction, _out_of_plane_direction);
}

template <bool is_ad>
Real
GeneralizedPlaneStrainTempl<is_ad>::computeQpOffDiagJacobianScalar(unsigned int svar_num)
{
  // Non-AD only; see computeScalarQpJacobian(). Guarded with if constexpr (rather than just
  // mooseAssert) because _grad_test/_grad_phi are AD-valued for is_ad, so the arithmetic below
  // does not type-check as returning a plain Real in that instantiation.
  if constexpr (is_ad)
    mooseError("In ADGeneralizedPlaneStrain, computeQpOffDiagJacobianScalar should not be called. "
               "Check computeOffDiagJacobianScalar implementation.");
  else
  {
    // Fills d-(this instance's field variable)-residual / d-kappa: only meaningful for kappa's own
    // column, and only on the subblock this instance's scalar variable is responsible for.
    if (svar_num != _kappa_var || !onOwnSubblock())
      return 0;

    if (_in_plane_component >= 0)
    {
      // Shift the in-plane component index (0 or 1, indexing into the tracked displacements) to
      // the actual global direction index, skipping over the out-of-plane direction. Mirrors
      // the removed GeneralizedPlaneStrainOffDiag::computeDispOffDiagJacobianScalar.
      unsigned int component = _in_plane_component;
      if (_out_of_plane_direction == 0)
        component += 1;
      else if (_out_of_plane_direction == 1 && component == 1)
        component += 1;

      return (*_Jacobian_mult)[_qp](
                 _out_of_plane_direction, _out_of_plane_direction, component, component) *
            _grad_test[_i][_qp](component);
    }
    else if (_temp_var && _var.number() == _temp_var->number())
    {
      Real factor = 0;
      for (const auto ies : index_range(_deigenstrain_dT))
        factor += ((*_Jacobian_mult)[_qp] * (*_deigenstrain_dT[ies])[_qp])(_out_of_plane_direction,
                                                                          _out_of_plane_direction);

      return factor * _test[_i][_qp];
    }

    return 0;
  }
}

template <bool is_ad>
Real
GeneralizedPlaneStrainTempl<is_ad>::computeScalarQpOffDiagJacobian(unsigned int jvar_num)
{
  // Non-AD only; see computeScalarQpJacobian(). Guarded with if constexpr; see the comment in
  // computeQpOffDiagJacobianScalar() above.
  if constexpr (is_ad)
    mooseError("In ADGeneralizedPlaneStrain, computeScalarQpOffDiagJacobian should not be called. "
               "Check computeOffDiagJacobian implementation.");
  else
  {
    // Fills d-kappa-residual / d-(this instance's field variable): only this instance's own
    // shape functions are available here, so this can only ever be nonzero for its own variable.
    if (jvar_num != _var.number() || !onOwnSubblock())
      return 0;

    if (_in_plane_component >= 0)
    {
      unsigned int component = _in_plane_component;
      if (_out_of_plane_direction == 0)
        component += 1;
      else if (_out_of_plane_direction == 1 && component == 1)
        component += 1;

      return (*_Jacobian_mult)[_qp](
                 _out_of_plane_direction, _out_of_plane_direction, component, component) *
            _grad_phi[_j][_qp](component);
    }
    else if (_temp_var && _var.number() == _temp_var->number())
    {
      Real factor = 0;
      for (const auto ies : index_range(_deigenstrain_dT))
        factor += ((*_Jacobian_mult)[_qp] * (*_deigenstrain_dT[ies])[_qp])(_out_of_plane_direction,
                                                                          _out_of_plane_direction);

      return factor * _phi[_j][_qp];
    }

    return 0;
  }
}

template <bool is_ad>
void
GeneralizedPlaneStrainTempl<is_ad>::computeScalarResidual()
{
  if constexpr (is_ad)
    mooseError("ADGeneralizedPlaneStrain fuses residual and Jacobian assembly via automatic "
              "differentiation; computeScalarResidual() should never be called.");
  else
  {
    if (!_reference_residual_excludes_pressure)
    {
      // Common case: reference residual is just abs() of the full residual, matching every other
      // MOOSE Kernel (and matching AD, which has no split-residual concept at all). This is a
      // behavior change from the removed non-AD GeneralizedPlaneStrain ScalarKernel, which
      // excluded the applied pressure from the reference residual by default; see the
      // 'reference_residual_excludes_pressure' backward-compatibility param below for the old
      // behavior.
      GenericKernelScalar<is_ad>::computeScalarResidual();
      return;
    }

    // kappa's row is only ever owned by the primary instance; the generic default above already
    // no-ops for non-primary instances via computeScalarQpResidual()'s own _primary guard, but
    // the split below re-reads _stress[_qp] directly rather than going through that already-
    // gated hook, so the guard must be explicit here too.
    if (!_primary)
      return;

    std::vector<Real> residual(_k_order);
    std::vector<Real> reference_residual(_k_order);
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    {
      for (_h = 0; _h < _k_order; _h++)
      {
        residual[_h] += _JxW[_qp] * _coord[_qp] * computeScalarQpResidual();
        if (onOwnSubblock())
          reference_residual[_h] +=
              std::abs(_JxW[_qp] * _coord[_qp] *
                      _stress[_qp](_out_of_plane_direction, _out_of_plane_direction));
      }
    }

    this->prepareVectorTag(
        _assembly, _kappa_var_ptr->number(), TaggingInterface::ResidualTagType::NonReference);
    for (const auto h : make_range(_k_order))
      _local_re(h) += residual[h];
    this->accumulateTaggedLocalResidual();

    this->prepareVectorTag(
        _assembly, _kappa_var_ptr->number(), TaggingInterface::ResidualTagType::Reference);
    for (const auto h : make_range(_k_order))
      _local_re(h) += reference_residual[h];
    this->accumulateTaggedLocalResidual();
  }
}

template class GeneralizedPlaneStrainTempl<false>;
template class GeneralizedPlaneStrainTempl<true>;
