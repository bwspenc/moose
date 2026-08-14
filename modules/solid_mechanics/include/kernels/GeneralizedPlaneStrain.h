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
#include "DerivativeMaterialInterface.h"
#include "SubblockIndexProvider.h"
#include "ADRankTwoTensorForward.h"
#include "RankFourTensorForward.h"

class Function;

/**
 * Assembles the contributions of the area integral performed over the elements to the residual
 * for the out-of-plane strain scalar variable, for both the automatic differentiation (AD) and
 * non-AD cases.
 *
 * A single AD instance (attached to one in-plane displacement) is sufficient because AD captures
 * the coupling between the out-of-plane strain scalar variable and the displacement/temperature
 * field variables automatically through the stress material property's AD chain. The non-AD case
 * has no such automatic differentiation, so the Action instead sets up one instance per in-plane
 * displacement (plus one for temperature, if coupled), each responsible for filling the
 * off-diagonal Jacobian coupling between the out-of-plane strain scalar variable and its own
 * attached field variable. Exactly one non-AD instance ("primary") additionally owns the scalar
 * variable's own residual and diagonal Jacobian contributions.
 */
template <bool is_ad>
class GeneralizedPlaneStrainTempl : public DerivativeMaterialInterface<GenericKernelScalar<is_ad>>
{
public:
  static InputParameters validParams();

  GeneralizedPlaneStrainTempl(const InputParameters & parameters);

protected:
  void initialSetup() override;
  GenericReal<is_ad> computeQpResidual() override;
  GenericReal<is_ad> computeScalarQpResidual() override;

  /// Diagonal Jacobian contribution for the out-of-plane strain scalar variable.
  /// Non-AD only: AD captures this through computeScalarQpResidual's AD chain instead.
  Real computeScalarQpJacobian() override;

  /// Fills d-(this instance's field variable)-residual / d-kappa.
  /// Non-AD only: AD captures this through computeScalarQpResidual's AD chain instead.
  Real computeQpOffDiagJacobianScalar(unsigned int svar_num) override;

  /// Fills d-kappa-residual / d-(this instance's field variable).
  /// Non-AD only: AD captures this through computeScalarQpResidual's AD chain instead.
  Real computeScalarQpOffDiagJacobian(unsigned int jvar_num) override;

  /**
   * Non-AD only backward-compatibility hook: when 'reference_residual_excludes_pressure' is set,
   * excludes the applied out-of-plane pressure from the 'Reference' residual tag contribution,
   * matching the exact behavior of the removed non-AD GeneralizedPlaneStrain ScalarKernel. There
   * is no AD equivalent of this method to override (ADKernelScalarBase assembles residual and
   * Jacobian together via automatic differentiation instead), so this is declared without
   * 'override'; it implicitly overrides KernelScalarBase::computeScalarResidual() only when
   * is_ad is false.
   */
  void computeScalarResidual();

  /// Base name of the material system
  const std::string _base_name;

  /// The stress tensor
  const GenericMaterialProperty<RankTwoTensor, is_ad> & _stress;

  /// Function defining the applied out-of-plane pressure
  const Function * const _out_of_plane_pressure_function;

  /// Material property defining the applied out-of-plane pressure
  const MaterialProperty<Real> & _out_of_plane_pressure_material;

  /// Factor applied to the out-of-plane pressure
  const Real _pressure_factor;

  /// Direction of the out-of-plane strain scalar variable
  unsigned int _out_of_plane_direction;

  /// Whether this instance owns the out-of-plane strain scalar variable's own residual and
  /// diagonal Jacobian contributions. Only meaningful for the non-AD case, where multiple
  /// instances (one per in-plane displacement, plus temperature) are coupled to the same scalar
  /// variable; exactly one of them must be primary to avoid contributing its residual/diagonal
  /// Jacobian more than once. Tautological (always true) for AD, which only ever has one instance.
  const bool _primary;

  /// A UserObject that carries the subblock ID for all elements
  const SubblockIndexProvider * const _subblock_id_provider;

  /// The index number of scalar_out_of_plane_strain this kernel acts on
  const unsigned int _scalar_var_id;

  /// Non-AD only: whether to exclude the applied out-of-plane pressure from the 'Reference'
  /// residual tag contribution. See computeScalarResidual().
  const bool _reference_residual_excludes_pressure;

  /// Non-AD only: nullptr for AD, which has no equivalent material property.
  const MaterialProperty<RankFourTensor> * const _Jacobian_mult;

  /// Non-AD only: the displacement variables, used to identify which in-plane component (if any)
  /// this instance's own variable corresponds to.
  std::vector<MooseVariable *> _disp_var;

  /// Non-AD only: the temperature variable, if coupled; nullptr otherwise.
  MooseVariable * const _temp_var;

  /// Non-AD only: derivative of each eigenstrain with respect to temperature.
  std::vector<const MaterialProperty<RankTwoTensor> *> _deigenstrain_dT;

  /// Non-AD only: the in-plane component (0 or 1) that this instance's own variable corresponds
  /// to, or -1 if this instance's own variable is not one of the tracked displacements.
  int _in_plane_component;

private:
  /// Whether the current element belongs to the subblock (region) that this instance's scalar
  /// variable is responsible for. Always true when no subblock_index_provider is given.
  bool onOwnSubblock() const;

  usingGenericKernelScalarMembers;
  using GenericKernelScalar<is_ad>::_grad_test;
  using GenericKernelScalar<is_ad>::_grad_phi;
  using GenericKernelScalar<is_ad>::_JxW;
  using GenericKernelScalar<is_ad>::_coord;
  using GenericKernelScalar<is_ad>::_qrule;
  using GenericKernelScalar<is_ad>::_current_elem;
  using GenericKernelScalar<is_ad>::_assembly;
  using GenericKernelScalar<is_ad>::_tid;
  using GenericKernelScalar<is_ad>::_sys;
  using GenericKernelScalar<is_ad>::_kappa_var_ptr;
};

typedef GeneralizedPlaneStrainTempl<false> GeneralizedPlaneStrain;
typedef GeneralizedPlaneStrainTempl<true> ADGeneralizedPlaneStrain;
