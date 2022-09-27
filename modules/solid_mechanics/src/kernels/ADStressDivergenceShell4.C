/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ADStressDivergenceShell4.h"

// MOOSE includes
#include "Assembly.h"
#include "Material.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "RankTwoTensor.h"
#include "NonlinearSystem.h"
#include "MooseMesh.h"
#include "ArbitraryQuadrature.h"
#include "DenseMatrix.h"

#include "libmesh/quadrature.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/string_to_enum.h"

registerMooseObject("TensorMechanicsApp", ADStressDivergenceShell4);

InputParameters
ADStressDivergenceShell4::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Quasi-static stress divergence kernel for Shell element");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "component",
      "component < 6",
      "An integer corresponding to the degree of freedom "
      "this kernel acts on. (0 for disp_x, "
      "1 for disp_y, 2 for disp_z, 3 for rot_x, 4 for rot_y, 5 for rot_z)");
  params.addRequiredParam<std::string>("through_thickness_order",
                                       "Quadrature order in out of plane direction");
  params.addParam<bool>(
      "large_strain", false, "Set to true to turn on finite strain calculations.");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<Real>(
      "penalty", 1e4, "Penalty parameter for out of plane stress");
  return params;
}

ADStressDivergenceShell4::ADStressDivergenceShell4(const InputParameters & parameters)
  : ADKernel(parameters),
    _component(getParam<unsigned int>("component")),
    _large_strain(getParam<bool>("large_strain")),
    _penalty(getParam<Real>("penalty")),
    _global_local_rotation(getADMaterialProperty<RankTwoTensor>("global_local_rotation"))
{
  _t_qrule = std::make_unique<QGauss>(
      1, Utility::string_to_enum<Order>(getParam<std::string>("through_thickness_order")));
  _t_weights = _t_qrule->get_weights();

  _stress.resize(_t_weights.size());
  _stress_old.resize(_t_weights.size());
  _B_mat.resize(_t_weights.size());
  if (_large_strain)
    _B_nl.resize(_t_weights.size());
  _J_map.resize(_t_weights.size());
  _gamma_x.resize(_t_weights.size());
  _gamma_y.resize(_t_weights.size());
  _gamma_z.resize(_t_weights.size());
  _gamma.resize(_t_weights.size());

  for (unsigned int i = 0; i < _t_weights.size(); ++i)
  {
    _stress[i] = &getADMaterialProperty<RankTwoTensor>("stress_t_points_" + std::to_string(i));
    _stress_old[i] =
        &getMaterialPropertyOldByName<RankTwoTensor>("stress_t_points_" + std::to_string(i));
    _B_mat[i] = &getADMaterialProperty<DenseMatrix<Real>>("B_t_points_" + std::to_string(i));
    if (_large_strain)
      _B_nl[i] = &getADMaterialProperty<DenseMatrix<Real>>("B_nl_t_points_" + std::to_string(i));

    _J_map[i] = &getADMaterialProperty<Real>("J_mapping_t_points_" + std::to_string(i));
    _gamma_x[i] =
        &getADMaterialProperty<Real>("gamma_test_x_t_points_" + std::to_string(i));
    _gamma_y[i] =
        &getADMaterialProperty<Real>("gamma_test_y_t_points_" + std::to_string(i));
    _gamma_z[i] =
        &getADMaterialProperty<Real>("gamma_test_z_t_points_" + std::to_string(i));
    _gamma[i] =
        &getADMaterialProperty<Real>("gamma_test_t_points_" + std::to_string(i));
  }
}

ADReal
ADStressDivergenceShell4::computeQpResidual()
{
  _q_weights = _qrule->get_weights();
  ADReal residual = 0.0;
  const unsigned int cart_comp = _component < 3 ? _component : _component - 3;

  for (_qp_z = 0; _qp_z < _t_weights.size(); ++_qp_z)
  {
    ADRankTwoTensor _local_global_rotation = _global_local_rotation[_qp].transpose();
    ADReal residual_global = 0.0;
//    std::cout<<"BWS local_rotation: "<<std::endl;
//    _local_global_rotation.printReal();
//    if (_component < 3)
    if (false)
    {
      ADReal residual_loc =
        (*_stress[_qp_z])[_qp](0, 0) * (*_B_mat[_qp_z])[_qp](0, _i + _component * 4) +
        (*_stress[_qp_z])[_qp](1, 1) * (*_B_mat[_qp_z])[_qp](1, _i + _component * 4) +
        2.0 * (*_stress[_qp_z])[_qp](0, 1) * (*_B_mat[_qp_z])[_qp](2, _i + _component * 4) +
        2.0 * (*_stress[_qp_z])[_qp](0, 2) * (*_B_mat[_qp_z])[_qp](3, _i + _component * 4) +
        2.0 * (*_stress[_qp_z])[_qp](1, 2) * (*_B_mat[_qp_z])[_qp](4, _i + _component * 4);

      if (_large_strain)
        residual_loc +=
          (*_stress_old[_qp_z])[_qp](0, 0) * (*_B_nl[_qp_z])[_qp](0, _i + _component * 4) +
          (*_stress_old[_qp_z])[_qp](1, 1) * (*_B_nl[_qp_z])[_qp](1, _i + _component * 4) +
          2.0 * (*_stress_old[_qp_z])[_qp](0, 1) * (*_B_nl[_qp_z])[_qp](2, _i + _component * 4) +
          2.0 * (*_stress_old[_qp_z])[_qp](0, 2) * (*_B_nl[_qp_z])[_qp](3, _i + _component * 4) +
          2.0 * (*_stress_old[_qp_z])[_qp](1, 2) * (*_B_nl[_qp_z])[_qp](4, _i + _component * 4);
      residual_global = residual_loc;
    }
    else
    {
      for (unsigned int comp5 = 0; comp5 < 5; ++comp5)
      {
        const unsigned int local_comp = comp5 < 3 ? comp5 : comp5 - 3;
        ADReal residual_loc =
          (*_stress[_qp_z])[_qp](0, 0) * (*_B_mat[_qp_z])[_qp](0, _i + comp5 * 4) +
          (*_stress[_qp_z])[_qp](1, 1) * (*_B_mat[_qp_z])[_qp](1, _i + comp5 * 4) +
          2.0 * (*_stress[_qp_z])[_qp](0, 1) * (*_B_mat[_qp_z])[_qp](2, _i + comp5 * 4) +
          2.0 * (*_stress[_qp_z])[_qp](0, 2) * (*_B_mat[_qp_z])[_qp](3, _i + comp5 * 4) +
          2.0 * (*_stress[_qp_z])[_qp](1, 2) * (*_B_mat[_qp_z])[_qp](4, _i + comp5 * 4);

        if (_large_strain)
          residual_loc +=
            (*_stress_old[_qp_z])[_qp](0, 0) * (*_B_nl[_qp_z])[_qp](0, _i + comp5 * 4) +
            (*_stress_old[_qp_z])[_qp](1, 1) * (*_B_nl[_qp_z])[_qp](1, _i + comp5 * 4) +
            2.0 * (*_stress_old[_qp_z])[_qp](0, 1) * (*_B_nl[_qp_z])[_qp](2, _i + comp5 * 4) +
            2.0 * (*_stress_old[_qp_z])[_qp](0, 2) * (*_B_nl[_qp_z])[_qp](3, _i + comp5 * 4) +
            2.0 * (*_stress_old[_qp_z])[_qp](1, 2) * (*_B_nl[_qp_z])[_qp](4, _i + comp5 * 4);
//        std::cout<<"BWS component: "<<_component<<" cart_comp: "<<cart_comp<<" comp5: "<<comp5<<" local_comp: "<<local_comp<<std::endl;
//        std::cout<<"BWS rloc: "<<MetaPhysicL::raw_value(residual_loc)<<
//          " rglob: "<<MetaPhysicL::raw_value(_local_global_rotation(cart_comp, local_comp) * residual_loc) <<std::endl;
        residual_global += _local_global_rotation(cart_comp, local_comp) * residual_loc;
      }
    }
//    if (_component > 2  && _i == _qp)
//      residual_global += _local_global_rotation(cart_comp, 2) * _penalty * (*_gamma[_qp_z])[_qp] / (_ad_JxW[_qp] * _ad_coord[_qp]);

    if(_component == 3 && _i == _qp)
      residual_global += _penalty * (*_gamma_x[_qp_z])[_qp];// / (_ad_JxW[_qp] * _ad_coord[_qp]);
    if(_component == 4 && _i == _qp)
      residual_global += _penalty * (*_gamma_y[_qp_z])[_qp];// / (_ad_JxW[_qp] * _ad_coord[_qp]);
    if(_component == 5 && _i == _qp)
      residual_global += _penalty * (*_gamma_z[_qp_z])[_qp];// / (_ad_JxW[_qp] * _ad_coord[_qp]);

    residual += residual_global * (*_J_map[_qp_z])[_qp] * _q_weights[_qp] * _t_weights[_qp_z] /
                (_ad_JxW[_qp] * _ad_coord[_qp]);
  }
  return residual;
}
