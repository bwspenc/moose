/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MechanicPDMaterial.h"
#include "NonlinearSystem.h"

// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<MechanicPDMaterial>()
{
  InputParameters params = validParams<PeridynamicMaterial>();
  params.addRequiredParam<std::vector<NonlinearVariableName> >("displacements", "Variables containing the displacements");
  params.addParam<std::string>("plane_stress", "Plane stress problem or not");
  params.addRequiredParam<Real>("youngs_modulus", "Young's modulus");
  params.addRequiredParam<Real>("poissons_ratio", "Poisson's ratio");
  params.addCoupledVar("strain_zz", "Variable containing the out-of-plane strain");
  params.addParam<NonlinearVariableName>("temp", "Variable containing the temperature");
  params.addRequiredParam<NonlinearVariableName>("bond_status", "Auxiliary variable for failure status of each bond");
  params.addParam<Real>("temp_ref", 273, "Reference temperature in K");
  params.addParam<Real>("thermal_expansion_coeff", 0.0, "Thermal expansion coefficient in 1/K");
  return params;
}

MechanicPDMaterial::MechanicPDMaterial(const InputParameters & parameters) :
  PeridynamicMaterial(parameters),
  _bond_elastic_strain(declareProperty<Real>("bond_elastic_strain")),
  _bond_critical_strain(declareProperty<Real>("bond_critical_strain")),
  _bond_critical_strain_old(declarePropertyOld<Real>("bond_critical_strain")),
  _elasticity_tensor(declareProperty<RankFourTensor>("elasticity_tensor")),
  _thermal_expansion(declareProperty<Real>("thermal_expansion")),
  _shape_tensor(declareProperty<RankTwoTensor>("shape_tensor")),
  _deformation_gradient(declareProperty<RankTwoTensor>("deformation_gradient")),
  _elastic_strain(declareProperty<RankTwoTensor>("elastic_strain")),
  _strain(declareProperty<RankTwoTensor>("strain")),
  _stress(declareProperty<RankTwoTensor>("stress")),
  _bond_force_ij(declareProperty<Real>("bond_force_ij")),
  _bond_force_i_j(declareProperty<Real>("bond_force_i_j")),
  _bond_dfdU_ij(declareProperty<Real>("bond_dfdU_ij")),
  _bond_dfdU_i_j(declareProperty<Real>("bond_dfdU_i_j")),
  _bond_dfdE_ij(declareProperty<Real>("bond_dfdE_ij")),
  _bond_dfdE_i_j(declareProperty<Real>("bond_dfdE_i_j")),
  _bond_dfdT_ij(declareProperty<Real>("bond_dfdT_ij")),
  _bond_dfdT_i_j(declareProperty<Real>("bond_dfdT_i_j")),
  _youngs_modulus(getParam<Real>("youngs_modulus")),
  _poissons_ratio(getParam<Real>("poissons_ratio")),
  _strain_zz_coupled(isCoupledScalar("strain_zz")),
  _strain_zz(coupledScalarValue("strain_zz")),
  _temp_var(isParamValid("temp") ? &_fe_problem.getVariable(_tid, getParam<NonlinearVariableName>("temp")) : NULL),
  _bond_status_var(&_fe_problem.getVariable(_tid, getParam<NonlinearVariableName>("bond_status"))),
  _temp_ref(getParam<Real>("temp_ref"))
{
  const std::vector<NonlinearVariableName> & nl_vnames(getParam<std::vector<NonlinearVariableName> >("displacements"));
  if (_pddim != nl_vnames.size())
    mooseError("Size of displacements vector is different from the mesh dimension!");

  for (unsigned int i = 0; i < _pddim; ++i)
    _disp_var.push_back(&_fe_problem.getVariable(_tid, nl_vnames[i]));

  // define shear, bulk moduli
  _shear_modulus = _youngs_modulus / 2.0 / (1.0 + _poissons_ratio);
  if (_pddim == 3 || isParamValid("plane_stress"))
  {
    _bulk_modulus = _youngs_modulus / _pddim / (1.0 - (_pddim - 1.0) * _poissons_ratio);
    _alpha = isParamValid("thermal_expansion_coeff") ? getParam<Real>("thermal_expansion_coeff") : 0;
  }
  else // plane strain case
  {
    _bulk_modulus = _youngs_modulus / 2.0 / (1.0 + _poissons_ratio) / (1.0 - 2.0 * _poissons_ratio);
   if (_strain_zz_coupled) // generalized plane strain
     _alpha = isParamValid("thermal_expansion_coeff") ? getParam<Real>("thermal_expansion_coeff") : 0;
   else
     _alpha = isParamValid("thermal_expansion_coeff") ? (1.0 + _poissons_ratio) * getParam<Real>("thermal_expansion_coeff") : 0;
  }

  // construct elasticity tensor
  std::vector<Real> iso_const(2);
  iso_const[0] = _youngs_modulus * _poissons_ratio / ((1.0 + _poissons_ratio) * (1.0 - 2.0 * _poissons_ratio));
  iso_const[1] = _youngs_modulus / (2.0 * (1.0 + _poissons_ratio));
  _Cijkl.fillFromInputVector(iso_const, RankFourTensor::symmetric_isotropic);

  setRandomResetFrequency(EXEC_INITIAL);
}

void
MechanicPDMaterial::initQpStatefulProperties()
{
  // Generate randomized critical stretch by Box-Muller method: randomized_critical_strain = small_random_numb    er + _input_critical_strain
//  double val = std::sqrt(2.0 * _Gc / _kappa / std::pow(_pddim, 2) / _current_elem_volume);
//  double val = std::sqrt(8.0 * 3.1415926 * _Gc / 27.0 / _E / 0.17697842); //3x
//  double val = std::sqrt(8.0 * 3.1415926 * _Gc / 27.0 / _E / 0.235971227); //4x
//  double val = std::sqrt(8.0 * 3.1415926 * _Gc / 27.0 / _E / 0.294964033); //5x
//  double val = 0.000446185; // 3x irregular
//  double val = 0.000389091; // 4x irregular
//  double val = 0.00034595; // 5x irregular
//  double val = 0.000592505; // 3x regular
//  double val = 0.000464875; // 4x regular
//  double val = 0.00038582; // 5x regular
//  double val = 0.000300093; // 3x irregular 3D
//  double val = 0.000233381; // 5x irregular 3D

  _bond_critical_strain_old[_qp] = (std::sqrt(- 2.0 * std::log(getRandomReal())) * std::cos(2.0 * 3.14159265358 * getRandomReal()) * 0.05 + 1.0) * 0.00034595;
  _bond_critical_strain[_qp] = _bond_critical_strain_old[_qp];
}

void
MechanicPDMaterial::computeQpStrain()
{
  _elasticity_tensor[_qp] = _Cijkl;
  _thermal_expansion[_qp] = _alpha;

  _bond_elastic_strain[_qp] = (_current_length / _origin_length - 1.0) - _alpha * ((_temp_i + _temp_j) / 2.0 - _temp_ref);
}

void
MechanicPDMaterial::computeNodalTemp()
{
  if(isParamValid("temp"))
  {
    const NumericVector<Number> & nsys_sln = *_nsys.currentSolution();
    _temp_i = nsys_sln(_current_elem->get_node(0)->dof_number(_nsys.number(), _temp_var->number(), 0));
    _temp_j = nsys_sln(_current_elem->get_node(1)->dof_number(_nsys.number(), _temp_var->number(), 0));
  }
  else
  {
    _temp_i = _temp_ref;
    _temp_j = _temp_ref;
  }
}

Real
MechanicPDMaterial::computeBondCurrentLength()
{
  const NumericVector<Number> & nsys_sln = *_nsys.currentSolution();
  RealGradient dxyz;
  dof_id_type dof0, dof1;

  for (unsigned int i = 0; i < _pddim; ++i)
  {
    dof0 = _current_elem->get_node(0)->dof_number(_nsys.number(), _disp_var[i]->number(), 0);
    dof1 = _current_elem->get_node(1)->dof_number(_nsys.number(), _disp_var[i]->number(), 0);
    dxyz(i) = (*_current_elem->get_node(1))(i) + nsys_sln(dof1) - (*_current_elem->get_node(0))(i) - nsys_sln(dof0);
  }

  return dxyz.norm();
}

void
MechanicPDMaterial::computeElasticStrainTensor()
{
  const NumericVector<Number> & nsys_sln = *_nsys.currentSolution();
  const NumericVector<Number> & aux_sln = *_aux.currentSolution();

  Node * node_i = _current_elem->get_node(0);
  Node * node_j = _current_elem->get_node(1);

  std::vector<dof_id_type> neighbors_i, bonds_i, neighbors_j, bonds_j;

  neighbors_i = _pdmesh.neighbors(node_i->id());
  bonds_i = _pdmesh.bonds(node_i->id());
  neighbors_j = _pdmesh.neighbors(node_j->id());
  bonds_j = _pdmesh.bonds(node_j->id());

  RankTwoTensor delta;
  delta.zero();
  delta(0, 0) = delta(1, 1) = delta(2, 2) = 1;

  // for node i
  // calculate the origin and current shape tensors at node i
  RankTwoTensor current_shape_i;
  current_shape_i.zero();
  _shape_tensor[0].zero();
  if (_pddim == 2)
    _shape_tensor[0](2, 2) = current_shape_i(2, 2) = 1;

  RealGradient origin_vector(_pddim), current_vector(_pddim);

  for (unsigned int k = 0; k < neighbors_i.size(); ++k)
  {
    Node * node_k = _mesh.nodePtr(neighbors_i[k]);
    double vol_k = _pdmesh.volume(neighbors_i[k]);
    for (unsigned int l = 0; l < _pddim; ++l)
    {
      origin_vector(l) = _pdmesh.coord(neighbors_i[k])(l) - _pdmesh.coord(node_i->id())(l);
      current_vector(l) = origin_vector(l) + nsys_sln(node_k->dof_number(_nsys.number(), _disp_var[l]->number(), 0)) - nsys_sln(node_i->dof_number(_nsys.number(), _disp_var[l]->number(), 0));
    }
    double origin_length = origin_vector.norm();

   // bond status for bond ik
    Elem * elem_k = _mesh.elemPtr(bonds_i[k]);
    dof_id_type bs_dof_ik = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status_ik = aux_sln(bs_dof_ik);

    for (unsigned int m = 0; m < _pddim; ++m)
      for (unsigned int n = 0; n < _pddim; ++n)
      {
        _shape_tensor[0](m, n) += vol_k * _horizon_i / origin_length * origin_vector(m) * origin_vector(n) * bond_status_ik;
        current_shape_i(m, n) += vol_k * _horizon_i / origin_length * current_vector(m) * origin_vector(n) * bond_status_ik;
      }
   }

   if (std::abs(_shape_tensor[0].det()) > 1e-30)
   {
    // inverse the origin shape tensor at node i
    _shape_tensor[0] = _shape_tensor[0].inverse();
    // calculate the deformation gradient tensor at node i
    _deformation_gradient[0] = current_shape_i * _shape_tensor[0];
   }
   else
   {
     _shape_tensor[0] = delta;
     _deformation_gradient[0] = delta;
   }

  // the green-lagrange strain tensor at node i
  _strain[0] = (_deformation_gradient[0].transpose() * _deformation_gradient[0] - delta) / 2.0;

  // the elastic strain tensor at node i
  if (_strain_zz_coupled)
    _strain[0](2, 2) = _strain_zz[0];

  // the elastic strain tensor at node j
  _elastic_strain[0] = _strain[0];
  _elastic_strain[0].addIa(- _alpha * (_temp_i - _temp_ref));

  // for node j
  // calculate the origin and current shape tensors at node j
  RankTwoTensor current_shape_j;
  current_shape_j.zero();
  _shape_tensor[1].zero();
  if (_pddim == 2)
    _shape_tensor[1](2, 2) = current_shape_j(2, 2) = 1;

  for (unsigned int k = 0; k < neighbors_j.size(); ++k)
  {
    Node * node_k = _mesh.nodePtr(neighbors_j[k]);
    double vol_k = _pdmesh.volume(neighbors_j[k]);
    for (unsigned int l = 0; l < _pddim; ++l)
    {
      origin_vector(l) = _pdmesh.coord(neighbors_j[k])(l) - _pdmesh.coord(node_j->id())(l);
      current_vector(l) = origin_vector(l) + nsys_sln(node_k->dof_number(_nsys.number(), _disp_var[l]->number(), 0)) - nsys_sln(node_j->dof_number(_nsys.number(), _disp_var[l]->number(), 0));
    }
    double origin_length = origin_vector.norm();

   // bond status for bond jk
    Elem * elem_k = _mesh.elemPtr(bonds_j[k]);
    dof_id_type bs_dof_jk = elem_k->dof_number(_aux.number(), _bond_status_var->number(), 0);
    Number bond_status_jk = aux_sln(bs_dof_jk);

    for (unsigned int m = 0; m < _pddim; ++m)
      for (unsigned int n = 0; n < _pddim; ++n)
      {
        _shape_tensor[1](m, n) += vol_k * _horizon_j / origin_length * origin_vector(m) * origin_vector(n) * bond_status_jk;
        current_shape_j(m, n) += vol_k * _horizon_j / origin_length * current_vector(m) * origin_vector(n) * bond_status_jk;
      }
   }

   if (std::abs(_shape_tensor[1].det()) > 1e-30)
   {
    // inverse the origin shape tensor at node j
    _shape_tensor[1] = _shape_tensor[1].inverse();
    // calculate the deformation gradient tensor at node j
    _deformation_gradient[1] = current_shape_j * _shape_tensor[1];
   }
   else
   {
     _shape_tensor[1] = delta;
     _deformation_gradient[1] = delta;
   }

  // the green-lagrange strain tensor at node j
  _strain[1] = (_deformation_gradient[1].transpose() * _deformation_gradient[1] - delta) / 2.0;

  if (_strain_zz_coupled)
    _strain[1](2, 2) = _strain_zz[0];

  // the elastic strain tensor at node j
  _elastic_strain[1] = _strain[1];
  _elastic_strain[1].addIa(- _alpha * (_temp_j - _temp_ref));
}

void
MechanicPDMaterial::computeStressTensor()
{
  _stress[0] = _Cijkl * _elastic_strain[0];
  _stress[1] = _Cijkl * _elastic_strain[1];
}
