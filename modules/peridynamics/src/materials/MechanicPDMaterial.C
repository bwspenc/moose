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
  params.addParam<std::string>("plane_strain", "Plane strain problem");
  params.addRequiredParam<Real>("youngs_modulus", "Young's modulus");
  params.addRequiredParam<Real>("poissons_ratio", "Poisson's ratio");

  params.addCoupledVar("temp", 273, "Temperature in K");
  params.addParam<Real>("temp_ref", 273, "Reference temperature in K");
  params.addParam<Real>("thermal_expansion_coeff", 0.0, "Thermal expansion coefficient in 1/K");
  return params;
}

MechanicPDMaterial::MechanicPDMaterial(const InputParameters & parameters) :
  PeridynamicMaterial(parameters),
  _bond_total_strain(declareProperty<Real>("bond_total_strain")),
  _bond_mechanic_strain(declareProperty<Real>("bond_mechanic_strain")),
  _bond_elastic_strain(declareProperty<Real>("bond_elastic_strain")),

  _youngs_modulus(getParam<Real>("youngs_modulus")),
  _poissons_ratio(getParam<Real>("poissons_ratio")),

  _temp(coupledValue("temp")),
  _temp_ref(getParam<Real>("temp_ref"))
{
  const std::vector<NonlinearVariableName> & nl_vnames(getParam<std::vector<NonlinearVariableName> >("displacements"));
  if (_pddim != nl_vnames.size())
    mooseError("Size of displacements vector is different from the mesh dimension!");
  // fetch nonlinear variables
  for (unsigned int i = 0; i < _pddim; ++i)
    _disp_var.push_back(&_fe_problem.getVariable(_tid, nl_vnames[i]));

  _shear_modulus = _youngs_modulus / 2.0 / (1.0 + _poissons_ratio);
  if (isParamValid("plane_strain"))
  {
    _bulk_modulus = _youngs_modulus / 2.0 / (1.0 + _poissons_ratio) / (1.0 - 2.0 * _poissons_ratio);
    _alpha = (1.0 + _poissons_ratio) * getParam<Real>("thermal_expansion_coeff");
  }
  else
  {
    _bulk_modulus = _youngs_modulus / _pddim / (1.0 - (_pddim - 1.0) * _poissons_ratio);
    _alpha = getParam<Real>("thermal_expansion_coeff");
  }
}

void
MechanicPDMaterial::initQpStatefulProperties()
{
  _bond_total_strain[_qp] = 0.0;
  _bond_mechanic_strain[_qp] = 0.0;
  _bond_elastic_strain[_qp] = 0.0;
}

void
MechanicPDMaterial::computeQpProperties()
{
  computeQpStrain();
  computeQpForce();
}

Real
MechanicPDMaterial::computeBondCurrentLength()
{
  // solution at the two end nodes
  const NumericVector<Number> & sol = *_nsys.currentSolution();
  std::vector<std::vector<Real> > disp(2);
  RealGradient dxyz;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < _pddim; ++j)
      disp[i].push_back(sol(_current_elem->get_node(i)->dof_number(_nsys.number(), _disp_var[j]->number(), 0)));
  // current length of a truss element
  for (unsigned int i = 0; i < _pddim; ++i)
    dxyz(i) = (*_current_elem->get_node(1))(i) + disp[1][i] - (*_current_elem->get_node(0))(i) - disp[0][i];

  return dxyz.norm();
}

void
MechanicPDMaterial::computeQpStrain()
{
  _bond_total_strain[_qp] = _current_length / _origin_length - 1.0;
  // bond temperature is taken as the average of two end nodes
  _bond_mechanic_strain[_qp] = _bond_total_strain[_qp] - _alpha * ((_temp[0] + _temp[1]) / 2.0 - _temp_ref);
  _bond_elastic_strain[_qp] = _bond_mechanic_strain[_qp];
}
