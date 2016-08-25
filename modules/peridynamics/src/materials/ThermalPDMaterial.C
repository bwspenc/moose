/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ThermalPDMaterial.h"
#include "Material.h"
#include "NonlinearSystem.h"
#include "Function.h"

template<>
InputParameters validParams<ThermalPDMaterial>()
{
  InputParameters params = validParams<PeridynamicMaterial>();
  params.addRequiredParam<NonlinearVariableName>("temp", "Variable containing the temperature");
  params.addParam<Real>("thermal_conductivity", 0, "The thermal conductivity value");
  params.addParam<FunctionName>("thermal_conductivity_function", "", "Thermal conductivity as a function of temperature");
  return params;
}

ThermalPDMaterial::ThermalPDMaterial(const InputParameters & parameters) :
  PeridynamicMaterial(parameters),
  _temp_var(&_fe_problem.getVariable(_tid, getParam<NonlinearVariableName>("temp"))),
  _thermal_conductivity(getParam<Real>("thermal_conductivity")),
  _thermal_conductivity_function(getParam<FunctionName>("thermal_conductivity_function") != "" ? &getFunction("thermal_conductivity_function") : NULL),
  _bond_response(declareProperty<Real>("bond_response")),
  _bond_drdT(declareProperty<Real>("bond_drdT"))
{
}

void
ThermalPDMaterial::computeQpForce()
{
  if (_thermal_conductivity_function)
  {
    Point p;
    _kappa = _thermal_conductivity_function->value((_temp_i + _temp_j) / 2.0, p);
  }
  else
    _kappa = _thermal_conductivity;

  double Kij = computeBondModulus();

  _bond_response[_qp] = Kij * (_temp_j - _temp_i) / _origin_length * _nv_i * _nv_j;
  _bond_drdT[_qp] = Kij / _origin_length * _nv_i * _nv_j;
}

void
ThermalPDMaterial::computeNodalTemp()
{
  const NumericVector<Number> & sol = *_nsys.currentSolution();
  _temp_i = sol(_current_elem->get_node(0)->dof_number(_nsys.number(), _temp_var->number(), 0));
  _temp_j = sol(_current_elem->get_node(1)->dof_number(_nsys.number(), _temp_var->number(), 0));
}
