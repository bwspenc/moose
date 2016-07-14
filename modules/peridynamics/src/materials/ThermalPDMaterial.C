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
  _thermal_conductivity(getParam<Real>("thermal_conductivity")),
  _thermal_conductivity_function(getParam<FunctionName>("thermal_conductivity_function") != "" ? &getFunction("thermal_conductivity_function") : NULL),
  _bond_response(declareProperty<Real>("bond_response")),
  _bond_drdT(declareProperty<Real>("bond_drdT"))
{
  NonlinearVariableName temp = parameters.get<NonlinearVariableName>("temp");
  _temp_var = &_fe_problem.getVariable(_tid, temp);
}

void
ThermalPDMaterial::initQpStatefulProperties()
{
}

void
ThermalPDMaterial::computeQpProperties()
{
  //obtain the temperature solution at the two nodes for each truss element
  const NumericVector<Number> & sol = *_nsys.currentSolution();
  unsigned int dof0(_current_elem->get_node(0)->dof_number(_nsys.number(), _temp_var->number(), 0));
  unsigned int dof1(_current_elem->get_node(1)->dof_number(_nsys.number(), _temp_var->number(), 0));
  // the temperature of the connecting bond is calculated as the avarage of the temperature of two end nodes, this value will be used for temperature (and possibly spatial location) dependent thermal conductivity and specific heat calculation
  double temp_avg = (sol(dof0) + sol(dof1)) / 2.0;
  if (_thermal_conductivity_function)
  {
    Point p;
    _kappa = _thermal_conductivity_function->value(temp_avg, p);
  }
  else
    _kappa = _thermal_conductivity;

  double Kij = computeBondModulus();
  if (std::abs(_bond_status[0] - 1.0) < 0.01)
  {
    _bond_response[_qp] = Kij * (sol(dof1) - sol(dof0)) / _origin_length * _nv_i * _nv_j;
    _bond_drdT[_qp] = Kij / _origin_length * _nv_i * _nv_j;
  }
  else
  {
    _bond_response[_qp] = Kij * (sol(dof1) - sol(dof0)) / _origin_length * _nv_i * _nv_j * 0.5;
    _bond_drdT[_qp] = Kij / _origin_length * _nv_i * _nv_j * 0.5;
  }
}
