/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CThermalPDMaterial.h"
#include "Material.h"
#include "Function.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<CThermalPDMaterial>()
{
  InputParameters params = validParams<PeridynamicMaterial>();
  params.addRequiredParam<NonlinearVariableName>("temp", "Variable containing the temperature");
  params.addParam<Real>("thermal_conductivity", 0, "The thermal conductivity value");
  params.addParam<FunctionName>("thermal_conductivity_function", "", "Thermal conductivity as a function of temperature");
  return params;
}

CThermalPDMaterial::CThermalPDMaterial(const InputParameters & parameters) :
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
CThermalPDMaterial::initQpStatefulProperties()
{
}

void
CThermalPDMaterial::computeQpProperties()
{
  //obtain the temperature solution at the two nodes for each truss element
  NonlinearSystem & nonlinear_sys = _fe_problem.getNonlinearSystem();
  const NumericVector<Number>& ghosted_solution = *nonlinear_sys.currentSolution();
  unsigned int temp_dof0(_current_elem->get_node(0)->dof_number(nonlinear_sys.number(), _temp_var->number(), 0));
  unsigned int temp_dof1(_current_elem->get_node(1)->dof_number(nonlinear_sys.number(), _temp_var->number(), 0));

  Real temp_node0 = ghosted_solution(temp_dof0);
  Real temp_node1 = ghosted_solution(temp_dof1);

// the temperature of the connecting bond is calculated as the avarage of the temperature of two end nodes, this value will be used for temperature (and possibly spatial location) dependent thermal conductivity and specific heat calculation
  Real temp_avg = (temp_node0 + temp_node1) / 2.0;

  double kappa;
  if (_thermal_conductivity_function)
  {
    Point p;
    kappa = _thermal_conductivity_function->value(temp_avg, p);
  }
  else
    kappa = _thermal_conductivity;

  if (_pddim == 2)
  {
    _bond_response[_qp] = 6.0 * kappa / (3.14159265358 * std::pow(_horizon, 3)) * (temp_node1 - temp_node0) / _origin_length * _nv_i * _nv_j;
    _bond_drdT[_qp] = 6.0 * kappa / (3.14159265358 * std::pow(_horizon, 3)) / _origin_length * _nv_i * _nv_j;
  }
  else if (_pddim == 3)
  {
    _bond_response[_qp] = 6.0 * kappa / (3.14159265358 * std::pow(_horizon, 4)) * (temp_node1 - temp_node0) / _origin_length * _nv_i * _nv_j;
    _bond_drdT[_qp] = 6.0 * kappa / (3.14159265358 * std::pow(_horizon, 4)) / _origin_length * _nv_i * _nv_j;
  }
}
