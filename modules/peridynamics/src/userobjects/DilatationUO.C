/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "DilatationUO.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<DilatationUO>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addRequiredParam<std::vector<NonlinearVariableName> >("displacements", "Variables containing the displacements");
  params.addParam<std::string>("plane_strain", "Plane strain problem");
  params.addRequiredParam<Real>("poissons_ratio", "Poisson's ratio");
  params.addParam<Real>("temp_ref", 273, "Reference temperature in K");
  params.addParam<Real>("thermal_expansion_coeff", 0.0, "Thermal expansion coefficient in 1/K");
  params.addCoupledVar("temp", 273, "Temperature in K");
  params.addCoupledVar("nodal_dilatation", "Auxiliary variable for dilatation for each node");
  params.addCoupledVar("bond_status", "Auxiliary variable for bond failure status");
  params.addCoupledVar("bond_contact", "Auxiliary variable for bond contact status");
  return params;
}

DilatationUO::DilatationUO(const InputParameters & parameters) :
  ElementUserObject(parameters),
  _aux(_fe_problem.getAuxiliarySystem()),
  _pdmesh(dynamic_cast<PeridynamicMesh &>(_mesh)),
  _pddim(_pdmesh.dim()),
  _temp_ref(getParam<Real>("temp_ref")),
  _temp(coupledValue("temp")),
  _nodal_dila_var(getVar("nodal_dilatation", 0)),
  _bond_status(coupledValue("bond_status")),
  _bond_contact(coupledValue("bond_contact"))
{
  const std::vector<NonlinearVariableName> & nl_vnames(getParam<std::vector<NonlinearVariableName> >("displacements"));
  if (_pddim != nl_vnames.size())
    mooseError("Size of displacements vector is different from the mesh dimension!");
  for (unsigned int i = 0; i < _pddim; ++i)
    _disp_var.push_back(&_fe_problem.getVariable(_tid, nl_vnames[i]));

  if (isParamValid("plane_strain"))
    _alpha = (1.0 + getParam<Real>("poissons_ratio")) * getParam<Real>("thermal_expansion_coeff");
  else
    _alpha = getParam<Real>("thermal_expansion_coeff");
}

void
DilatationUO::initialize()
{
  std::vector<std::string> zero_vars;
  zero_vars.push_back("nodal_dilatation");
  _aux.zeroVariables(zero_vars);
}

void
DilatationUO::execute()
{
  NonlinearSystem & non_sys = _fe_problem.getNonlinearSystem();
  const NumericVector<Number> & current_sol = *non_sys.currentSolution();

  std::vector<std::vector<Real> > disp(2);
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < _pddim; ++j)
      disp[i].push_back(current_sol(_current_elem->get_node(i)->dof_number(non_sys.number(), _disp_var[j]->number(), 0)));
  // origin length of a truss element
  RealGradient dxyz;
  for (unsigned int i = 0; i < _pddim; ++i)
    dxyz(i) = (*_current_elem->get_node(1))(i) - (*_current_elem->get_node(0))(i);
  double origin_length = dxyz.norm();
  // current length of a truss element
  for (unsigned int i = 0; i < _pddim; ++i)
    dxyz(i) += disp[1][i] - disp[0][i];
  double current_length = dxyz.norm();
  double val = current_length / origin_length - 1.0 - _alpha * ((_temp[0] + _temp[1]) / 2.0 - _temp_ref);

  dof_id_type node0 = _current_elem->get_node(0)->id();
  dof_id_type node1 = _current_elem->get_node(1)->id();
  // dof index for each nodal aux variable
  long int nd_dof0 = _current_elem->get_node(0)->dof_number(_aux.number(), _nodal_dila_var->number(), 0);
  long int nd_dof1 = _current_elem->get_node(1)->dof_number(_aux.number(), _nodal_dila_var->number(), 0);

  if (std::abs(_bond_status[0] - 1.0) < 0.01)
  {
    _aux.solution().add(nd_dof0, val * _pdmesh.volume(node1));
    _aux.solution().add(nd_dof1, val * _pdmesh.volume(node0));
  }
  else
  {
    _aux.solution().add(nd_dof0, 0.0);
    _aux.solution().add(nd_dof1, 0.0);
  }
}

void
DilatationUO::threadJoin(const UserObject & uo)
{
}

void
DilatationUO::finalize()
{
  _aux.solution().close();
}
