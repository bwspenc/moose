/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PeridynamicsAction.h"
#include "Factory.h"
#include "FEProblem.h"
#include "Conversion.h"

template<>
InputParameters validParams<PeridynamicsAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription("Set up peridynamic stress divergence kernels. Default is BPD. For SPD, partial Jacobian is the default");
  params.addRequiredParam<std::vector<NonlinearVariableName> >("displacements", "The nonlinear displacement variables for the problem");
  params.addParam<std::string>("state_based_formulation", "Whether to use SPD or not");
  params.addParam<std::string>("full_jacobian", "whether to use full SPD jacobian or not");
  params.addParam<NonlinearVariableName>("temp", "The temperature variable");
  params.addParam<NonlinearVariableName>("strain_zz", "The strain_zz variable");
  params.addParam<bool>("use_displaced_mesh", true, "Whether to use displaced mesh in the kernels");
  params.addParam<std::vector<SubdomainName> >("block", "The list of ids of the blocks (subdomain) that the stress divergence kernel will be applied to");
  params.addParam<std::vector<AuxVariableName> >("save_in", "The displacement residuals");
  params.addParam<std::vector<AuxVariableName> >("diag_save_in", "The displacement diagonal preconditioner terms");
  return params;
}

PeridynamicsAction::PeridynamicsAction(const InputParameters & params) :
  Action(params)
{
}

void
PeridynamicsAction::act()
{
  std::vector<NonlinearVariableName>displacements = getParam<std::vector<NonlinearVariableName> >("displacements");
  unsigned int _ndisp = displacements.size();

  std::vector<VariableName>coupled_displacements;
  for (unsigned int i = 0; i < _ndisp; ++i)
    coupled_displacements.push_back(displacements[i]);

  std::vector<AuxVariableName>save_in = getParam<std::vector<AuxVariableName> >("save_in");
  if (isParamValid("save_in") && save_in.size() != _ndisp)
    mooseError("Number of save_in variables should equal to the number of displacement variables " << _ndisp);

  std::vector<AuxVariableName>diag_save_in = getParam<std::vector<AuxVariableName> >("diag_save_in");
  if (isParamValid("diag_save_in") && diag_save_in.size() != _ndisp)
    mooseError("Number of diag_save_in variables should equal to the number of displacement variables " << _ndisp);

  std::string type;
  if (isParamValid("state_based_formulation"))
    type = "StressDivergenceBPD";
  else
    type = "StressDivergenceSPD";

  InputParameters params = _factory.getValidParams(type);
  params.set<std::vector<VariableName> >("displacements") = coupled_displacements;

  if (isParamValid("state_based_formulation") && isParamValid("full_jacobian"))
    params.set<std::string>("full_jacobian") = getParam<std::string>("full_jacobian");

  if (isParamValid("temp"))
    params.addCoupledVar("temp", "The temperature variable");

  if (isParamValid("strain_zz"))
    params.addCoupledVar("strain_zz", "The strain_zz variable");

  params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");

  // check whether this kernel is restricted to certain block?
  if (isParamValid("block"))
    params.set<std::vector<SubdomainName> >("block") = getParam<std::vector<SubdomainName> >("block");

  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    std::string kernel_name = "Peridynamics_" + Moose::stringify(i);

    params.set<unsigned int>("component") = i;
    params.set<NonlinearVariableName>("variable") = displacements[i];
    if (isParamValid("save_in"))
      params.set<std::vector<AuxVariableName> >("save_in") = std::vector<AuxVariableName>(1, save_in[i]);
    if (isParamValid("diag_save_in"))
      params.set<std::vector<AuxVariableName> >("diag_save_in") = std::vector<AuxVariableName>(1, diag_save_in[i]);

    addkernel(kernel_name, params);
  }
}

void
PeridynamicsAction::addkernel(const std::string & name,  InputParameters & params)
{
  if (isParamValid("state_based_formulation"))
    _problem->addKernel("StressDivergenceSPD", name, params);
  else
    _problem->addKernel("StressDivergenceBPD", name, params);
}
