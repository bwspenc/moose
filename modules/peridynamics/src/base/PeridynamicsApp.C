/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

// Default
#include "PeridynamicsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

// Actions
#include "PeridynamicsAction.h"

// AuxKernels
#include "BondCriticalStrainAux.h"
#include "BondContactStrainAux.h"
#include "BondStatusAux.h"
#include "BondContactAux.h"
#include "FailureIndexAux.h"
#include "NodalVolSumAux.h"
#include "NodalDilatationAux.h"

// UserObjects
#include "FailureIndexUO.h"
#include "PDNodalUO.h"
#include "DilatationUO.h"

// Materials
#include "CElasticBPDMaterial.h"
#include "CElasticSPDMaterial.h"
#include "VElasticBPDMaterial.h"
#include "VElasticSPDMaterial.h"
#include "CThermalPDMaterial.h"
#include "VThermalPDMaterial.h"

//Kernels
#include "StressDivergencePD.h"
#include "HeatConductionPD.h"
#include "HeatSourcePD.h"

template<>
InputParameters validParams<PeridynamicsApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

PeridynamicsApp::PeridynamicsApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  PeridynamicsApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  PeridynamicsApp::associateSyntax(_syntax, _action_factory);
}

PeridynamicsApp::~PeridynamicsApp()
{
}

// External entry point for dynamic application loading
extern "C" void PeridynamicsApp__registerApps() { PeridynamicsApp::registerApps(); }
void
PeridynamicsApp::registerApps()
{
  registerApp(PeridynamicsApp);
}

// External entry point for dynamic object registration
extern "C" void PeridynamicsApp__registerObjects(Factory & factory) { PeridynamicsApp::registerObjects(factory); }
void
PeridynamicsApp::registerObjects(Factory & factory)
{
  registerAux(FailureIndexAux);
  registerAux(BondStatusAux);
  registerAux(BondContactAux);
  registerAux(BondCriticalStrainAux);
  registerAux(BondContactStrainAux);
  registerAux(NodalVolSumAux);
  registerAux(NodalDilatationAux);

  registerMaterial(CThermalPDMaterial);
  registerMaterial(VThermalPDMaterial);
  registerMaterial(CElasticBPDMaterial);
  registerMaterial(CElasticSPDMaterial);
  registerMaterial(VElasticBPDMaterial);
  registerMaterial(VElasticSPDMaterial);

  registerKernel(HeatConductionPD);
  registerKernel(HeatSourcePD);
  registerKernel(StressDivergencePD);

  registerUserObject(FailureIndexUO);
  registerUserObject(PDNodalUO);
  registerUserObject(DilatationUO);
}

// External entry point for dynamic syntax association
extern "C" void PeridynamicsApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { PeridynamicsApp::associateSyntax(syntax, action_factory); }
void
PeridynamicsApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  syntax.registerActionSyntax("PeridynamicsAction", "Kernels/Peridynamics");

  registerAction(PeridynamicsAction, "add_kernel");
}
