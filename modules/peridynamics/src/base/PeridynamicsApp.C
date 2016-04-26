#include "PeridynamicsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

#include "BondCriticalStrainPDAux.h"
#include "BondStatusPDAux.h"
#include "FailureIndexPD.h"
#include "FailureIndexPDAux.h"
#include "HeatConductionMaterialPD.h"
#include "HeatConductionPD.h"
#include "HeatSourcePD.h"
#include "LinearIsotropicMaterialPD.h"
#include "StressDivergencePD.h"

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
  registerAux(FailureIndexPDAux);
  registerAux(BondStatusPDAux);
  registerAux(BondCriticalStrainPDAux);

  registerMaterial(HeatConductionMaterialPD);
  registerMaterial(LinearIsotropicMaterialPD);

  registerKernel(HeatConductionPD);
  registerKernel(HeatSourcePD);
  registerKernel(StressDivergencePD);

  registerUserObject(FailureIndexPD);
}

// External entry point for dynamic syntax association
extern "C" void PeridynamicsApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { PeridynamicsApp::associateSyntax(syntax, action_factory); }
void
PeridynamicsApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
