//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearViscoelasticityManager.h"
#include "libmesh/quadrature.h"

registerMooseObject("TensorMechanicsApp", LinearViscoelasticityManager);

InputParameters
LinearViscoelasticityManager::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addClassDescription("Manages the updating of the semi-implicit "
                             "single-step first-order finite difference "
                             "time-stepping scheme");
  params.addRequiredParam<std::string>("viscoelastic_model",
                                       "name of the LinearViscoelasticityBase object to manage");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_TIMESTEP_BEGIN, EXEC_TIMESTEP_END};
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

LinearViscoelasticityManager::LinearViscoelasticityManager(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _viscoelastic_model_name(getParam<std::string>("viscoelastic_model")),
    _viscoelastic_model(nullptr)
{
}

void
LinearViscoelasticityManager::execute()
{
  if (_mi_feproblem.getCurrentExecuteOnFlag() == EXEC_TIMESTEP_BEGIN)
  {
    for (unsigned int _qp = 0; _qp < _qrule->n_points(); ++_qp)
      _viscoelastic_model->recomputeQpApparentProperties(_qp);
  }
}

void
LinearViscoelasticityManager::initialize()
{
  std::shared_ptr<MaterialBase> test =
      _mi_feproblem.getMaterial(_viscoelastic_model_name, _material_data_type, _mi_tid, true);

  if (!test)
    mooseError(_viscoelastic_model_name + " does not exist");

  _viscoelastic_model = std::dynamic_pointer_cast<LinearViscoelasticityBase>(test);

  if (!_viscoelastic_model)
    mooseError(_viscoelastic_model_name + " is not a LinearViscoelasticityBase object");
}
