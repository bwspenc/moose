//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "USCFluxLimitedDirichletBC.h"
#include "Function.h"

registerMooseObject("MooseApp", USCFluxLimitedDirichletBC);

defineLegacyParams(USCFluxLimitedDirichletBC);

InputParameters
USCFluxLimitedDirichletBC::validParams()
{
  InputParameters params = DirichletBCBase::validParams();
  params.addRequiredParam<FunctionName>("function", "The voltage function.");
  params.addRequiredParam<PostprocessorName>("switch_postprocessor", "PP that switches BCs");
  params.addClassDescription(
      "Imposes the essential boundary condition $u=g(t,\\vec{x})$, where $g$ "
      "is a (possibly) time and space-dependent MOOSE Function.");
  return params;
}

USCFluxLimitedDirichletBC::USCFluxLimitedDirichletBC(const InputParameters & parameters)
  : DirichletBCBase(parameters),
  _func(getFunction("function")),
  _switch_postprocessor(getPostprocessorValue("switch_postprocessor"))
{
}

Real
USCFluxLimitedDirichletBC::computeQpValue()
{
  if (_switch_postprocessor == 1)
    return _func.value(_t, *_current_node);
  else
    mooseError("Shouldn't get here");
}

bool
USCFluxLimitedDirichletBC::shouldApply()
{
  return (_switch_postprocessor == 1);
}
