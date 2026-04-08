//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * Skeleton nodal auxiliary kernel that updates an auxiliary variable from its old value and the
 * current/old values of a coupled variable.
 */
class OldValueCoupledNodalAux : public AuxKernel
{
public:
  static InputParameters validParams();

  OldValueCoupledNodalAux(const InputParameters & parameters);

protected:
  Real computeValue() override;

  /// Old value of the auxiliary variable this kernel computes
  const VariableValue & _u_old;

  /// Current value of the coupled variable
  const VariableValue & _coupled_value;

  /// Old value of the coupled variable
  const VariableValue & _coupled_value_old;
};
