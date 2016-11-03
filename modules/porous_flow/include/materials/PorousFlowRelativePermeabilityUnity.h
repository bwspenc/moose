/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWRELATIVEPERMEABILITYUNITY_H
#define POROUSFLOWRELATIVEPERMEABILITYUNITY_H

#include "PorousFlowRelativePermeabilityBase.h"

class PorousFlowRelativePermeabilityUnity;

template<>
InputParameters validParams<PorousFlowRelativePermeabilityUnity>();

/**
 * This class simply sets relative permeability = 1 at the nodes
 */
class PorousFlowRelativePermeabilityUnity : public PorousFlowRelativePermeabilityBase
{
public:
  PorousFlowRelativePermeabilityUnity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
};

#endif //POROUSFLOWRELATIVEPERMEABILITYUNITY_H
