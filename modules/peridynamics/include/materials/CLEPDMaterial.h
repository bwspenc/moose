/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CLEPDMATERIAL_H
#define CLEPDMATERIAL_H

#include "MechanicPDMaterial.h"

class CLEPDMaterial;

template<>
InputParameters validParams<CLEPDMaterial>();

class CLEPDMaterial : public MechanicPDMaterial
{
public:
  CLEPDMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpStrain();
  virtual void computeQpForce();
};

#endif //CLEPDMATERIAL_H
