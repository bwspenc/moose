/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef VTHERMALPDMATERIAL_H
#define VTHERMALPDMATERIAL_H

#include "ThermalPDMaterial.h"

class VThermalPDMaterial;

template<>
InputParameters validParams<VThermalPDMaterial>();

class VThermalPDMaterial : public ThermalPDMaterial
{
public:
  VThermalPDMaterial(const InputParameters & parameters);

protected:
  virtual Real computeBondModulus();
};

#endif //VTHERMALPDMATERIAL_H
