/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CTHERMALPDMATERIAL_H
#define CTHERMALPDMATERIAL_H

#include "ThermalPDMaterial.h"

class CThermalPDMaterial;

template<>
InputParameters validParams<CThermalPDMaterial>();

class CThermalPDMaterial : public ThermalPDMaterial
{
public:
  CThermalPDMaterial(const InputParameters & parameters);

protected:
  virtual Real computeBondModulus();
};

#endif //CTHERMALPDMATERIAL_H
