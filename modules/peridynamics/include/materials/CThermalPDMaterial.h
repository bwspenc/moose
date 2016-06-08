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

#include "PeridynamicMaterial.h"

class CThermalPDMaterial;
class Function;

template<>
InputParameters validParams<CThermalPDMaterial>();

class CThermalPDMaterial : public PeridynamicMaterial
{
public:
  CThermalPDMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  const Real _thermal_conductivity;
  Function * _thermal_conductivity_function;

  MaterialProperty<Real> & _bond_response;
  MaterialProperty<Real> & _bond_drdT;

  MooseVariable * _temp_var;
};

#endif //CTHERMALPDMATERIAL_H
