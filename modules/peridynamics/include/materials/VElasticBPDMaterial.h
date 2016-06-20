/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef VELASTICBPDMATERIAL_H
#define VELASTICBPDMATERIAL_H

#include "ElasticBPDMaterial.h"

class VElasticBPDMaterial;

template<>
InputParameters validParams<VElasticBPDMaterial>();

class VElasticBPDMaterial : public ElasticBPDMaterial
{
public:
  VElasticBPDMaterial(const InputParameters & parameters);

protected:
  virtual Real computeBondModulus();
};

#endif //VELASTICBPDMATERIAL_H
