/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CELASTICBPDMATERIAL_H
#define CELASTICBPDMATERIAL_H

#include "ElasticBPDMaterial.h"

class CElasticBPDMaterial;

template<>
InputParameters validParams<CElasticBPDMaterial>();

class CElasticBPDMaterial : public ElasticBPDMaterial
{
public:
  CElasticBPDMaterial(const InputParameters & parameters);

protected:
  virtual Real computeBondModulus();
};

#endif //CELASTICBPDMATERIAL_H
