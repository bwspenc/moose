/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CELASTICSPDMATERIAL_H
#define CELASTICSPDMATERIAL_H

#include "ElasticSPDMaterial.h"

class CElasticSPDMaterial;

template<>
InputParameters validParams<CElasticSPDMaterial>();

class CElasticSPDMaterial : public ElasticSPDMaterial
{
public:
  CElasticSPDMaterial(const InputParameters & parameters);

protected:
  virtual Real computeBondModulus();
};

#endif //CELASTICSPDMATERIAL_H
