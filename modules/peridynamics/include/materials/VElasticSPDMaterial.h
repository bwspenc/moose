/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef VELASTICSPDMATERIAL_H
#define VELASTICSPDMATERIAL_H

#include "ElasticSPDMaterial.h"

class VElasticSPDMaterial;

template<>
InputParameters validParams<VElasticSPDMaterial>();

class VElasticSPDMaterial : public ElasticSPDMaterial
{
public:
  VElasticSPDMaterial(const InputParameters & parameters);

protected:
  virtual Real computeBondModulus();
};

#endif //VELASTICSPDMATERIAL_H
