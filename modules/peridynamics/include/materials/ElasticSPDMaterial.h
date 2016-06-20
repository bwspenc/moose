/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELASTICSPDMATERIAL_H
#define ELASTICSPDMATERIAL_H

#include "MechanicPDMaterial.h"

class ElasticSPDMaterial;

template<>
InputParameters validParams<ElasticSPDMaterial>();

class ElasticSPDMaterial : public MechanicPDMaterial
{
public:
  ElasticSPDMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpStrain();
  virtual void computeQpForce();
};

#endif //ELASTICSPDMATERIAL_H
