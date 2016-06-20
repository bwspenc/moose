/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                      Peridynamics                            */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELASTICBPDMATERIAL_H
#define ELASTICBPDMATERIAL_H

#include "MechanicPDMaterial.h"

class ElasticBPDMaterial;

template<>
InputParameters validParams<ElasticBPDMaterial>();

class ElasticBPDMaterial : public MechanicPDMaterial
{
public:
  ElasticBPDMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpStrain();
  virtual void computeQpForce();
};

#endif //ELASTICBPDMATERIAL_H
