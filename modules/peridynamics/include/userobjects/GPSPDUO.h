/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef GPSPDUO_H
#define GPSPDUO_H

#include "ElementUserObject.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "PeridynamicMesh.h"

class GPSPDUO;

template<>
InputParameters validParams<GPSPDUO>();

class GPSPDUO : public ElementUserObject
{
public:
  GPSPDUO(const InputParameters & parameters);
  virtual ~GPSPDUO() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & uo){};
  virtual void finalize(){};
  virtual Real returnResidual() const;
  virtual Real returnJacobian() const;

protected:
  const MaterialProperty<RankFourTensor> & _Cijkl;
  const MaterialProperty<RankTwoTensor> & _stress;

  PeridynamicMesh & _pdmesh;

private:
  Real _residual;

  Real _jacobian;
};

#endif // GPSPDUO_H
