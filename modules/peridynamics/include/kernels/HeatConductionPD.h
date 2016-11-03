/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                       Peridynamics                           */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef HEATCONDUCTIONPD_H
#define HEATCONDUCTIONPD_H

#include "Kernel.h"
#include "PeridynamicMesh.h"

class HeatConductionPD;

template<>
InputParameters validParams<HeatConductionPD>();

class HeatConductionPD : public Kernel
{
public:
  HeatConductionPD(const InputParameters & parameters);
  virtual ~HeatConductionPD();

protected:
  virtual void computeResidual();
  virtual Real computeQpResidual() {return 0;}
  virtual void computeJacobian();

  const MaterialProperty<Real> & _bond_response;
  const MaterialProperty<Real> & _bond_drdT;

  AuxiliarySystem & _aux;

  const NumericVector<Number> & _aux_sln;

  MooseVariable * _bond_status_var;

  PeridynamicMesh & _pdmesh;
};

#endif //HEATCONDUCTIONPD_H
