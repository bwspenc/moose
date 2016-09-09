/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                        Peridynamics                          */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef GPSPDDIAG_H
#define GPSPDDIAG_H

#include "ScalarKernel.h"
#include "GPSPDUO.h"

class GPSPDDiag;

template<>
InputParameters validParams<GPSPDDiag>();

class GPSPDDiag : public ScalarKernel
{
public:
  GPSPDDiag(const InputParameters & parameters);

protected:
  virtual void reinit(){};
  virtual void computeResidual();
  virtual void computeJacobian();

  const GPSPDUO & _gps_uo;
};
#endif //GPSPDDIAG_H
