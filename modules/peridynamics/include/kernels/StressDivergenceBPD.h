/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STRESSDIVERGENCEBPD_H
#define STRESSDIVERGENCEBPD_H

#include "Kernel.h"
#include "PeridynamicMesh.h"

class StressDivergenceBPD;

template<>
InputParameters validParams<StressDivergenceBPD>();

class StressDivergenceBPD : public Kernel
{
public:
  StressDivergenceBPD(const InputParameters & parameters);
  virtual ~StressDivergenceBPD() {}

protected:
  virtual void initialSetup();
  virtual void computeResidual();
  virtual Real computeQpResidual() {return 0;}
  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _bond_force_ij;
  const MaterialProperty<Real> & _bond_dfdU_ij;
  const MaterialProperty<Real> & _bond_dfdT_ij;

private:
  AuxiliarySystem & _aux;

  const NumericVector<Number> & _aux_sln;

  const unsigned int _component;

  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

  const bool _temp_coupled;
  const unsigned int _temp_var;

  MooseVariable * _bond_status_var;

  PeridynamicMesh & _pdmesh;

  const std::vector<RealGradient> * _orientation;
};

#endif //STRESSDIVERGENCEBPD_H