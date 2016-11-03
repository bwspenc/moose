/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                         Peridynamics                         */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef GPSPDOFFDIAG_H
#define GPSPDOFFDIAG_H

#include "Kernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "PeridynamicMesh.h"

class GPSPDOffDiag;

template<>
InputParameters validParams<GPSPDOffDiag>();

class GPSPDOffDiag : public Kernel
{
public:
  GPSPDOffDiag(const InputParameters & parameters);
  virtual ~GPSPDOffDiag() {}

protected:
  virtual void computeResidual(){};
  virtual Real computeQpResidual(){return 0;}
  virtual void computeJacobian(){};
  virtual void computeOffDiagJacobian(unsigned int jvar){};
  virtual void computeOffDiagJacobianScalar(unsigned int jvar);
  virtual void computeDispFullOffDiagJacobianScalar(unsigned int component, unsigned int jvar);
  virtual void computeDispPartialOffDiagJacobianScalar(unsigned int component, unsigned int jvar);
  virtual void computeTempOffDiagJacobianScalar(unsigned int jvar);

  const MaterialProperty<Real> & _bond_dfdE_ij;
  const MaterialProperty<Real> & _bond_dfdE_i_j;

  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<RankFourTensor> & _Cijkl;
  const MaterialProperty<RankTwoTensor> & _shape;
  const MaterialProperty<RankTwoTensor> & _dgrad;

private:
  AuxiliarySystem & _aux;
  const NumericVector<Number> & _aux_sln;

  NonlinearSystem & _nsys;

  MooseVariable * _temp_var;

  const unsigned int _strain_zz_var;

  MooseVariable * _bond_status_var;

  PeridynamicMesh & _pdmesh;
  const unsigned int _pddim;

  std::vector<MooseVariable *> _disp_var;
};

#endif //GPSPDOFFDIAG_H
