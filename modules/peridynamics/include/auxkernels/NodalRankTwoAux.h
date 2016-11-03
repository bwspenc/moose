/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                        Peridynamics                          */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef NODALRANKTWOAUX_H
#define NODALRANKTWOAUX_H

#include "AuxKernel.h"
#include "PeridynamicMesh.h"
#include "NonlinearSystem.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class NodalRankTwoAux;

template<>
InputParameters validParams<NodalRankTwoAux>();

class NodalRankTwoAux : public AuxKernel
{
public:
  NodalRankTwoAux(const InputParameters & parameters);
  virtual ~NodalRankTwoAux() {}

protected:
  virtual Real computeValue();
  virtual RankTwoTensor computeNodalStrain();
  virtual RankTwoTensor computeNodalStress();

  FEProblem & _fe_problem;
  NonlinearSystem & _nsys;
  PeridynamicMesh & _pdmesh;
  const unsigned int _pddim;

  std::string _rank_two_tensor;
  double _youngs_modulus;
  double _poissons_ratio;
  MooseEnum _quantity_type;

  unsigned int _i;
  unsigned int _j;

  const Point _point1;
  const Point _point2;
  Point _direction;

  std::vector<MooseVariable *> _disp_var;

  RankFourTensor _Cijkl;
};

#endif //NODALRANKTWOAUX_H
