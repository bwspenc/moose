/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef FILEMESHPD_H
#define FILEMESHPD_H

#include "PeridynamicMesh.h"

class FileMeshPD;

template<>
InputParameters validParams<FileMeshPD>();

class FileMeshPD : public PeridynamicMesh
{
public:
  FileMeshPD(const InputParameters & parameters);

  virtual ~FileMeshPD();

  virtual MooseMesh & clone() const;

  virtual void buildMesh();
};

#endif /* FILEMESHPD_H */
