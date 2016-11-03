/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*                          Peridynamics                        */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef GENERATEDMESHPD_H
#define GENERATEDMESHPD_H

#include "PeridynamicMesh.h"

class GeneratedMeshPD;

template<>
InputParameters validParams<GeneratedMeshPD>();

class GeneratedMeshPD : public PeridynamicMesh
{
public:
  GeneratedMeshPD(const InputParameters & parameters);
  virtual ~GeneratedMeshPD();

  virtual MooseMesh & clone() const;

  virtual void buildMesh();

  virtual void build2DRectangular(UnstructuredMesh & mesh, BoundaryInfo & boundary_info);
  virtual void build2DCircular(UnstructuredMesh & mesh, BoundaryInfo & boundary_info);
  virtual void build3DRectangular(UnstructuredMesh & mesh, BoundaryInfo & boundary_info);
  virtual void build3DCylindrical(UnstructuredMesh & mesh, BoundaryInfo & boundary_info);

protected:
  // Number of elements in x, or R direction
  unsigned int _nx, _shape;
  // domain size in x, y, z, R direction
  double _xmin, _ymin, _zmin, _xmax, _ymax, _zmax, _R;
};

#endif /* GENERATEDMESHPD_H */
