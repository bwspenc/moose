//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "XFEMCutElem.h"
#include "EFAElement3D.h"

using namespace libMesh;

namespace libMesh
{
class MeshBase;
class Elem;
class Node;
}

class XFEMCutElem3D : public XFEMCutElem
{
public:
  /**
   * Construct an XFEM cut element for a three-dimensional host element.
   *
   * @param elem Element on which the cut element is built.
   * @param CEMelem Fragment representing the cut element topology.
   * @param n_qpoints Number of quadrature points used for integration.
   * @param n_sides Number of sides in the host element.
   */
  XFEMCutElem3D(Elem * elem,
                const EFAElement3D * const CEMelem,
                unsigned int n_qpoints,
                unsigned int n_sides);
  ~XFEMCutElem3D();

private:
  EFAElement3D _efa_elem3d; // 3D EFAelement
  virtual Point getNodeCoordinates(EFANode * node, MeshBase * displaced_mesh = nullptr) const;

public:
  virtual void computePhysicalVolumeFraction();
  virtual void computePhysicalFaceAreaFraction(unsigned int side);
  virtual void computeMomentFittingWeights();
  virtual Point getCutPlaneOrigin(unsigned int plane_id, MeshBase * displaced_mesh = nullptr) const;
  virtual Point getCutPlaneNormal(unsigned int plane_id, MeshBase * displaced_mesh = nullptr) const;
  virtual void
  getCrackTipOriginAndDirection(unsigned tip_id, Point & origin, Point & direction) const;
  virtual void getFragmentFaces(std::vector<std::vector<Point>> & frag_faces,
                                MeshBase * displaced_mesh = nullptr) const;
  virtual const EFAElement * getEFAElement() const;
  virtual unsigned int numCutPlanes() const;
  virtual void getIntersectionInfo(unsigned int plane_id,
                                   Point & normal,
                                   std::vector<Point> & intersectionPoints,
                                   MeshBase * displaced_mesh = nullptr) const;
};
