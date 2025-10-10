//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <vector>

#include "MooseTypes.h"
#include "XFEM.h"

using namespace libMesh;

namespace libMesh
{
class MeshBase;
class Elem;
class Node;
class QBase;
}
class EFANode;
class EFAElement;

class XFEMCutElem
{
public:
  /**
   * Construct the base XFEM cut element wrapper.
   *
   * @param elem Element on which the cut element is built.
   * @param n_qpoints Number of quadrature points used for integration.
   * @param n_sides Number of sides in the host element.
   */
  XFEMCutElem(Elem * elem, unsigned int n_qpoints, unsigned int n_sides);
  virtual ~XFEMCutElem();

protected:
  unsigned int _n_nodes;
  unsigned int _n_qpoints;
  unsigned int _n_sides;
  std::vector<Node *> _nodes;
  std::vector<Point> _qp_points;
  std::vector<Real> _qp_weights;
  Real _elem_volume;
  std::vector<Real> _elem_side_area;
  Real _physical_volfrac;
  std::vector<Real> _physical_areafrac;
  bool _have_weights;
  std::vector<bool> _have_face_weights;
  /// quadrature weights from volume fraction and moment fitting
  std::vector<Real> _new_weights;
  /// face quadrature weights from surface area fraction
  std::vector<std::vector<Real>> _new_face_weights;
  virtual Point getNodeCoordinates(EFANode * node, MeshBase * displaced_mesh = nullptr) const = 0;

public:
  void setQuadraturePointsAndWeights(const std::vector<Point> & qp_points,
                                     const std::vector<Real> & qp_weights);
  /**
   * Compute the volume fraction of the element fragment.
   */
  virtual void computePhysicalVolumeFraction() = 0;

  /**
   * Get the volume fraction of the element fragment.
   *
   * @return Physical volume fraction stored for the fragment.
   */
  Real getPhysicalVolumeFraction() const;

  /**
   * Compute the surface area fraction of an element side.
   *
   * @param side Side index of the element whose area fraction is requested.
   */
  virtual void computePhysicalFaceAreaFraction(unsigned int side) = 0;

  /**
   * Get the surface area fraction of an element side.
   *
   * @param side Side index of the element whose area fraction is requested.
   * @return Physical surface area fraction of the element side.
   */
  Real getPhysicalFaceAreaFraction(unsigned int side) const;

  virtual void computeMomentFittingWeights() = 0;
  Real getMomentFittingWeight(unsigned int i_qp) const;
  virtual Point getCutPlaneOrigin(unsigned int plane_id,
                                  MeshBase * displaced_mesh = nullptr) const = 0;
  virtual Point getCutPlaneNormal(unsigned int plane_id,
                                  MeshBase * displaced_mesh = nullptr) const = 0;
  virtual void
  getCrackTipOriginAndDirection(unsigned tip_id, Point & origin, Point & direction) const = 0;
  virtual void getFragmentFaces(std::vector<std::vector<Point>> & frag_faces,
                                MeshBase * displaced_mesh = nullptr) const = 0;
  virtual const EFAElement * getEFAElement() const = 0;
  virtual unsigned int numCutPlanes() const = 0;
  void getWeightMultipliers(MooseArray<Real> & weights,
                            QBase * qrule,
                            Xfem::XFEM_QRULE xfem_qrule,
                            const MooseArray<Point> & q_points);
  void getFaceWeightMultipliers(MooseArray<Real> & face_weights,
                                QBase * qrule,
                                Xfem::XFEM_QRULE xfem_qrule,
                                const MooseArray<Point> & q_points,
                                unsigned int side);

  /**
   * Compute integration weights for the cut element.
   *
   * @param qrule Standard MOOSE quadrature rule.
   * @param xfem_qrule Integration scheme for the cut element.
   * @param q_points Quadrature points for the element.
   */
  void computeXFEMWeights(QBase * qrule,
                          Xfem::XFEM_QRULE xfem_qrule,
                          const MooseArray<Point> & q_points);

  /**
   * Compute face integration weights for the cut element side.
   *
   * @param qrule Standard MOOSE face quadrature rule.
   * @param xfem_qrule Integration scheme for the cut element (surface area fraction only).
   * @param q_points Quadrature points for the element side.
   * @param side Side index of the element.
   */
  void computeXFEMFaceWeights(QBase * qrule,
                              Xfem::XFEM_QRULE xfem_qrule,
                              const MooseArray<Point> & q_points,
                              unsigned int side);
  bool isPointPhysical(const Point & p) const;
  virtual void getIntersectionInfo(unsigned int plane_id,
                                   Point & normal,
                                   std::vector<Point> & intersectionPoints,
                                   MeshBase * displaced_mesh = nullptr) const = 0;
};
