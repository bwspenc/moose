//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"

class XFEM;

/**
 * Coupled auxiliary value
 */
class XFEMMaterialStateMarkerBase : public ElementUserObject
{
public:
  /**
   * Build parameters describing the base XFEM material state marker user object.
   *
   * @return Input parameters defining the material state marking behavior.
   */
  static InputParameters validParams();

  XFEMMaterialStateMarkerBase(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

protected:
  /**
   * Determine whether the current element should be cut by a new crack.
   *
   * @param direction Normal direction of the crack when the element cracks.
   * @return True if the element nucleates a crack.
   */
  virtual bool doesElementCrack(RealVectorValue & direction);

private:
  MooseMesh & _mesh;
  std::vector<BoundaryID> _initiation_boundary_ids;
  bool _secondary_cracks;
  std::shared_ptr<XFEM> _xfem;
  std::map<unsigned int, RealVectorValue> _marked_elems;
  std::set<unsigned int> _marked_frags;
  std::map<unsigned int, unsigned int> _marked_elem_sides;
};
