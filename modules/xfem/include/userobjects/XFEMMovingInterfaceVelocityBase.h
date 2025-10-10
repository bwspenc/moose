//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DiscreteElementUserObject.h"
#include "NodeValueAtXFEMInterface.h"

class XFEMMovingInterfaceVelocityBase : public DiscreteElementUserObject
{
public:
  /**
   * Build the parameter set for XFEM moving interface velocity user objects.
   *
   * @return Input parameters defining how interface velocities are computed.
   */
  static InputParameters validParams();

  XFEMMovingInterfaceVelocityBase(const InputParameters & parameters);
  virtual ~XFEMMovingInterfaceVelocityBase() {}

  /**
   * Initialize state prior to evaluating interface velocities during the current timestep.
   */
  virtual void initialize() override;

  /**
   * Compute the interface velocity for a node on the cut mesh.
   *
   * @param node_id Identifier of the node whose velocity is being computed.
   * @param normal Normal direction at the node defining positive motion.
   * @return Interface velocity magnitude projected along the supplied normal.
   */
  virtual Real computeMovingInterfaceVelocity(dof_id_type node_id,
                                              RealVectorValue normal) const = 0;

  /**
   * Compute total number of nodes that are used to define an interface
   */
  unsigned int numberNodes() const { return _value_at_interface_uo->numberNodes(); }

protected:
  /// Pointer to NodeValueAtXFEMInterface object
  const NodeValueAtXFEMInterface * _value_at_interface_uo;
};
