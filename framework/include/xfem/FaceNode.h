/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef FACENODE_H
#define FACENODE_H

class EFANode;

class FaceNode
{
public:

  FaceNode(EFANode* node, double xi, double eta);
  FaceNode(const FaceNode & other_face_node);

  ~FaceNode();

private:

  EFANode * _node;
  double _xi;
  double _eta;

public:

  EFANode * getNode();
  double getParametricCoordinates(unsigned int i);
  void switchNode(EFANode* new_old, EFANode* old_node);
};

#endif
