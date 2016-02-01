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

#ifndef ELEMENTFRAGMENTALGORITHM_H
#define ELEMENTFRAGMENTALGORITHM_H

#include <cstddef>
#include <iostream>
#include <sstream>

#include <vector>
#include <map>
#include <set>
#include <limits>

class EFANode;
class EFAElement;

class ElementFragmentAlgorithm
{
public:

  /**
   * Constructor
   **/
  ElementFragmentAlgorithm();

  ~ElementFragmentAlgorithm();

private:
  //unsigned int MaxElemId;
  std::map< unsigned int, EFANode*> _permanent_nodes;
  std::map< unsigned int, EFANode*> _embedded_nodes;
  std::map< unsigned int, EFANode*> _temp_nodes;
  std::map< unsigned int, EFAElement*> _elements;
//  std::map< std::set< EFAnode* >, std::set< EFAelement* > > _merged_edge_map;
  std::set< EFAElement*> _crack_tip_elements;
  std::vector< EFANode* > _new_nodes;
  std::vector< EFAElement* > _child_elements;
  std::vector< EFAElement* > _parent_elements;
  std::map< EFANode*, std::set< EFAElement *> > _inverse_connectivity;

public:

  unsigned int add2DElements( std::vector< std::vector<unsigned int> > &quads );
  EFAElement* add2DElement( std::vector<unsigned int> quad, unsigned int id );
  EFAElement* add3DElement( std::vector<unsigned int> quad, unsigned int id );

  void updateEdgeNeighbors();
  void initCrackTipTopology();
  void addElemEdgeIntersection(unsigned int elemid, unsigned int edgeid, double position);
  void addFragEdgeIntersection(unsigned int elemid, unsigned int frag_edge_id, double position);
  void addElemFaceIntersection(unsigned int elemid, unsigned int faceid,
                               std::vector<unsigned int> edgeid, std::vector<double> position);
  void addFragFaceIntersection(unsigned int ElemID, unsigned int FragFaceID,
                               std::vector<unsigned int> FragFaceEdgeID, std::vector<double> position);

  void updatePhysicalLinksAndFragments();

  void updateTopology(bool mergeUncutVirtualEdges=true);
  void reset();
  void clearAncestry();
  void restoreFragmentInfo(EFAElement * const elem, const EFAElement * const from_elem);

  void createChildElements();
  void connectFragments(bool mergeUncutVirtualEdges);

  void sanityCheck();
  void updateCrackTipElements();
  void printMesh();
  void error(const std::string &error_string);

  const std::vector<EFAElement*> &getChildElements(){return _child_elements;};
  const std::vector<EFAElement*> &getParentElements(){return _parent_elements;};
  const std::vector<EFANode*> &getNewNodes(){return _new_nodes;};
  const std::set<EFAElement*> &getCrackTipElements(){return _crack_tip_elements;};
  const std::map<unsigned int, EFANode*> &getPermanentNodes(){return _permanent_nodes;}
  const std::map<unsigned int, EFANode*> &getTempNodes(){return _temp_nodes;}
  const std::map<unsigned int, EFANode*> &getEmbeddedNodes(){return _embedded_nodes;}
  EFAElement* getElemByID(unsigned int id);
  unsigned int getElemIdByNodes(unsigned int * node_id);
  void clearPotentialIsolatedNodes();
};

#endif // #ifndef ELEMENTFRAGMENTALGORITHM_H
