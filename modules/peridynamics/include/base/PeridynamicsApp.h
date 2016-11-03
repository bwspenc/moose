#ifndef PERIDYNAMICSAPP_H
#define PERIDYNAMICSAPP_H

#include "MooseApp.h"

class PeridynamicsApp;

template<>
InputParameters validParams<PeridynamicsApp>();

class PeridynamicsApp : public MooseApp
{
public:
  PeridynamicsApp(InputParameters parameters);
  virtual ~PeridynamicsApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* PERIDYNAMICSAPP_H */
