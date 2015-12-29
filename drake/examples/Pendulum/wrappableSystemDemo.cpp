//
// Created by Twan Koolen on 12/24/15.
//

#include "WrappableSystem.h"
#include "Pendulum.h"
#include <iostream>

using namespace std;
using namespace Drake;

int main() {

  auto pendulum = make_shared<Pendulum>();
  shared_ptr<WrappableSystem> pendulum_wrapped = make_shared<WrappableSystem>(pendulum);

  auto controller = make_shared<PendulumEnergyShapingController>(*pendulum);
  shared_ptr<WrappableSystem> controller_wrapped = make_shared<WrappableSystem>(controller);

  shared_ptr<WrappableSystem> cascade_wrapped = cascade(pendulum_wrapped, controller_wrapped);

  PendulumState<double> x;
  x.theta = 0.2; x.thetadot = 0.3;
  PendulumInput<double> u;
  u.tau = 1.0;

  auto xdot_wrapped = pendulum_wrapped->dynamics(0.0, toEigen(x), toEigen(u));
  cout << xdot_wrapped.transpose() << endl;

  auto cascade_output = cascade_wrapped->output(0.0, toEigen(x), toEigen(u));
  cout << cascade_output.transpose() << endl;

  return 0;
}
