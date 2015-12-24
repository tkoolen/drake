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
  shared_ptr<WrappableSystemInterface> pendulum_wrapped = make_shared<WrappableSystem>(pendulum);

  auto controller = make_shared<PendulumEnergyShapingController>(*pendulum);
  shared_ptr<WrappableSystemInterface> controller_wrapped = make_shared<WrappableSystem>(controller);

  PendulumState<double> x;
  x.theta = 0.0; x.thetadot = 0.0;
  PendulumInput<double> u;
  u.tau = 1.0;

  auto xdot_wrapped = pendulum_wrapped->dynamics(0.0, toEigen(x), toEigen(u));
  cout << xdot_wrapped.transpose() << endl;

  return 0;
}
