//
// Created by Twan Koolen on 12/24/15.
//

#ifndef DRAKE_WRAPPABLESYSTEM_H
#define DRAKE_WRAPPABLESYSTEM_H

#include "Core.h"
#include "System.h"
#include <functional>
#include <memory>


namespace Drake {
  class WrappableSystem {
  public:
    template<typename ScalarType> using InputVector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;
    template<typename ScalarType> using StateVector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;
    template<typename ScalarType> using OutputVector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

    template<typename SystemType>
    WrappableSystem(const std::shared_ptr<SystemType> &system_ptr) :
        dynamics_double([=](const double& t, const StateVector<double> &x, const InputVector<double> &u) {
          auto x_sys = createStateVector<double>(*system_ptr);
          x_sys = x;
          typename SystemType::template InputVector<double> u_sys;
          u_sys = u;
          auto xdot_sys = system_ptr->dynamics(t, x_sys, u_sys);
          return toEigen(xdot_sys);
        }),
        output_double([=](double t, const StateVector<double> &x, const InputVector<double> &u) {
          auto x_sys = createStateVector<double>(*system_ptr);
          x_sys = x;
          typename SystemType::template InputVector<double> u_sys;
          u_sys = u;
          return toEigen(system_ptr->output(t, x_sys, u_sys));
        }),
        isTimeVarying([=] { return system_ptr->isTimeVarying(); }),
        isDirectFeedthrough([=] { return system_ptr->isDirectFeedthrough(); }),
        getNumStates([=] {return Drake::getNumStates(*system_ptr); }),
        getNumInputs([=] {return Drake::getNumInputs(*system_ptr); }),
        getNumOutputs([=] {return Drake::getNumOutputs(*system_ptr); })
    {
      // empty
    }

    StateVector<double> dynamics(const double& t, const StateVector<double>& x, const InputVector<double>& u) const {
      return dynamics_double(t, x, u);
    }

    OutputVector<double> output(const double& t, const StateVector<double>& x, const InputVector<double>& u) const {
      return output_double(t, x, u);
    }

    std::function<bool()> isTimeVarying;
    std::function<bool()> isDirectFeedthrough;

    std::function<size_t()> getNumStates;
    std::function<size_t()> getNumInputs;
    std::function<size_t()> getNumOutputs;

  private:
    std::function<InputVector<double>(const double&, const StateVector<double> &, const InputVector<double> &)> dynamics_double;
    std::function<OutputVector<double>(const double&, const StateVector<double> &, const InputVector<double> &)> output_double;
  };

  template <>
  struct CascadeReturnType<WrappableSystem, WrappableSystem> {
    typedef std::shared_ptr<WrappableSystem> type;
  };

  template <>
  std::shared_ptr<WrappableSystem> cascade(const std::shared_ptr<WrappableSystem> &sys1, const std::shared_ptr<WrappableSystem> &sys2) {
    auto casc_temp = std::make_shared<CascadeSystem<WrappableSystem, WrappableSystem>>(sys1, sys2);
    return std::make_shared<WrappableSystem>(casc_temp);
  };
}

#endif //DRAKE_WRAPPABLESYSTEM_H
