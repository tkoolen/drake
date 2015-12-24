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
//  namespace internal {
//    template <typename System, typename Scalar>
//    class WrappableSystemDynamics {
//    private:
//      typename System::template StateVector<Scalar> x_sys;
//      typename System::template InputVector<Scalar> u_sys;
//      std::function<typename System::template StateVector<Scalar> >
//
//    public:
//      StateVector<Scalar> dynamics(Scalar t, const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& x, const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& u) {
//
//      }
//    };
//  }

  class WrappableSystemInterface {
  public:
    template<typename ScalarType> using InputVector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;
    template<typename ScalarType> using StateVector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;
    template<typename ScalarType> using OutputVector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

    virtual ~WrappableSystemInterface() {};
    virtual StateVector<double> dynamics(const double& t, const StateVector<double>& x, const InputVector<double>& u) const = 0;
    virtual bool isTimeVarying() const = 0;
    virtual bool isDirectFeedthrough() const = 0;
  };

  class WrappableSystem : public WrappableSystemInterface {
  public:
    template<typename SystemType>
    WrappableSystem(const std::shared_ptr<SystemType> &system_ptr) :
        dynamics_double([=](const double& t, const StateVector<double> &x, const InputVector<double> &u) {
          typename SystemType::template StateVector<double> x_sys;
          x_sys = x;
          typename SystemType::template InputVector<double> u_sys;
          u_sys = u;
          return toEigen(system_ptr->dynamics(t, x_sys, u_sys));
        }),
//      output([=](double t, const StateVector<double> &x) { return system_ptr->output(t, x); }),
        is_time_varying([=] { return system_ptr->isTimeVarying(); }),
        is_direct_feedthrough([=] { return system_ptr->isDirectFeedthrough(); }) {
      // empty
    }

    virtual ~WrappableSystem() {};

    virtual StateVector<double> dynamics(const double& t, const StateVector<double>& x, const InputVector<double>& u) const override {
      return dynamics_double(t, x, u);
    }

    virtual bool isTimeVarying() const override {
      return is_time_varying();
    }

    virtual bool isDirectFeedthrough() const override {
      return is_direct_feedthrough();
    }

  private:
    std::function<InputVector<double>(const double&, const StateVector<double> &, const InputVector<double> &)> dynamics_double;
    std::function<OutputVector<double>(const StateVector<double> &)> output;

    // sys info
    std::function<bool()> is_time_varying;
    std::function<bool()> is_direct_feedthrough;
  };

  template <>
  struct CascadeReturnType<WrappableSystemInterface, WrappableSystemInterface> {
    typedef std::shared_ptr<WrappableSystemInterface> type;
  };

  template<>
  std::shared_ptr<WrappableSystemInterface> cascade(const std::shared_ptr<WrappableSystemInterface> &sys1, const std::shared_ptr<WrappableSystemInterface> &sys2) {
    auto casc_temp = std::make_shared<CascadeSystem<WrappableSystemInterface, WrappableSystemInterface>>(sys1, sys2);
    return std::make_shared<WrappableSystem>(casc_temp);
  };
}

#endif //DRAKE_WRAPPABLESYSTEM_H
