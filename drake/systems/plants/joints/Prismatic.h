#ifndef DRAKE_PRISMATIC_H
#define DRAKE_PRISMATIC_H


#include "drake/systems/plants/joints/FixedAxisOneDoF.h"

namespace Drake {

template <typename J>
class Prismatic : public FixedAxisOneDoF<J> {
 public:
  using FixedAxisOneDoF<J>::getJointAxis;

  Prismatic(const Eigen::Matrix<J, 3, 1>& translation_axis) : FixedAxisOneDoF<J>(createSpatialAxis(translation_axis)) { }

  template <typename DerivedQ>
  Transform3D<Promote<J, typename DerivedQ::Scalar>> jointTransform(const Eigen::MatrixBase<DerivedQ>& q) const {
    using Q = typename DerivedQ::Scalar;
    using T = Promote<J, Q>;

    Transform3D<T> ret;
    ret.linear().setIdentity();
    const auto& joint_axis = getJointAxis();
    ret.translation() = Convert<T>()(q[0]) * ConvertMatrix<T>()(joint_axis.template bottomRows<3>());
    ret.makeAffine();
    return ret;
  }

 private:
  static SpatialVector<J> createSpatialAxis(const Eigen::Matrix<J, 3, 1>& translation_axis) {
    auto ret = SpatialVector<J>::Zero().eval();
    ret.template bottomRows<3>() = translation_axis;
    return ret;
  }

};

}

#endif //DRAKE_PRISMATIC_H
