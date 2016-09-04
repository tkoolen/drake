#pragma once

// TODO: find a better place for this

#include "drake/systems/lcm/lcm_and_vector_base_translator.h"
#include "drake/systems/plants/RigidbodyTree.h"
#include "drake/common/eigen_types.h"
#include "lcmtypes/bot_core/robot_state_t.hpp"

namespace drake {

class RobotStateTranslator : public systems::lcm::LcmAndVectorBaseTranslator {
 public:
  RobotStateTranslator(const RigidBodyTree& tree);

  virtual void TranslateLcmToVectorBase(
      const void* lcm_message_bytes, int lcm_message_length,
      systems::VectorBase<double>* vector_base) const override;

  virtual void TranslateVectorBaseToLcm(
      const systems::VectorBase<double>& vector_base,
      std::vector<uint8_t>* lcm_message_bytes) const override;

 private:
  static const RigidBody* GetFloatingBody(const RigidBodyTree& tree);

  // TODO: move to better place.
  template <typename DestScalar, typename Derived>
  static void EigenVectorToStdVector(std::vector<DestScalar> &dest,
                         const Eigen::MatrixBase<Derived>& src) {
    static_assert(Derived::ColsAtCompileTime == 1, "src must be a vector");
    dest.clear();
    dest.resize(static_cast<size_t>(src.size()));
    for (Eigen::Index i = 0; i < src.size(); i++) {
      dest.push_back(static_cast<DestScalar>(src[i]));
    }
  }

  Eigen::Isometry3d EvalFloatingBodyPose(
      const Eigen::Ref<const Eigen::VectorXd>& q) const;

  TwistVector<double> EvalFloatingBodyTwistInBodyFrame(
      const Eigen::Ref<const Eigen::VectorXd>& q,
      const Eigen::Ref<const Eigen::VectorXd>& v) const;

  int num_floating_joint_positions() const;

  int num_floating_joint_velocities() const;

  const RigidBodyTree& tree_;
  const RigidBody* const floating_body_;
  mutable bot_core::robot_state_t message_;

};

}  // drake
