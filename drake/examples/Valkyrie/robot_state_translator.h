#pragma once

// TODO: find a better place for this

#include "drake/common/eigen_types.h"
#include "drake/systems/lcm/lcm_and_vector_base_translator.h"
#include "drake/systems/plants/RigidbodyTree.h"
#include "lcmtypes/bot_core/robot_state_t.hpp"

namespace drake {

class RobotStateTranslator: public systems::lcm::LcmAndVectorBaseTranslator {
 public:
  RobotStateTranslator(const RigidBodyTree &tree);

  virtual void TranslateLcmToVectorBase(
      const void *lcm_message_bytes, int lcm_message_length,
      systems::VectorBase<double> *vector_base) const override;

  virtual void TranslateVectorBaseToLcm(
      const systems::VectorBase<double> &vector_base,
      std::vector<uint8_t> *lcm_message_bytes) const override;

 private:
  void InitializeMessage() const;

  Eigen::Isometry3d EvalFloatingBodyPose(
      const Eigen::Ref<const Eigen::VectorXd> &q) const;

  TwistVector<double> EvalFloatingBodyTwistInBodyFrame(
      const Eigen::Ref<const Eigen::VectorXd> &q,
      const Eigen::Ref<const Eigen::VectorXd> &v) const;

  TwistVector<double> TransformTwistFromBodyFrameToWorldAlignedBodyFrame(
      Eigen::Isometry3d &floating_body_to_world,
      const TwistVector<double> &twist_in_body) const;

  int num_floating_joint_positions() const;

  int num_floating_joint_velocities() const;
  const RigidBodyTree &tree_;
  const RigidBody *const floating_body_;

  mutable bot_core::robot_state_t message_;
  const RigidBodyTree &CheckPreConditions(const RigidBodyTree &tree);
};

}  // drake
