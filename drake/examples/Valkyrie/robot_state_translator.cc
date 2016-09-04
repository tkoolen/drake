#include "drake/examples/valkyrie/robot_state_translator.h"
#include "drake/util/lcmUtil.h"

using drake::systems::lcm::LcmAndVectorBaseTranslator;
using Eigen::Dynamic;
using Eigen::Index;
using Eigen::Matrix;
using Eigen::Isometry3d;

namespace drake {

RobotStateTranslator::RobotStateTranslator(const RigidBodyTree& tree)
    : LcmAndVectorBaseTranslator(tree.number_of_positions() +
                                 tree.number_of_velocities()),
      tree_(tree),
      floating_body_(GetFloatingBody(tree)) {
  message_.num_joints = static_cast<int16_t>(tree_.number_of_velocities() -
      num_floating_joint_velocities());

}

void RobotStateTranslator::TranslateLcmToVectorBase(
    const void* lcm_message_bytes, int lcm_message_length,
    systems::VectorBase<double>* vector_base) const {}

void RobotStateTranslator::TranslateVectorBaseToLcm(
    const systems::VectorBase<double>& vector_base,
    std::vector<uint8_t>* lcm_message_bytes) const {
  //  message_.utime = TODO

  auto x = vector_base.get_value();
  auto q = x.head(tree_.number_of_positions());
  auto v = x.tail(tree_.number_of_velocities());

  Isometry3d pose = EvalFloatingBodyPose(q);
  EncodePose(pose, message_.pose);

  TwistVector<double> twist = EvalFloatingBodyTwistInBodyFrame(q, v);
  EncodeTwist(twist, message_.twist);

  // set the joint_position, joint_velocity, joint_effort parts
  auto q_non_floating = q.tail(q.size() - num_floating_joint_positions());
  EigenVectorToStdVector(message_.joint_position, q_non_floating);

  int num_floating_joint_velocities =
  auto v_non_floating = v.tail(v.size() - num_floating_joint_velocities());
  EigenVectorToStdVector(message_.joint_velocity, v_non_floating);


  const int lcm_message_length = message_.getEncodedSize();
  lcm_message_bytes->resize(static_cast<size_t>(lcm_message_length));
  message_.encode(lcm_message_bytes, 0, lcm_message_length);
}

int RobotStateTranslator::num_floating_joint_velocities() const {
  int num_floating_joint_velocities =
      floating_body_ ? floating_body_->getJoint().getNumVelocities() : 0;
  return num_floating_joint_velocities;
}

int RobotStateTranslator::num_floating_joint_positions() const {
  int num_floating_joint_positions =
      floating_body_ ? floating_body_->getJoint().getNumPositions() : 0;
  return num_floating_joint_positions;
}

TwistVector<double> RobotStateTranslator::EvalFloatingBodyTwistInBodyFrame(
    const Eigen::Ref<const Eigen::VectorXd>& q,
    const Eigen::Ref<const Eigen::VectorXd>& v) const {
  TwistVector<double> twist;
  if (floating_body_) {
    const auto& joint = floating_body_->getJoint();

    int q_start_index = floating_body_->get_position_start_index();
    auto q_body = q.middleRows(q_start_index, joint.getNumPositions());
    int v_start_index = floating_body_->get_velocity_start_index();
    auto v_body = v.middleRows(v_start_index, joint.getNumVelocities());

    Matrix<double, kTwistSize, Dynamic, 0, kTwistSize,
           DrakeJoint::MAX_NUM_VELOCITIES>
        motion_subspace;
    Matrix<double, Dynamic, Dynamic>* dmotion_subspacedq = nullptr;
    joint.motionSubspace(q_body, motion_subspace, dmotion_subspacedq);
    twist = motion_subspace * v_body;
  } else {
    twist.setZero();
  }
  return twist;
}

Eigen::Isometry3d RobotStateTranslator::EvalFloatingBodyPose(
    const Eigen::Ref<const Eigen::VectorXd>& q) const {
  Isometry3d pose;
  if (floating_body_) {
    const auto& joint = floating_body_->getJoint();
    int q_start_index = floating_body_->get_position_start_index();
    auto q_body = q.middleRows(q_start_index, joint.getNumPositions());
    pose = joint.jointTransform(q_body);
  } else {
    pose.setIdentity();
  }
  return pose;
}

const RigidBody* RobotStateTranslator::GetFloatingBody(
    const RigidBodyTree& tree) {
  int floating_joint_count = 0;
  const RigidBody* ret = nullptr;
  for (const auto& body_ptr : tree.bodies) {
    if (body_ptr->hasParent() && body_ptr->getJoint().isFloating()) {
      floating_joint_count++;
      ret = body_ptr.get();
    }
  }

  // Make sure that there's not more than one floating joint.
  if (floating_joint_count > 1) {
    DRAKE_ABORT_MSG(
        "Can't handle robots with more than one floating joint due to the "
        "constraints imposed by the robot_state_t LCM type.");
  }

  // Make sure the floating joint positions and velocities are at the start
  // of the position and velocity vectors.
  if (ret) {
    if (ret->get_position_start_index() != 0 ||
        ret->get_velocity_start_index() != 0) {
      DRAKE_ABORT_MSG(
          "Floating joint positions and velocities must be at the head of the "
          "position and velocity vectors.");
    }
  }
  return ret;
}

}  // drake
