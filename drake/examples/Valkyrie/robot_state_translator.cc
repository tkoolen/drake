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
      tree_(CheckPreConditions(tree)),
      floating_body_(tree.bodies[1]->getJoint().isFloating() ? tree.bodies[1]
          .get() : nullptr) {
  InitializeMessage();
}
void RobotStateTranslator::InitializeMessage() const {
  // To match usage of robot_state_t throughout OpenHumanoids code, set
  // joint_names field to position coordinate names.
  int non_floating_joint_position_start_index = num_floating_joint_positions();
  for (int i = non_floating_joint_position_start_index; i < tree_
      .number_of_positions(); i++) {
    message_.joint_name.push_back(tree_.getPositionName(i));
  }
  message_.num_joints = static_cast<int16_t>(message_.joint_name.size());
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

  Isometry3d floating_body_to_world = EvalFloatingBodyPose(q);
  EncodePose(floating_body_to_world, message_.pose);

  TwistVector<double> twist_in_body = EvalFloatingBodyTwistInBodyFrame(q, v);

  // To match usage of robot_state_t throughout OpenHumanoids code, transform
  // twist in body frame to a frame that has the same orientation as world
  // frame, but the same origin as the floating body frame.
  TwistVector<double> twist_in_world_aligned_body =
      TransformTwistFromBodyFrameToWorldAlignedBodyFrame(floating_body_to_world,
                                                         twist_in_body);
  EncodeTwist(twist_in_world_aligned_body, message_.twist);

  // set the joint_position, joint_velocity, joint_effort parts
  auto q_non_floating = q.tail(q.size() - num_floating_joint_positions());
  eigenVectorToStdVector(q_non_floating, message_.joint_position);

  auto v_non_floating = v.tail(v.size() - num_floating_joint_velocities());
  eigenVectorToStdVector(v_non_floating, message_.joint_velocity);

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

TwistVector<double>
RobotStateTranslator::TransformTwistFromBodyFrameToWorldAlignedBodyFrame(
    Isometry3d &floating_body_to_world,
    const TwistVector<double> &twist_in_body) const {
  Isometry3d floating_body_to_world_aligned_floating_body;
  floating_body_to_world_aligned_floating_body.linear() = floating_body_to_world
      .linear();
  floating_body_to_world_aligned_floating_body.translation().setZero();
  floating_body_to_world_aligned_floating_body.makeAffine();
  TwistVector<double> twist_in_world_aligned_body = transformSpatialMotion
      (floating_body_to_world_aligned_floating_body, twist_in_body);
  return twist_in_world_aligned_body;
}

const RigidBodyTree& RobotStateTranslator::CheckPreConditions(
    const RigidBodyTree& tree) {
  if (tree.get_number_of_bodies() < 2) {
    DRAKE_ABORT_MSG("This class assumes at least one non-world body.");
  }

  bool floating_joint_found = false;
  for (const auto& body_ptr : tree.bodies) {
    if (body_ptr->hasParent()) {
      const auto& joint = body_ptr->getJoint();
      if (joint.isFloating()) {
        if (floating_joint_found) {
          DRAKE_ABORT_MSG("robot_state_t assumes at most one floating joint.");
        }
        floating_joint_found = true;

        if (body_ptr != tree.bodies[1]) {
          DRAKE_ABORT_MSG(
              "This class assumes that the first non-world body is the "
              "floating body.");
        }

        if (body_ptr->get_position_start_index() != 0 ||
            body_ptr->get_velocity_start_index() != 0) {
          DRAKE_ABORT_MSG(
              "This class assumes that floating joint positions and are at the "
              "head of the position and velocity vectors.");
        }
      } else {
        if (joint.getNumPositions() > 1 || joint.getNumVelocities() > 1) {
          DRAKE_ABORT_MSG(
              "robot_state_t assumes non-floating joints to be "
              "1-DoF or fixed.");
        }
      }
    }
  }

  return tree;
}

}  // drake
