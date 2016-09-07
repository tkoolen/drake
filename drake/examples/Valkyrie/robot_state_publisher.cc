#include "robot_state_publisher.h"
#include "drake/common/constants.h"
#include "drake/util/lcmUtil.h"

using Eigen::Dynamic;
using Eigen::Index;
using Eigen::Matrix;
using Eigen::Isometry3d;
using drake::systems::VectorBase;
using drake::systems::kVectorValued;
using drake::systems::kContinuousSampling;

namespace drake {
namespace examples {

RobotStatePublisher::RobotStatePublisher(const RigidBodyTree& tree,
                                         int num_actuators,
                                         const std::string& channel,
                                         lcm::LCM* lcm)
    : tree_(CheckPreConditions(tree)),
      floating_body_(tree.bodies[1]->getJoint().isFloating()
                         ? tree.bodies[1].get()
                         : nullptr),
      channel_(channel),
      lcm_(lcm) {
  int num_states = tree_.number_of_positions() + tree_.number_of_velocities();

  int input_num = 0;
  DeclareInputPort(kVectorValued, num_states, kContinuousSampling);
  state_port_index_ = input_num++;
  DeclareInputPort(kVectorValued, num_actuators, kContinuousSampling);
  torque_port_index_ = input_num++;
  DeclareInputPort(kVectorValued, kTwistSize, kContinuousSampling);
  foot_wrench_port_indices_[Side::LEFT] = input_num++;
  DeclareInputPort(kVectorValued, kTwistSize, kContinuousSampling);
  foot_wrench_port_indices_[Side::RIGHT] = input_num;

  InitializeMessage();
}

RobotStatePublisher::~RobotStatePublisher() {}

std::string RobotStatePublisher::get_name() const {
  return "RobotStatePublisher::" + channel_;
}

void RobotStatePublisher::DoPublish(
    const systems::ContextBase<double>& context) const {
  const VectorBase<double>* const state_vector =
      context.get_vector_input(state_port_index_);

  // TODO: message.utime

  auto x = state_vector->get_value();  // TODO: need to check for nullness?
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

  // TODO: extract out into its own method
  const int lcm_message_length = message_.getEncodedSize();
  message_bytes_.resize(static_cast<size_t>(lcm_message_length));
  message_.encode(message_bytes_.data(), 0, lcm_message_length);
  lcm_->publish(channel_, message_bytes_.data(),
                static_cast<unsigned int>(message_bytes_.size()));
}

void RobotStatePublisher::InitializeMessage() const {
  // To match usage of robot_state_t throughout OpenHumanoids code, set
  // joint_names field to position coordinate names.
  int non_floating_joint_position_start_index = num_floating_joint_positions();
  for (int i = non_floating_joint_position_start_index;
       i < tree_.number_of_positions(); i++) {
    message_.joint_name.push_back(tree_.getPositionName(i));
  }
  message_.num_joints = static_cast<int16_t>(message_.joint_name.size());
}

int RobotStatePublisher::num_floating_joint_positions() const {
  int num_floating_joint_positions =
      floating_body_ ? floating_body_->getJoint().getNumPositions() : 0;
  return num_floating_joint_positions;
}

int RobotStatePublisher::num_floating_joint_velocities() const {
  int num_floating_joint_velocities =
      floating_body_ ? floating_body_->getJoint().getNumVelocities() : 0;
  return num_floating_joint_velocities;
}

Eigen::Isometry3d RobotStatePublisher::EvalFloatingBodyPose(
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

TwistVector<double> RobotStatePublisher::EvalFloatingBodyTwistInBodyFrame(
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
RobotStatePublisher::TransformTwistFromBodyFrameToWorldAlignedBodyFrame(
    Isometry3d& floating_body_to_world,
    const TwistVector<double>& twist_in_body) const {
  Isometry3d floating_body_to_world_aligned_floating_body;
  floating_body_to_world_aligned_floating_body.linear() =
      floating_body_to_world.linear();
  floating_body_to_world_aligned_floating_body.translation().setZero();
  floating_body_to_world_aligned_floating_body.makeAffine();
  TwistVector<double> twist_in_world_aligned_body = transformSpatialMotion(
      floating_body_to_world_aligned_floating_body, twist_in_body);
  return twist_in_world_aligned_body;
}

const RigidBodyTree& RobotStatePublisher::CheckPreConditions(
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

}  // examples
}  // drake
