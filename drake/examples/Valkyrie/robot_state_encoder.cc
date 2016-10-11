#include "drake/examples/Valkyrie/robot_state_encoder.h"
#include "drake/common/constants.h"
#include "drake/systems/framework/system_port_descriptor.h"
#include "drake/systems/plants/rigid_body_plant/kinematics_results.h"
#include "drake/util/drakeUtil.h"
#include "drake/util/lcmUtil.h"

using std::make_unique;
using std::move;
using std::unique_ptr;

using bot_core::robot_state_t;

using Eigen::Dynamic;
using Eigen::Index;
using Eigen::Matrix;
using Eigen::Isometry3d;
using Eigen::Translation3d;

namespace drake {
namespace systems {

RobotStateEncoder::RobotStateEncoder(const RigidBodyTree& tree)
    : tree_(CheckPreConditions(tree)),
      floating_body_(tree.bodies[1]->getJoint().is_floating()
                         ? tree.bodies[1].get()
                         : nullptr),
      lcm_message_port_(
          DeclareOutputPort(kAbstractValued, 1, kContinuousSampling)),
      kinematics_results_port_(
          DeclareInputPort(kAbstractValued, 1, kContinuousSampling)),
      effort_ports_(DeclareEffortInputPorts()),
      foot_wrench_ports_(DeclareWrenchInputPorts()),
      hand_wrench_ports_(DeclareWrenchInputPorts()) {
  set_name("RobotStateEncoder");
}

RobotStateEncoder::~RobotStateEncoder() {}

void RobotStateEncoder::EvalOutput(const Context<double>& context,
                                   SystemOutput<double>* output) const {
  auto& message = output->GetMutableData(lcm_message_port_.get_index())
                      ->GetMutableValue<robot_state_t>();
  message.utime = static_cast<int64_t>(context.get_time() * 1e6);
  SetStateAndEfforts(message, context);
  SetForceTorque(message, context);
}

std::unique_ptr<SystemOutput<double>> RobotStateEncoder::AllocateOutput(
    const Context<double>& context) const {
  auto output = make_unique<LeafSystemOutput<double>>();

  auto data = make_unique<Value<robot_state_t>>(robot_state_t());
  output->add_port(move(data));

  return std::unique_ptr<SystemOutput<double>>(output.release());
}

const SystemPortDescriptor<double>& RobotStateEncoder::get_lcm_message_port()
    const {
  return lcm_message_port_;
}

const SystemPortDescriptor<double>& RobotStateEncoder::get_kinematics_results_port() const {
  return kinematics_results_port_;
}

const SystemPortDescriptor<double>& RobotStateEncoder::get_effort_port(
    const RigidBodyActuator& actuator) const {
  return *effort_ports_.at(&actuator);
}

const SystemPortDescriptor<double>&
RobotStateEncoder::get_foot_contact_wrench_port(const Side& side) const {
  return *foot_wrench_ports_.at(side);
}

const SystemPortDescriptor<double>&
RobotStateEncoder::get_hand_contact_wrench_port(const Side& side) const {
  return *hand_wrench_ports_.at(side);
}

std::map<const RigidBodyActuator*, const SystemPortDescriptor<double>*>
RobotStateEncoder::DeclareEffortInputPorts() {
  std::map<const RigidBodyActuator*, const SystemPortDescriptor<double>*> ret;

  // Currently, all RigidBodyActuators are assumed to be one-dimensional.
  const int actuator_effort_length = 1;
  for (const auto& actuator : tree_.actuators) {
    ret[&actuator] = &DeclareInputPort(kVectorValued, actuator_effort_length,
                                       kContinuousSampling);
  }
  return ret;
}

std::map<Side, const SystemPortDescriptor<double>*>
RobotStateEncoder::DeclareWrenchInputPorts() {
  std::map<Side, const SystemPortDescriptor<double>*> ret;
  for (const auto& side : Side::values) {
    ret[side] =
        &DeclareInputPort(kVectorValued, kTwistSize, kContinuousSampling);
  }
  return ret;
}

void RobotStateEncoder::SetStateAndEfforts(
    robot_state_t& message, const Context<double>& context) const {
  const auto& kinematics_results =
      EvalAbstractInput(context, kinematics_results_port_.get_index())
          ->GetValue<KinematicsResults<double>>();

  // Pose of floating body with respect to world.
  Isometry3d floating_body_to_world =
      kinematics_results.get_pose_in_world(*floating_body_);
  EncodePose(floating_body_to_world, message.pose);

  // Twist of floating body with respect to world.
  TwistVector<double> floating_body_twist_in_world =
      kinematics_results.get_twist_with_respect_to_world(*floating_body_);
  // To match usage of robot_state_t throughout OpenHumanoids code, transform
  // twist in body frame to a frame that has the same orientation as world
  // frame, but the same origin as the floating body frame.
  Isometry3d world_to_world_aligned_body(
      Translation3d(-floating_body_to_world.translation()));
  TwistVector<double> floating_body_twist_in_in_world_aligned_body =
      transformSpatialMotion(world_to_world_aligned_body,
                             floating_body_twist_in_world);
  EncodeTwist(floating_body_twist_in_in_world_aligned_body, message.twist);

  // Joint names, positions, velocities, and efforts.
  // Note: the order of the actuators in the rigid body tree determines the
  // order of the joint_name, joint_position, joint_velocity, and
  // joint_effort fields.
  message.joint_name.clear();
  message.joint_position.clear();
  message.joint_velocity.clear();
  message.joint_effort.clear();
  for (const auto& actuator : tree_.actuators) {
    const auto& body = *actuator.body_;
    int effort_port_index = effort_ports_.at(&actuator)->get_index();

    // To match usage of robot_state_t throughout OpenHumanoids code, set
    // joint_names field to position coordinate names.
    int position_index = body.get_position_start_index();
    message.joint_name.push_back(tree_.get_position_name(position_index));

    auto position =
        static_cast<float>(kinematics_results.get_joint_position(body)[0]);
    auto velocity =
        static_cast<float>(kinematics_results.get_joint_velocity(body)[0]);
    auto effort = static_cast<float>(
        EvalVectorInput(context, effort_port_index)->GetAtIndex(0));

    message.joint_position.push_back(position);
    message.joint_velocity.push_back(velocity);
    message.joint_effort.push_back(effort);
  }
  message.num_joints = static_cast<int16_t>(message.joint_name.size());
}

void RobotStateEncoder::SetForceTorque(robot_state_t& message,
                                       const Context<double>& context) const {
  auto& force_torque = message.force_torque;
  const int kTorqueXIndex = 0;
  const int kTorqueYIndex = 1;
  const int kForceZIndex = 5;
  {
    auto left_foot_wrench =
        EvalVectorInput(context, foot_wrench_ports_.at(Side::LEFT)->get_index())
            ->get_value();
    force_torque.l_foot_force_z =
        static_cast<float>(left_foot_wrench[kForceZIndex]);
    force_torque.l_foot_torque_x =
        static_cast<float>(left_foot_wrench[kTorqueXIndex]);
    force_torque.l_foot_torque_y =
        static_cast<float>(left_foot_wrench[kTorqueYIndex]);
  }
  {
    auto right_foot_wrench =
        EvalVectorInput(context,
                        foot_wrench_ports_.at(Side::RIGHT)->get_index())
            ->get_value();
    force_torque.r_foot_force_z =
        static_cast<float>(right_foot_wrench[kForceZIndex]);
    force_torque.r_foot_torque_x =
        static_cast<float>(right_foot_wrench[kTorqueXIndex]);
    force_torque.r_foot_torque_y =
        static_cast<float>(right_foot_wrench[kTorqueYIndex]);
  }
  {
    auto left_hand_wrench =
        EvalVectorInput(context, hand_wrench_ports_.at(Side::LEFT)->get_index())
            ->get_value();
    eigenVectorToCArray(left_hand_wrench.head<kSpaceDimension>(),
                        force_torque.l_hand_torque);
    eigenVectorToCArray(left_hand_wrench.tail<kSpaceDimension>(),
                        force_torque.l_hand_torque);
  }
  {
    auto right_hand_wrench =
        EvalVectorInput(context,
                        hand_wrench_ports_.at(Side::RIGHT)->get_index())
            ->get_value();
    eigenVectorToCArray(right_hand_wrench.head<kSpaceDimension>(),
                        force_torque.r_hand_torque);
    eigenVectorToCArray(right_hand_wrench.tail<kSpaceDimension>(),
                        force_torque.r_hand_torque);
  }
}

int RobotStateEncoder::num_floating_joint_positions() const {
  int num_floating_joint_positions =
      floating_body_ ? floating_body_->getJoint().get_num_positions() : 0;
  return num_floating_joint_positions;
}

int RobotStateEncoder::num_floating_joint_velocities() const {
  int num_floating_joint_velocities =
      floating_body_ ? floating_body_->getJoint().get_num_velocities() : 0;
  return num_floating_joint_velocities;
}

Eigen::Isometry3d RobotStateEncoder::EvalFloatingBodyPose(
    const Eigen::Ref<const Eigen::VectorXd>& q) const {
  Isometry3d pose;
  if (floating_body_) {
    const auto& joint = floating_body_->getJoint();
    int q_start_index = floating_body_->get_position_start_index();
    auto q_body = q.middleRows(q_start_index, joint.get_num_positions());
    pose = joint.jointTransform(q_body);
  } else {
    pose.setIdentity();
  }
  return pose;
}

TwistVector<double> RobotStateEncoder::EvalFloatingBodyTwistInBodyFrame(
    const Eigen::Ref<const Eigen::VectorXd>& q,
    const Eigen::Ref<const Eigen::VectorXd>& v) const {
  TwistVector<double> twist;
  if (floating_body_) {
    const auto& joint = floating_body_->getJoint();

    int q_start_index = floating_body_->get_position_start_index();
    auto q_body = q.middleRows(q_start_index, joint.get_num_positions());
    int v_start_index = floating_body_->get_velocity_start_index();
    auto v_body = v.middleRows(v_start_index, joint.get_num_velocities());

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
RobotStateEncoder::TransformTwistFromBodyFrameToWorldAlignedBodyFrame(
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

const RigidBodyTree& RobotStateEncoder::CheckPreConditions(
    const RigidBodyTree& tree) {
  if (tree.get_num_bodies() < 2) {
    DRAKE_ABORT_MSG("This class assumes at least one non-world body.");
  }

  bool floating_joint_found = false;
  for (const auto& body_ptr : tree.bodies) {
    if (body_ptr->has_parent_body()) {
      const auto& joint = body_ptr->getJoint();
      if (joint.is_floating()) {
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
        if (joint.get_num_positions() > 1 || joint.get_num_velocities() > 1) {
          DRAKE_ABORT_MSG(
              "robot_state_t assumes non-floating joints to be "
              "1-DoF or fixed.");
        }
      }
    }
  }

  return tree;
}

}  // systems
}  // drake
