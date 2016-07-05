#pragma once

#include <stdexcept>
#include <lcm/lcm-cpp.hpp>
#include <Eigen/Dense>
#include "drake/systems/System.h"
#include "drake/systems/plants/RigidBodyTree.h"
#include "lcmtypes/bot_core/robot_state_t.hpp"
#include "drake/util/drakeUtil.h"
#include "drake/util/lcmUtil.h"

/*
 * NOTE(tkoolen): the bot_core::robot_state_t is somewhat strange.
 * My complaints:
 * * the joint_name field is a misnomer; historically this has been set to the
 *   joint position names
 * * it's not necessary to publish all the joint names in every message
 * * the joint_effort and force_torque fields should not be part of the message
 *   in my opinion
 * * pose and twist should not be separate, but should be combined with
 *   joint_position and joint_velocity respectively, into configuration and
 *   velocity
 * * number of positions is assumed equal to number of velocities
 */

template<template<typename> class RobotStateVector>
class RobotStatePublisher {
 public:
  template<typename ScalarType>
  using StateVector = NullVector<ScalarType>;
  template<typename ScalarType>
  using OutputVector = RobotStateVector<ScalarType>;
  template<typename ScalarType>
  using InputVector = RobotStateVector<ScalarType>;

  RobotStatePublisher(std::shared_ptr<lcm::LCM> lcm,
                      const std::string &channel,
                      const std::string &urdf_filename,
                      const DrakeJoint::FloatingBaseType floating_base_type)
      : lcm_(lcm), channel_(channel),
        tree_(new RigidBodyTree(urdf_filename, floating_base_type)),
        floating_body_(GetFloatingBody(tree_)) {

    // compute number of non-floating joints
    int num_non_floating_joints = 0;
    for (const auto &body_ptr : tree_->bodies) {
      if (body_ptr->hasParent()) {
        const auto &joint = body_ptr->getJoint();
        if (!joint.isFloating()) {
          num_non_floating_joints++;

          // make sure all non-floating joints are 1-DoF
          if (joint.getNumPositions() != 1 ||
              joint.getNumVelocities() != 1) {
            std::cerr
                << "Can't handle robots for non-floating joints are not all "
                    "1-DoF due to constraints imposed by the robot_state_t LCM "
                    "type."
                << std::endl;
            std::abort(-1);
          }
        }
      }
    }

    // set num_joints, resize vectors
    msg_.num_joints = num_non_floating_joints;
    msg_.joint_name.reserve(num_non_floating_joints);
    msg_.joint_position.resize(num_non_floating_joints);
    msg_.joint_velocity.resize(num_non_floating_joints);
    msg_.joint_effort.resize(num_non_floating_joints);

    // set joint_names
    int position_index = floating_body_ ? floating_body_->position_num_start +
        floating_body_->getJoint().getNumPositions() : 0;
    for (; position_index < tree_->getNumPositions(); position_index++) {
      msg_.joint_position.push_back(tree_->getPositionName(position_index));
    }
  }

  StateVector<double> dynamics(const double &t, const StateVector<double> &x,
                               const InputVector<double> &u) const {
    return StateVector<double>();
  }

  OutputVector<double> output(const double &t, const StateVector<double> &x,
                              const InputVector<double> &u) const {
    using namespace Eigen;

    msg.timestamp = static_cast<int64_t>(t * 1000.0);

    auto uvec = toEigen(u);
    auto q = uvec.head(tree_->number_of_positions());
    auto v = uvec.tail(tree_->number_of_velocities());

    // set the pose and twist parts
    Isometry3d pose;
    Matrix<double, 6, 1> twist;
    if (floating_body_) {
      auto &joint = floating_body_->getJoint();

      int q_start_index = floating_body_->position_num_start;
      auto q_body = q.middleRows(q_start_index, joint.getNumPositions());
      int v_start_index = floating_body_->velocity_num_start;
      auto v_body = v.middleRows(v_start_index, joint.getNumVelocities());

      pose = floating_joint_->jointTransform(q_body);
      Eigen::Matrix<Scalar, TWIST_SIZE, Eigen::Dynamic, 0, TWIST_SIZE,
                    DrakeJoint::MAX_NUM_VELOCITIES> motion_subspace;
      Matrix<Scalar, Dynamic, Dynamic> *dmotion_subspacedq = nullptr;
      joint.motionSubspace(q_body, motion_subspace, dmotion_subspacedq);
      twist = motion_subspace * v_body;
    }
    else {
      pose.setIdentity();
      twist.setZero();
    }
    encodePosition3d(pose, msg_.pose);
    encodeTwist(pose, msg_.twist);

    // set the joint_position, joint_velocity, joint_effort parts
    int num_floating_joint_positions =
        floating_body_? floating_body_->getJoint().getNumPositions() : 0;
    int num_floating_joint_velocities =
        floating_body_? floating_body_->getJoint().getNumVelocities() : 0;

    auto qJoint = q.tail(tree_->number_of_positions() - )
    eigenVectorToStdVector()

    lcm->publish(channel_, &msg_);


    // TODO: set efforts
    // TODO: set force_torque

    return u;  // pass the output through
  }

 private:
  static const RigidBody *GetFloatingBody(const RigidBodyTree &tree) {
    int floating_joint_count = 0;
    const RigidBody *ret = nullptr;
    for (const auto &body_ptr : tree.bodies) {
      if (body_ptr->hasParent() && body_ptr->getJoint().isFloating()) {
        floating_joint_count++;
        ret = body_ptr.get();
      }
    }

    // make sure that there's not more than one floating joint
    if (floating_joint_count > 1) {
      std::cerr
          << "Can't handle robots with more than one floating joint due to the "
              "constraints imposed by the robot_state_t LCM type."
          << std::endl;
      std::abort(-1);
    }

    // make sure the floating joint positions and velocities are at the start
    // of the position and velocity vectors
    if (ret) {
      if (ret->position_num_start != 0 || ret->velocity_num_start != 0) {
        std::cerr
            << "Floating joint positions and velocities must be at the head of "
                "the position and velocity vectors."
            << std::endl;
        std::abort(-1);
      }
    }

    return ret;
  }

  const std::shared_ptr<lcm::LCM> lcm_;
  const std::string channel_;
  mutable std::shared_ptr<RigidBodyTree> tree_;
  const RigidBody *const floating_body_;
  mutable bot_core::robot_state_t msg_;
};

