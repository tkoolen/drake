#pragma once

#include <map>

#include "drake/systems/plants/RigidBodyTree.h"
#include "drake/systems/plants/rigid_body_plant/kinematics_results.h"
#include "drake/common/drake_export.h"
#include "drake/systems/framework/leaf_system.h"

namespace drake {
namespace systems {

template <typename T>
class DRAKE_EXPORT SimpleJointSensors final : public systems::LeafSystem<T> {
 public:
  using PositionMap = std::map<const DrakeJoint*, VectorX<T>>;
  using VelocityMap = std::map<const DrakeJoint*, VectorX<T>>;

  using System<T>::DeclareInputPort;
  using System<T>::DeclareOutputPort;
  using System<T>::set_name;
  using System<T>::EvalAbstractInput;

  SimpleJointSensors(const RigidBodyTree& tree,
                     const std::vector<const DrakeJoint*>& joints);

  ~SimpleJointSensors() override {}

  // Disable copy and assign.
  SimpleJointSensors(const SimpleJointSensors&) = delete;

  SimpleJointSensors& operator=(const SimpleJointSensors&) = delete;

  void EvalOutput(const Context<T>& context,
                  SystemOutput<T>* output) const override;

  std::unique_ptr<SystemOutput<T>> AllocateOutput(
      const Context<T>& context) const override;

  const SystemPortDescriptor<T>& get_kinematics_results_port() const;

  const SystemPortDescriptor<T>& get_position_port() const;

  const SystemPortDescriptor<T>& get_velocity_port() const;

 private:
  const RigidBodyTree& tree_;

  // TODO(tkoolen): turn this into std::vector<const DrakeJoint*> joints_;
  std::vector<const RigidBody*> bodies_;

  // Input Port.
  const SystemPortDescriptor<T>& kinematics_results_port_;

  // Output ports.
  const SystemPortDescriptor<T>& position_port_;
  const SystemPortDescriptor<T>& velocity_port_;

  static std::vector<const RigidBody*> FindSuccessorBodies(
      const RigidBodyTree& tree, std::vector<const DrakeJoint*> joints);
};

}  // systems
}  // drake
