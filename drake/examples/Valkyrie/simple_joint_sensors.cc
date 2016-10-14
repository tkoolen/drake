#include "simple_joint_sensors.h"

using std::unique_ptr;
using std::make_unique;
using std::move;
using std::vector;
using std::map;

using drake::systems::kAbstractValued;
using drake::systems::kVectorValued;
using drake::systems::kContinuousSampling;

namespace drake {
namespace systems {

template <typename T>
SimpleJointSensors<T>::SimpleJointSensors(
    const RigidBodyTree& tree, const vector<const DrakeJoint*>& joints)
    : tree_(tree),
      bodies_(FindSuccessorBodies(tree, joints)),
      kinematics_results_port_(
          DeclareInputPort(kAbstractValued, 1, kContinuousSampling)),
      position_port_(DeclareOutputPort(kAbstractValued, 1, kContinuousSampling)),
      velocity_port_(
          DeclareOutputPort(kAbstractValued, 1, kContinuousSampling)) {
  set_name("SimpleJointSensors");
}

template <typename T>
void SimpleJointSensors<T>::EvalOutput(
    const systems::Context<T>& context,
    systems::SystemOutput<T>* output) const {
  using systems::AbstractValue;
  using systems::KinematicsResults;
//  using SimpleJointSensors<T>::PositionMap;
//  using SimpleJointSensors<T>::VelocityMap;


  const AbstractValue* input =
      EvalAbstractInput(context, kinematics_results_port_.get_index());
  const auto& results = input->GetValue<KinematicsResults<T>>();
  AbstractValue* position_port_data = output->GetMutableData(position_port_.get_index());
  PositionMap& positions = position_port_data->GetMutableValue<PositionMap>();
  AbstractValue *velocity_port_data =
      output->GetMutableData(velocity_port_.get_index());
  VelocityMap& velocities = velocity_port_data->GetMutableValue<VelocityMap>();

  positions.clear();
  velocities.clear();
  for (const auto& body_ptr : bodies_) {
    positions[&body_ptr->getJoint()] = results.get_joint_position(*body_ptr);
    velocities[&body_ptr->getJoint()] = results.get_joint_velocity(*body_ptr);
  }
}

template <typename T>
std::unique_ptr<systems::SystemOutput<T>> SimpleJointSensors<T>::AllocateOutput(
    const systems::Context<T>& context) const {
  using systems::LeafSystemOutput;

  auto output = make_unique<LeafSystemOutput<T>>();

  PositionMap positions;
  VelocityMap velocities;
  for (const auto& body_ptr : bodies_) {
    const auto& joint = body_ptr->getJoint();
    positions.insert({&joint, VectorX<T>(joint.get_num_positions())});
    velocities.insert({&joint, VectorX<T>(joint.get_num_velocities())});
  }

  // TODO(tkoolen): this currently only works because the order of the port
  // declarations in the constructor matches the order in which add_port is
  // called...

  output->add_port(make_unique<Value<PositionMap>>(positions));
  output->add_port(make_unique<Value<VelocityMap>>(velocities));

  return std::unique_ptr<SystemOutput<T>>(output.release());
}

template <typename T>
vector<const RigidBody*> SimpleJointSensors<T>::FindSuccessorBodies(
    const RigidBodyTree& tree, vector<const DrakeJoint*> joints) {
  vector<const RigidBody*> ret;
  for (const auto& joint_ptr : joints) {
    for (const auto& body_ptr : tree.bodies) {
      if (&body_ptr->getJoint() == joint_ptr) {
        ret.push_back(body_ptr.get());
        break;
      }
      throw std::runtime_error("Joint not found.");
    }
  }
  return ret;
}

template <typename T>
const SystemPortDescriptor<T>&
SimpleJointSensors<T>::get_kinematics_results_port() const {
  return kinematics_results_port_;
}

template <typename T>
const SystemPortDescriptor<T>& SimpleJointSensors<T>::get_position_port()
    const {
  return position_port_;
}

template <typename T>
const SystemPortDescriptor<T>& SimpleJointSensors<T>::get_velocity_port()
    const {
  return velocity_port_;
}

// Explicit instantiations.
template class DRAKE_EXPORT SimpleJointSensors<double>;

}  // systems
}  // drake
