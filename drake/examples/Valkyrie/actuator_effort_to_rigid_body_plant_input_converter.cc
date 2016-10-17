#include "actuator_effort_to_rigid_body_plant_input_converter.h"

namespace drake {
namespace systems {

template <typename T>
ActuatorEffortToRigidBodyPlantInputConverter<T>::
    ActuatorEffortToRigidBodyPlantInputConverter(
        const std::vector<const RigidBodyActuator*>& ordered_actuators)
    : ordered_actuators_(ordered_actuators),
      effort_ports_(DeclareEffortInputPorts(ordered_actuators)),
      rigid_body_plant_input_port_(DeclareOutputPort(
          kVectorValued, static_cast<int>(ordered_actuators.size()),
          kContinuousSampling)) {
  set_name("ActuatorEffortToRigidBodyPlantInputConverter");
}

template <typename T>
void ActuatorEffortToRigidBodyPlantInputConverter<T>::EvalOutput(
    const Context<T>& context, SystemOutput<T>* output) const {
  int index = 0;
  for (const auto& actuator : ordered_actuators_) {
    int port_index = effort_ports_.at(actuator)->get_index();
    T effort = EvalVectorInput(context, port_index)->get_value()[0];
    output->GetMutableVectorData(rigid_body_plant_input_port_.get_index())
        ->SetAtIndex(index++, effort);
  }
}

template <typename T>
const std::map<const RigidBodyActuator*, const SystemPortDescriptor<T>*>
ActuatorEffortToRigidBodyPlantInputConverter<T>::DeclareEffortInputPorts(
    const std::vector<const RigidBodyActuator*>& ordered_actuators) {
  std::map<const RigidBodyActuator*, const SystemPortDescriptor<T>*> ret;

  // Currently, all RigidBodyActuators are assumed to be one-dimensional.
  const int effort_length = 1;
  for (const auto& actuator : ordered_actuators) {
    ret[actuator] =
        &DeclareInputPort(kVectorValued, effort_length, kContinuousSampling);
  }
  return ret;
}

template <typename T>
const SystemPortDescriptor<double>&
ActuatorEffortToRigidBodyPlantInputConverter<T>::get_effort_input_port(
    const RigidBodyActuator& actuator) {
  return *effort_ports_.at(&actuator);
}

// Explicit instantiations.
template class DRAKE_EXPORT
    ActuatorEffortToRigidBodyPlantInputConverter<double>;

}  // systems
}  // drake
