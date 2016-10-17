#pragma once

#include <vector>

#include "drake/common/drake_export.h"
#include "drake/systems/plants/rigid_body_actuator.h"
#include "drake/systems/framework/leaf_system.h"

namespace drake {
namespace systems {

template <typename T>
class DRAKE_EXPORT ActuatorEffortToRigidBodyPlantInputConverter : public LeafSystem<T> {
 public:
  using System<T>::DeclareInputPort;
  using System<T>::DeclareOutputPort;
  using System<T>::set_name;
  using System<T>::EvalVectorInput;

  ActuatorEffortToRigidBodyPlantInputConverter(
      const std::vector<const RigidBodyActuator*>& ordered_actuators);

  ~ActuatorEffortToRigidBodyPlantInputConverter() override {}

  // Disable copy and assign.
  ActuatorEffortToRigidBodyPlantInputConverter(
      const ActuatorEffortToRigidBodyPlantInputConverter&) = delete;

  ActuatorEffortToRigidBodyPlantInputConverter& operator=(
      const ActuatorEffortToRigidBodyPlantInputConverter&) = delete;

  void EvalOutput(const Context<T>& context,
                  SystemOutput<T>* output) const override;

  const SystemPortDescriptor<double> &
  get_effort_input_port(const RigidBodyActuator &actuator);
 private:
  std::vector<const RigidBodyActuator*> ordered_actuators_;
  const std::map<const RigidBodyActuator*, const SystemPortDescriptor<T>*>
      effort_ports_;
  const SystemPortDescriptor<T>& rigid_body_plant_input_port_;
  const std::map<const RigidBodyActuator *, const SystemPortDescriptor<T> *>
  DeclareEffortInputPorts(const std::vector<const RigidBodyActuator *> &vector);
};


} // systems
} // drake
