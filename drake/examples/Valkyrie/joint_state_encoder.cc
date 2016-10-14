#include "joint_state_encoder.h"
#include "lcmtypes/bot_core/joint_state_t.hpp"

using std::unique_ptr;
using std::make_unique;
using std::move;

using bot_core::joint_state_t;

namespace drake {
namespace systems {

JointStateEncoder::JointStateEncoder()
    : position_port_(
          DeclareInputPort(kAbstractValued, 1, kContinuousSampling)),
      velocity_port_(
          DeclareInputPort(kAbstractValued, 1, kContinuousSampling)),
      effort_port_(DeclareInputPort(kAbstractValued, 1, kContinuousSampling)),
      message_port_(
          DeclareOutputPort(kAbstractValued, 1, kContinuousSampling)) {}

void JointStateEncoder::EvalOutput(const Context<double>& context,
                                   SystemOutput<double>* output) const {
  const auto& positions =
      EvalAbstractInput(context, position_port_.get_index())
          ->GetValue<PositionMap>();
  const auto& velocities =
      EvalAbstractInput(context, velocity_port_.get_index())
          ->GetValue<VelocityMap>();
  const auto& efforts = EvalAbstractInput(context, effort_port_.get_index())
                            ->GetValue<EffortMap>();
  auto& message = output->GetMutableData(message_port_.get_index())
                      ->GetMutableValue<joint_state_t>();

  // TODO(tkoolen): assumes that positions, velocities, and efforts associated
  // with a joint are scalar-valued.
  for (const auto& actuator_and_effort : efforts) {
    const RigidBodyActuator& actuator = *actuator_and_effort.first;
    const auto& joint = actuator.body_->getJoint();

    message.joint_name.push_back(joint.get_name());
    message.joint_effort.push_back(
        static_cast<float>(actuator_and_effort.second));
    message.joint_position.push_back(
        static_cast<float>(positions.at(&joint)[0]));
    message.joint_velocity.push_back(
        static_cast<float>(velocities.at(&joint)[0]));
  }
}

std::unique_ptr<SystemOutput<double>> JointStateEncoder::AllocateOutput(
    const Context<double>& context) const {
  using systems::LeafSystemOutput;

  auto output = make_unique<LeafSystemOutput<double>>();
  auto message = make_unique<Value<joint_state_t>>(joint_state_t());

  output->add_port(move(message));
  return std::unique_ptr<SystemOutput<double>>(output.release());
}

const SystemPortDescriptor<double>& JointStateEncoder::get_position_port()
    const {
  return position_port_;
}

const SystemPortDescriptor<double>& JointStateEncoder::get_velocity_port()
    const {
  return velocity_port_;
}

const SystemPortDescriptor<double>& JointStateEncoder::get_effort_port()
    const {
  return effort_port_;
}

const SystemPortDescriptor<double>& JointStateEncoder::get_lcm_message_port()
    const {
  return message_port_;
}

}  // systems
}  // drake
