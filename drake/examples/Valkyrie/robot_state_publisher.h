#pragma once

#include <lcm/lcm-cpp.hpp>

#include "drake/drakeRobotStatePublisher_export.h"
#include "drake/systems/framework/leaf_system.h"
#include "drake/systems/framework/context.h"
#include "drake/systems/plants/RigidbodyTree.h"
#include "lcmtypes/bot_core/robot_state_t.hpp"

namespace drake {
namespace examples { // TODO: right namespace?

class DRAKEROBOTSTATEPUBLISHER_EXPORT RobotStatePublisher : public systems::LeafSystem<double> {

  /**
   * A constructor.
   *
   * @param[in] channel The LCM channel on which to publish.
   *
   * @param[in] translator_dictionary A dictionary for obtaining the appropriate
   * translator for a particular LCM channel.
   *
   * @param[in] lcm A pointer to the LCM subsystem.
   */
  RobotStatePublisher(const std::string& channel, ::lcm::LCM* lcm);


  ~RobotStatePublisher() override;

  // Disable copy and assign.
  RobotStatePublisher(const RobotStatePublisher&) = delete;

  RobotStatePublisher& operator=(const RobotStatePublisher&) = delete;

  std::string get_name() const override;

  /**
   * Takes the VectorBase from the input port of the context and publishes
   * it onto an LCM channel.
   */
  void DoPublish(const systems::ContextBase<double>& context) const override;

  /**
   * This System has no output ports so EvalOutput() does nothing.
   */
  void EvalOutput(const systems::ContextBase<double>& context,
                  systems::SystemOutput<double>* output) const override {}

 private:
  // The channel on which to publish LCM messages.
  const std::string channel_;

  // A pointer to the LCM subsystem.
  ::lcm::LCM* lcm_;

  // LCM Message.
  mutable bot_core::robot_state_t message_;

  // Data to send.
  mutable std::vector<uint8_t> message_bytes_;

  // TODO: ports
};

} // examples
} // drake
