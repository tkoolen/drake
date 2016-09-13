#pragma once

#include "drake/drakeRobotCommandToDesiredEffortConverter_export.h"
#include "drake/systems/framework/context.h"
#include "drake/systems/framework/leaf_system.h"

// TODO(tkoolen): namespace
// TODO(tkoolen): export
namespace drake {

class DRAKEROBOTCOMMANDTODESIREDEFFORTCONVERTER_EXPORT
    RobotCommandToDesiredEffortConverter : public systems::LeafSystem<double> {
 public:
  RobotCommandToDesiredEffortConverter();

  ~RobotCommandToDesiredEffortConverter() override;

  // Disable copy and assign.
  RobotCommandToDesiredEffortConverter(
      const RobotCommandToDesiredEffortConverter&) = delete;

  RobotCommandToDesiredEffortConverter& operator=(
      const RobotCommandToDesiredEffortConverter&) = delete;

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
};

}  // drake
