
#include <memory>
#include <drake/systems/framework/primitives/pass_through.h>

#include "drake/common/drake_path.h"
#include "drake/common/text_logging.h"
#include "drake/examples/Valkyrie/robot_state_translator.h"
#include "drake/systems/analysis/simulator.h"
#include "drake/systems/framework/diagram_builder.h"
#include "drake/systems/framework/primitives/constant_vector_source.h"
#include "drake/systems/lcm/lcm_publisher_system.h"
#include "drake/systems/lcm/lcm_subscriber_system.h"
#include "drake/systems/lcm/translator_between_lcmt_drake_signal.h"
#include "drake/systems/plants/parser_urdf.h"
#include "drake/systems/plants/rigid_body_system/rigid_body_plant.h"

namespace drake {
namespace examples {
namespace valkyrie {

int do_main(int argc, const char* argv[]) {
  using std::move;
  using std::make_unique;
  using systems::RigidBodyPlant;
  using systems::lcm::LcmPublisherSystem;
  using systems::lcm::LcmSubscriberSystem;
  using systems::lcm::TranslatorBetweenLcmtDrakeSignal;
  using systems::DiagramContext;
  using Eigen::VectorXd;

  drake::log()->set_level(spdlog::level::trace);

  // Create RigidBodyTree.
  auto tree_ptr = make_unique<RigidBodyTree>();
  drake::parsers::urdf::AddModelInstanceFromUrdfFile(
      drake::GetDrakePath() +
          "/examples/Valkyrie/urdf/urdf/"
          "valkyrie_A_sim_drake_one_neck_dof_wide_ankle_rom.urdf",
      DrakeJoint::ROLLPITCHYAW, nullptr /* weld to frame */, tree_ptr.get());

  // Instantiate a RigidBodyPlant from the RigidBodyTree.
  RigidBodyPlant<double> plant(move(tree_ptr));

  // LCM communication.
  auto lcm = std::make_unique<lcm::LCM>();

  // Create torque command subscriber.
  TranslatorBetweenLcmtDrakeSignal torque_translator(plant.get_input_size());
  auto torque_command_source = make_unique<LcmSubscriberSystem>(
      "ROBOT_TORQUE", torque_translator, lcm.get());
//  VectorXd tau = VectorXd::Zero(plant.get_input_size());
//  torque_command_source = make_unique<drake::systems::ConstantVectorSource<double>(tau);

  // Placeholder for actuator dynamics.
  systems::PassThrough<double> actuators(plant.get_input_size());

  // Create robot_state_t publisher.
  const RobotStateTranslator robot_state_translator(
      plant.get_multibody_world());
  auto robot_state_publisher = make_unique<LcmPublisherSystem>(
      "EST_ROBOT_STATE", robot_state_translator, lcm.get());

  // Wire things up.
  // TODO: use syntactic sugar from cars-simulator2 branch
  auto builder = std::make_unique<systems::DiagramBuilder<double>>();
  builder->Connect(torque_command_source->get_output_port(0), actuators
      .get_input_port(0));
  builder->Connect(actuators.get_output_port(0), plant.get_input_port(0));
  builder->Connect(plant.get_output_port(0),
                   robot_state_publisher->get_input_port(0));
  auto diagram = builder->Build();

  auto lcm_receive_thread =
      std::make_unique<systems::lcm::LcmReceiveThread>(lcm.get());

  // Create simulator.
  auto simulator = std::make_unique<systems::Simulator<double>>(*diagram);
  simulator->Initialize();
  auto context = simulator->get_mutable_context();

  // Set initial state.
  auto plant_context = diagram->GetMutableSubsystemContext(context, &plant);
  plant.SetZeroConfiguration(plant_context);

  while (true) {
    const double time = context->get_time();
    SPDLOG_TRACE(drake::log(), "Time is now {}", time);
    simulator->StepTo(time + 0.01);
  }

  return 0;
}

}  // valkyrie
}  // examples
}  // drake

int main(int argc, const char* argv[]) {
  return drake::examples::valkyrie::do_main(argc, argv);
}