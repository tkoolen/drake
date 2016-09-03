
#include <memory>

#include "drake/systems/analysis/simulator.h"
#include "drake/common/drake_path.h"
#include "drake/systems/framework/diagram_builder.h"
#include "drake/common/text_logging.h"
#include "drake/systems/lcm/lcm_publisher_system.h"
#include "drake/systems/lcm/lcm_subscriber_system.h"
#include "drake/systems/plants/parser_urdf.h"
#include "drake/systems/plants/rigid_body_system/rigid_body_plant.h"

namespace drake {
namespace examples {
namespace valkyrie {

int do_main(int argc, const char* argv[]) {
  using std::move;
  using std::make_unique;
  using drake::systems::RigidBodyPlant;
  using systems::lcm::LcmPublisherSystem;
  using systems::lcm::LcmSubscriberSystem;

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

  // Set up LCM communication.
  auto lcm = std::make_unique<lcm::LCM>();
  const RobotStateTranslator robot_state_translator;
  auto robot_state_publisher = make_unique<LcmPublisherSystem>
      ("EST_ROBOT_STATE", robot_state_translator, lcm.get());

  // Wire things up.
  auto builder = std::make_unique<systems::DiagramBuilder<double>>();
  // TODO: more connections
  builder->Connect(plant.get_output_port(0), robot_state_publisher
      ->get_input_port(0));
  auto diagram = builder->Build();

  auto lcm_receive_thread = std::make_unique<systems::lcm::LcmReceiveThread>(
      lcm.get());

  auto simulator = std::make_unique<systems::Simulator<double>>(*diagram);
  simulator->Initialize();
  while (true) {
    const double time = simulator->get_context().get_time();
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
