#include "drake/util/lcmUtil.h"
#include "drake/util/drakeUtil.h"
#include "drake/util/drakeGeometryUtil.h"

using namespace Eigen;

void encodeVector3d(
    const Eigen::Ref<const Eigen::Vector3d>& vec, bot_core::vector_3d_t& msg) {
  msg.x = vec[0];
  msg.y = vec[1];
  msg.z = vec[2];
}

void encodeQuaternion(
    const Eigen::Ref<const Eigen::Vector4d>& vec, bot_core::quaternion_t& msg) {
  msg.w = vec[0];
  msg.x = vec[0];
  msg.y = vec[0];
  msg.z = vec[0];
}

void encodePose(const Eigen::Isometry3d& pose, bot_core::position_3d_t& msg) {
  auto rotation = rotmat2quat(pose.linear());
  encodeQuaternion(rotation, msg.rotation);
  encodeVector3d(pose.translation(), msg.translation);
}

void encodeTwist(const Eigen::Ref<const Eigen::Matrix<double, 6, 1>>& twist,
    bot_core::twist_t& msg) {

}

void encodePolynomial(const Polynomial<double>& polynomial,
                      drake::lcmt_polynomial& msg) {
  eigenVectorToStdVector(polynomial.getCoefficients(), msg.coefficients);
  msg.num_coefficients = polynomial.getNumberOfCoefficients();
}

Polynomial<double> decodePolynomial(const drake::lcmt_polynomial& msg) {
  Map<const VectorXd> coefficients(msg.coefficients.data(),
                                   msg.coefficients.size());
  return Polynomial<double>(coefficients);
}

void encodePiecewisePolynomial(
    const PiecewisePolynomial<double>& piecewise_polynomial,
    drake::lcmt_piecewise_polynomial& msg) {
  msg.num_segments = piecewise_polynomial.getNumberOfSegments();
  msg.num_breaks = piecewise_polynomial.getNumberOfSegments() + 1;
  msg.breaks = piecewise_polynomial.getSegmentTimes();
  msg.polynomial_matrices.resize(piecewise_polynomial.getNumberOfSegments());
  for (int i = 0; i < piecewise_polynomial.getNumberOfSegments(); ++i) {
    encodePolynomialMatrix<Eigen::Dynamic, Eigen::Dynamic>(
        piecewise_polynomial.getPolynomialMatrix(i),
        msg.polynomial_matrices[i]);
  }
}

PiecewisePolynomial<double> decodePiecewisePolynomial(
    const drake::lcmt_piecewise_polynomial& msg) {
  typedef PiecewisePolynomial<double>::PolynomialMatrix PolynomialMatrix;
  std::vector<PolynomialMatrix> polynomial_matrices;
  for (size_t i = 0; i < msg.polynomial_matrices.size(); ++i) {
    polynomial_matrices.push_back(
        decodePolynomialMatrix<Dynamic, Dynamic>(msg.polynomial_matrices[i]));
  }
  return PiecewisePolynomial<double>(polynomial_matrices, msg.breaks);
}

void verifySubtypeSizes(drake::lcmt_support_data& support_data) {
  // Check for errors in sizes of variable-length fields.
  if (support_data.contact_pts.size() != 3) {
    throw std::runtime_error("contact_pts must have 3 rows");
  }
  for (int j = 0; j < 3; ++j) {
    if (static_cast<int32_t>(support_data.contact_pts[j].size()) !=
        support_data.num_contact_pts) {
      std::stringstream msg;
      msg << "num_contact_pts must match the size of each row of contact_pts."
          << std::endl;
      msg << "num_contact_pts: " << support_data.num_contact_pts
          << ", support_data.contact_pts[" << j
          << "].size(): " << support_data.contact_pts[j].size() << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
}

void verifySubtypeSizes(drake::lcmt_qp_controller_input& qp_input) {
  // Check (and try to fix) errors in the sizes of the variable-length fields in
  // our message
  if (static_cast<int32_t>(qp_input.support_data.size()) !=
      qp_input.num_support_data) {
    std::cerr << "WARNING: num support data doesn't match" << std::endl;
    qp_input.num_support_data = qp_input.support_data.size();
  }
  if (static_cast<int32_t>(qp_input.body_motion_data.size()) !=
      qp_input.num_tracked_bodies) {
    std::cerr << "WARNING: num tracked bodies doesn't match" << std::endl;
    qp_input.num_tracked_bodies = qp_input.body_motion_data.size();
  }
  if (static_cast<int32_t>(qp_input.body_wrench_data.size()) !=
      qp_input.num_external_wrenches) {
    std::cerr << "WARNING: num external wrenches doesn't match" << std::endl;
    qp_input.num_external_wrenches = qp_input.body_wrench_data.size();
  }
  if (static_cast<int32_t>(qp_input.joint_pd_override.size()) !=
      qp_input.num_joint_pd_overrides) {
    std::cerr << "WARNING: num joint pd override doesn't match" << std::endl;
    qp_input.num_joint_pd_overrides = qp_input.joint_pd_override.size();
  }
  for (int i = 0; i < qp_input.num_support_data; ++i) {
    verifySubtypeSizes(qp_input.support_data[i]);
  }
}
