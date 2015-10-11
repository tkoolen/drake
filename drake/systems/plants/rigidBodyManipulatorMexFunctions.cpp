#include "rigidBodyManipulatorMexFunctions.h"
#include "RigidBodyManipulator.h"
#include "standardMexConversions.h"
#include <typeinfo>

using namespace std;
using namespace Eigen;

/**
 * fromMex specializations
 */
RigidBodyManipulator &fromMex(const mxArray *source, RigidBodyManipulator *) {
  return *static_cast<RigidBodyManipulator *>(getDrakeMexPointer(source));
}

std::set<int> fromMex(const mxArray *source, std::set<int> *) {
  // for robotnum. this is kind of a weird one, but OK for now
  std::set<int> robotnum_set;
  int num_robot = static_cast<int>(mxGetNumberOfElements(source));
  double *robotnum = mxGetPrSafe(source);
  for (int i = 0; i < num_robot; i++) {
    robotnum_set.insert((int) robotnum[i] - 1);
  }
  return robotnum_set;
}

template<typename Scalar>
KinematicsCache<Scalar> &fromMex(const mxArray *mex, KinematicsCache<Scalar> *) {
  if (!mxIsClass(mex, "DrakeMexPointer")) {
    throw MexToCppConversionError("Expected DrakeMexPointer containing KinematicsCache");
  }
  auto name = mxGetStdString(mxGetPropertySafe(mex, "name"));
  if (name != typeid(KinematicsCache<Scalar>).name()) {
    ostringstream buf;
    buf << "Expected KinematicsCache of type " << typeid(KinematicsCache<Scalar>).name() << ", but got " << name;
    throw MexToCppConversionError(buf.str());
  }
  return *static_cast<KinematicsCache<Scalar> *>(getDrakeMexPointer(mex));
}

/**
 * toMex specializations
 */
void toMex(const KinematicPath &path, mxArray *dest[], int nlhs) {
  if (nlhs > 0) {
    dest[0] = stdVectorToMatlab(path.body_path);
  }
  if (nlhs > 1) {
    dest[1] = stdVectorToMatlab(path.joint_path);
  }
  if (nlhs > 2) {
    dest[2] = stdVectorToMatlab(path.joint_direction_signs);
  }
}

/**
 * make_function
 * Note that a completely general make_function implementation is not possible due to ambiguities, but this works for all of the cases in this file
 * Inspired by from http://stackoverflow.com/a/21740143/2228557
 */

//plain function pointers
template<typename... Args, typename ReturnType>
auto make_function(ReturnType(*p)(Args...))
-> std::function<ReturnType(Args...)> { return {p}; }

//nonconst member function pointers
template<typename... Args, typename ReturnType, typename ClassType>
auto make_function(ReturnType(ClassType::*p)(Args...))
-> std::function<ReturnType(ClassType&, Args...)> { return {p}; }

//const member function pointers
template<typename... Args, typename ReturnType, typename ClassType>
auto make_function(ReturnType(ClassType::*p)(Args...) const)
-> std::function<ReturnType(const ClassType&, Args...)> { return {p}; }

typedef AutoDiffScalar<VectorXd> AutoDiffDynamicSize;
//typedef AutoDiffScalar<Matrix<double, Dynamic, 1, 0, 72, 1> AutoDiffFixedMaxSize;

/**
 * Mex function implementations
 */
void centerOfMassJacobianDotTimesVmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&RigidBodyManipulator::centerOfMassJacobianDotTimesV<double>);
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::centerOfMassJacobianDotTimesV<AutoDiffDynamicSize>);
//  auto func_trigpoly = make_function(&RigidBodyManipulator::centerOfMassJacobianDotTimesV<TrigPolyd>);
  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic); //, func_trigpoly);
}

void centerOfMassmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&RigidBodyManipulator::centerOfMass<double>);
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::centerOfMass<AutoDiffDynamicSize>);
//  auto func_trigpoly = make_function(&RigidBodyManipulator::centerOfMass<TrigPolyd>);
  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic); //, func_trigpoly);
}

void centerOfMassJacobianmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&RigidBodyManipulator::centerOfMassJacobian<double>);
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::centerOfMassJacobian<AutoDiffDynamicSize>);
//  auto func_trigpoly = make_function(&RigidBodyManipulator::centerOfMassJacobian<TrigPolyd>);
  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic); //, func_trigpoly);
}

void centroidalMomentumMatrixDotTimesvmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&RigidBodyManipulator::centroidalMomentumMatrixDotTimesV<double>);
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::centroidalMomentumMatrixDotTimesV<AutoDiffDynamicSize>);
//  auto func_trigpoly = make_function(&RigidBodyManipulator::centroidalMomentumMatrixDotTimesV<TrigPolyd>);
  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic); //, func_trigpoly);
}

void centroidalMomentumMatrixmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&RigidBodyManipulator::centroidalMomentumMatrix<double>);
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::centroidalMomentumMatrix<AutoDiffDynamicSize>);
//  auto func_trigpoly = make_function(&RigidBodyManipulator::centroidalMomentumMatrix<TrigPolyd>);
  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic); //, func_trigpoly);
}

template <typename DerivedQ, typename DerivedV>
void doKinematicsTemp(const RigidBodyManipulator &model, KinematicsCache<typename DerivedQ::Scalar> &cache, const MatrixBase<DerivedQ> &q, const MatrixBase<DerivedV> &v, bool compute_JdotV) {
  // temporary solution. Explicit doKinematics calls will not be necessary in the near future.
  if (v.size() == 0 && model.num_velocities != 0)
    cache.initialize(q);
  else
    cache.initialize(q, v);
  model.doKinematics(cache, compute_JdotV);
}

void doKinematicsmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&doKinematicsTemp<Map<const VectorXd>, Map<const VectorXd>>);

  typedef Matrix<AutoDiffScalar<VectorXd>, Dynamic, 1> VectorXAutoDiff;
  auto func_autodiff_dynamic = make_function(&doKinematicsTemp<VectorXAutoDiff, VectorXAutoDiff>);

  typedef Matrix<TrigPolyd, Dynamic, 1> VectorXTrigPoly;
  auto func_trigpoly = make_function(&doKinematicsTemp<VectorXTrigPoly, VectorXTrigPoly>);
  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic, func_trigpoly);
}

void findKinematicPathmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func = make_function(&RigidBodyManipulator::findKinematicPath);
  mexCallFunction(func, nlhs, plhs, nrhs, prhs);
}

void forwardJacDotTimesVmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  typedef Map<const Matrix3Xd> DerivedPointsDouble;
  auto func_double = make_function(&RigidBodyManipulator::forwardJacDotTimesV<DerivedPointsDouble>);

  typedef Matrix<AutoDiffDynamicSize, 3, Dynamic> DerivedPointsAutoDiff;
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::forwardJacDotTimesV<DerivedPointsAutoDiff>);

//  typedef Matrix<TrigPolyd, 3, Dynamic> DerivedPointsTrigPoly;
//  auto func_trigpoly = make_function(&RigidBodyManipulator::forwardJacDotTimesV<DerivedPointsTrigPoly>);

  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic); //, func_trigpoly);
}

void forwardKinmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  typedef Map<const Matrix3Xd> DerivedPointsDouble;
  auto func_double = make_function(&RigidBodyManipulator::forwardKin<DerivedPointsDouble>);

  typedef Matrix<AutoDiffDynamicSize, 3, Dynamic> DerivedPointsAutoDiff;
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::forwardKin<DerivedPointsAutoDiff>);

  // can't add this because of rotmat2quat, rotmat2rpy
//  typedef Matrix<TrigPolyd, 3, Dynamic> DerivedPointsTrigPoly;
//  auto func_trigpoly = make_function(&RigidBodyManipulator::forwardKin<DerivedPointsTrigPoly>);

  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic);
}

void forwardKinJacobianmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  typedef Map<const Matrix3Xd> DerivedPointsDouble;
  auto func_double = make_function(&RigidBodyManipulator::forwardKinJacobian<DerivedPointsDouble>);

  typedef Matrix<AutoDiffDynamicSize, 3, Dynamic> DerivedPointsAutoDiff;
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::forwardKinJacobian<DerivedPointsAutoDiff>);

//  typedef Matrix<TrigPolyd, 3, Dynamic> DerivedPointsTrigPoly;
//  auto func_trigpoly = make_function(&RigidBodyManipulator::forwardKinJacobian<DerivedPointsTrigPoly>);

  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic); //, func_trigpoly);
}

void forwardKinPositionGradientmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&RigidBodyManipulator::forwardKinPositionGradient<double>);
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::forwardKinPositionGradient<AutoDiffDynamicSize>);
//  auto func_trigpoly = make_function(&RigidBodyManipulator::forwardKinPositionGradient<TrigPolyd>);

  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic); //, func_trigpoly);
}

void geometricJacobianDotTimesVmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&RigidBodyManipulator::geometricJacobianDotTimesV<double>);
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::geometricJacobianDotTimesV<AutoDiffDynamicSize>);
  auto func_trigpoly = make_function(&RigidBodyManipulator::geometricJacobianDotTimesV<TrigPolyd>);

  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic, func_trigpoly);
}

template <typename Scalar>
GradientVar<Scalar, TWIST_SIZE, Eigen::Dynamic> geometricJacobianTemp(const RigidBodyManipulator &model, const KinematicsCache<Scalar> &cache, int base_body_or_frame_ind, int end_effector_body_or_frame_ind, int expressed_in_body_or_frame_ind, int gradient_order, bool in_terms_of_qdot) {
  // temporary solution. Gross v_or_qdot_indices pointer will be gone soon.
  return model.geometricJacobian(cache, base_body_or_frame_ind, end_effector_body_or_frame_ind, expressed_in_body_or_frame_ind, gradient_order, in_terms_of_qdot, nullptr);
};

void geometricJacobianmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&geometricJacobianTemp<double>);
  auto func_autodiff_dynamic = make_function(&geometricJacobianTemp<AutoDiffDynamicSize>);
  auto func_trigpoly = make_function(&geometricJacobianTemp<TrigPolyd>);
  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic, func_trigpoly);
}

void massMatrixmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&RigidBodyManipulator::massMatrix<double>);
  auto func_autodiff_dynamic = make_function(&RigidBodyManipulator::massMatrix<AutoDiffScalar<VectorXd>>);
  auto func_trigpoly = make_function(&RigidBodyManipulator::massMatrix<TrigPolyd>);
  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic, func_trigpoly);
}

template <typename Scalar, typename DerivedF, typename DerivedDF>
GradientVar<Scalar, Eigen::Dynamic, 1>  dynamicsBiasTermTemp(const RigidBodyManipulator &model, KinematicsCache<Scalar> &cache, const MatrixBase<DerivedF> &f_ext_value, const MatrixBase<DerivedDF> &f_ext_gradient, int gradient_order) {
  // temporary solution. GradientVar will disappear, obviating the need for the extra argument. integer body indices will be handled differently.

  unordered_map<const RigidBody *, GradientVar<Scalar, 6, 1> > f_ext;

  if (f_ext_value.size() > 0) {
    assert(f_ext_value.cols() == model.bodies.size());
    if (gradient_order > 0) {
      assert(f_ext_gradient.rows() == f_ext_value.size());
      assert(f_ext_gradient.cols() == model.num_positions + model.num_velocities);
    }

    for (DenseIndex i = 0; i < f_ext_value.cols(); i++) {
      GradientVar<Scalar, TWIST_SIZE, 1> f_ext_gradientvar(TWIST_SIZE, 1, model.num_positions + model.num_velocities, gradient_order);
      f_ext_gradientvar.value() = f_ext_value.col(i);

      if (gradient_order > 0) {
        f_ext_gradientvar.gradient().value() = f_ext_gradient.template middleRows<TWIST_SIZE>(TWIST_SIZE * i);
      }

      f_ext.insert({model.bodies[i].get(), f_ext_gradientvar});
    }
  }

  return model.dynamicsBiasTerm(cache, f_ext, gradient_order);
};

void dynamicsBiasTermmex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto func_double = make_function(&dynamicsBiasTermTemp<double, Map<const Matrix<double, 6, Dynamic> >, Map<const Matrix<double, Dynamic, Dynamic> > >);
  auto func_autodiff_dynamic = make_function(&dynamicsBiasTermTemp<AutoDiffDynamicSize, Matrix<AutoDiffDynamicSize, 6, Dynamic>, Matrix<AutoDiffDynamicSize, Dynamic, Dynamic> >);
  auto func_trigpoly = make_function(&dynamicsBiasTermTemp<TrigPolyd, Matrix<TrigPolyd, 6, Dynamic>, Matrix<TrigPolyd, Dynamic, Dynamic> >);

  mexTryToCallFunctions(nlhs, plhs, nrhs, prhs, func_double, func_autodiff_dynamic, func_trigpoly);
}
