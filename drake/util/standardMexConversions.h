//
// Created by Twan Koolen on 10/10/15.
//

#ifndef DRAKE_STANDARDMEXCONVERSIONS_H
#define DRAKE_STANDARDMEXCONVERSIONS_H

#include "drakeMexUtil.h"

/**
 * fromMex specializations
 */

int fromMex(const mxArray* source, int*) {
  if (mxGetM(source) != 1 || mxGetN(source) != 1)
    throw MexToCppConversionError("Expected scalar.");
  return static_cast<int>(mxGetScalar(source));
}

bool fromMex(const mxArray* source, bool*) {
  if (!mxIsLogicalScalar(source))
    throw MexToCppConversionError("Expected logical.");
  return mxGetLogicals(source)[0];
}

/*
 * note: leaving Options, MaxRows, and MaxCols out as template parameters (leaving them to their defaults)
 * results in an internal compiler error on MSVC. See https://connect.microsoft.com/VisualStudio/feedback/details/1847159
 */
template <int Rows, int Cols, int Options, int MaxRows, int MaxCols>
Eigen::Map<const Eigen::Matrix<double, Rows, Cols, Options, MaxRows, MaxCols>> fromMex(const mxArray *mex, Eigen::MatrixBase<Eigen::Map<const Eigen::Matrix<double, Rows, Cols, Options, MaxRows, MaxCols>>> *) {
  if (!mxIsNumeric(mex)) {
    throw MexToCppConversionError("Expected a numeric array");
  }
  return matlabToEigenMap<Rows, Cols>(mex);
}

/*
 * note: leaving Options, MaxRows, and MaxCols out as template parameters (leaving them to their defaults)
 * results in an internal compiler error on MSVC. See https://connect.microsoft.com/VisualStudio/feedback/details/1847159
 */
template<typename DerType, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
Eigen::Matrix<Eigen::AutoDiffScalar<DerType>, Rows, Cols, Options, MaxRows, MaxCols> fromMex(
    const mxArray *mex, Eigen::MatrixBase<Eigen::Matrix<Eigen::AutoDiffScalar<DerType>, Rows, Cols, Options, MaxRows, MaxCols>> *) {
  if (!mxIsClass(mex, "TaylorVar"))
    throw MexToCppConversionError("Expected an array of TaylorVar");

  return taylorVarToEigen<Rows, Cols>(mex);
}

template<typename CoefficientType, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
Eigen::Matrix<TrigPoly<CoefficientType>, Rows, Cols, Options, MaxRows, MaxCols> fromMex(
    const mxArray *mex, Eigen::MatrixBase<Eigen::Matrix<TrigPoly<CoefficientType>, Rows, Cols, Options, MaxRows, MaxCols>> *) {
  if (!mxIsClass(mex, "TrigPoly"))
    throw MexToCppConversionError("Expected an array of TrigPoly");

  return trigPolyToEigen(mex);
}

/**
 * toMex specializations
 */

void toMex(const bool& source, mxArray* dest[], int nlhs) {
  dest[0] = mxCreateLogicalScalar(source);
}

/*
 * note: leaving Options, MaxRows, and MaxCols out as template parameters (leaving them to their defaults)
 * results in an internal compiler error on MSVC. See https://connect.microsoft.com/VisualStudio/feedback/details/1847159
 */
template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void toMex(const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> &source, mxArray *dest[], int nlhs) {
  if (nlhs > 0)
    dest[0] = eigenToMatlabGeneral(source);
};

template<typename Scalar, int Rows, int Cols>
void toMex(const GradientVar<Scalar, Rows, Cols> &source, mxArray *dest[], int nlhs, bool top_level = true) {
  if (top_level) {
    // check number of output arguments
    if (nlhs > source.maxOrder() + 1) {
      std::ostringstream buf;
      buf << nlhs << " output arguments desired, which is more than the maximum number: " << source.maxOrder() + 1 << ".";
      mexErrMsgTxt(buf.str().c_str());
    }
  }

  if (nlhs != 0) {
    // set an output argument
    toMex(source.value(), dest, nlhs);

    // recurse
    if (source.hasGradient()) {
      toMex(source.gradient(), &dest[1], nlhs - 1, false);
    }
  }
};

#endif //DRAKE_STANDARDMEXCONVERSIONS_H
