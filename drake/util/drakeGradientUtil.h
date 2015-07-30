#ifndef DRAKEGRADIENTUTIL_H_
#define DRAKEGRADIENTUTIL_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <array>
#include <cassert>

#undef DLLEXPORT
#if defined(WIN32) || defined(WIN64)
  #if defined(drakeGradientUtil_EXPORTS)
    #define DLLEXPORT __declspec( dllexport )
  #else
    #define DLLEXPORT __declspec( dllimport )
  #endif
#else
  #define DLLEXPORT
#endif

template<std::size_t Size>
std::array<int, Size> intRange(int start)
{
  std::array<int, Size> ret;
  for (unsigned int i = 0; i < Size; i++) {
    ret[i] = i + start;
  }
  return ret;
}

/*
 * Recursively defined template specifying a matrix type of the correct size for a gradient of a matrix function with respect to Nq variables, of any order.
 */
template<typename Derived, int Nq, int DerivativeOrder = 1>
struct Gradient {
  typedef typename Eigen::Matrix<typename Derived::Scalar, ((Derived::SizeAtCompileTime == Eigen::Dynamic || Nq == Eigen::Dynamic) ? Eigen::Dynamic : Gradient<Derived, Nq, DerivativeOrder - 1>::type::SizeAtCompileTime), Nq> type;
};

/*
 * Base case for recursively defined gradient template.
 */
template<typename Derived, int Nq>
struct Gradient<Derived, Nq, 1> {
  typedef typename Eigen::Matrix<typename Derived::Scalar, Derived::SizeAtCompileTime, Nq> type;
};

/*
 * Output type of matGradMultMat
 */
template<typename DerivedA, typename DerivedB, typename DerivedDA>
struct MatGradMultMat {
  typedef typename Eigen::Matrix<typename DerivedA::Scalar, (DerivedA::RowsAtCompileTime == Eigen::Dynamic || DerivedB::ColsAtCompileTime == Eigen::Dynamic ? Eigen::Dynamic : DerivedA::RowsAtCompileTime * DerivedB::ColsAtCompileTime), DerivedDA::ColsAtCompileTime> type;
};

/*
 * Output type of matGradMult
 */
template<typename DerivedDA, typename DerivedB>
struct MatGradMult {
  typedef typename Eigen::Matrix<typename DerivedDA::Scalar, (DerivedDA::RowsAtCompileTime == Eigen::Dynamic || DerivedB::SizeAtCompileTime == Eigen::Dynamic ? Eigen::Dynamic : DerivedDA::RowsAtCompileTime / DerivedB::RowsAtCompileTime * DerivedB::ColsAtCompileTime), DerivedDA::ColsAtCompileTime> type;
};

/*
 * Output type and array types of getSubMatrixGradient for std::arrays specifying rows and columns
 */
template<int QSubvectorSize, typename Derived, std::size_t NRows, std::size_t NCols>
struct GetSubMatrixGradientArray {
  typedef typename Eigen::Matrix<typename Derived::Scalar, (NRows * NCols), ((QSubvectorSize == Eigen::Dynamic) ? Derived::ColsAtCompileTime : QSubvectorSize)> type;
};

/*
 * Output type of getSubMatrixGradient for a single element of the input matrix
 */
template<int QSubvectorSize, typename Derived>
struct GetSubMatrixGradientSingleElement {
  typedef typename Eigen::Block<const Derived, 1, ((QSubvectorSize == Eigen::Dynamic) ? Derived::ColsAtCompileTime : QSubvectorSize)> type;
};
 
/*
 * Profile results: looks like return value optimization works; a version that sets a reference
 * instead of returning a value is just as fast and this is cleaner.
 */
template <typename Derived>
typename Derived::PlainObject transposeGrad(
    const Eigen::MatrixBase<Derived>& dX, typename Derived::Index rows_X)
{
  typename Derived::PlainObject dX_transpose(dX.rows(), dX.cols());
  typename Derived::Index numel = dX.rows();
  typename Derived::Index index = 0;
  for (int i = 0; i < numel; i++) {
    dX_transpose.row(i) = dX.row(index);
    index += rows_X;
    if (index >= numel) {
      index = (index % numel) + 1;
    }
  }
  return dX_transpose;
}

template <typename DerivedA, typename DerivedB, typename DerivedDA, typename DerivedDB>
typename MatGradMultMat<DerivedA, DerivedB, DerivedDA>::type
matGradMultMat(
    const Eigen::MatrixBase<DerivedA>& A,
    const Eigen::MatrixBase<DerivedB>& B,
    const Eigen::MatrixBase<DerivedDA>& dA,
    const Eigen::MatrixBase<DerivedDB>& dB) {
  assert(dA.cols() == dB.cols());

  typename MatGradMultMat<DerivedA, DerivedB, DerivedDA>::type ret(A.rows() * B.cols(), dA.cols());

  for (int col = 0; col < B.cols(); col++) {
    auto block = ret.template block<DerivedA::RowsAtCompileTime, DerivedDA::ColsAtCompileTime>(col * A.rows(), 0, A.rows(), dA.cols());

    // A * dB part:
    block.noalias() = A * dB.template block<DerivedA::ColsAtCompileTime, DerivedDA::ColsAtCompileTime>(col * A.cols(), 0, A.cols(), dA.cols());

    for (int row = 0; row < B.rows(); row++) {
      // B * dA part:
      block.noalias() += B(row, col) * dA.template block<DerivedA::RowsAtCompileTime, DerivedDA::ColsAtCompileTime>(row * A.rows(), 0, A.rows(), dA.cols());
    }
  }
  return ret;

  // much slower and requires eigen/unsupported:
//  return Eigen::kroneckerProduct(Eigen::MatrixXd::Identity(B.cols(), B.cols()), A) * dB + Eigen::kroneckerProduct(B.transpose(), Eigen::MatrixXd::Identity(A.rows(), A.rows())) * dA;
}

template<typename DerivedDA, typename DerivedB>
typename MatGradMult<DerivedDA, DerivedB>::type
matGradMult(const Eigen::MatrixBase<DerivedDA>& dA, const Eigen::MatrixBase<DerivedB>& B) {
  assert(B.rows() == 0 ? dA.rows() == 0 : dA.rows() % B.rows() == 0);
  typename DerivedDA::Index A_rows = B.rows() == 0 ? 0 : dA.rows() / B.rows();
  const int A_rows_at_compile_time = (DerivedDA::RowsAtCompileTime == Eigen::Dynamic || DerivedB::RowsAtCompileTime == Eigen::Dynamic) ?
                                     Eigen::Dynamic :
                                     static_cast<int>(DerivedDA::RowsAtCompileTime / DerivedB::RowsAtCompileTime);

  typename MatGradMult<DerivedDA, DerivedB>::type ret(A_rows * B.cols(), dA.cols());
  ret.setZero();
  for (int col = 0; col < B.cols(); col++) {
    auto block = ret.template block<A_rows_at_compile_time, DerivedDA::ColsAtCompileTime>(col * A_rows, 0, A_rows, dA.cols());
    for (int row = 0; row < B.rows(); row++) {
      // B * dA part:
      block.noalias() += B(row, col) * dA.template block<A_rows_at_compile_time, DerivedDA::ColsAtCompileTime>(row * A_rows, 0, A_rows, dA.cols());
    }
  }
  return ret;
}

// TODO: could save copies once http://eigen.tuxfamily.org/bz/show_bug.cgi?id=329 is fixed
template<typename Derived>
DLLEXPORT Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> getSubMatrixGradient(
    const Eigen::MatrixBase<Derived>& dM, const std::vector<int>& rows, const std::vector<int>& cols,
    typename Derived::Index M_rows, int q_start = 0, typename Derived::Index q_subvector_size = -1) {
  if (q_subvector_size < 0) {
    q_subvector_size = dM.cols() - q_start;
  }
  Eigen::MatrixXd dM_submatrix(rows.size() * cols.size(), q_subvector_size);
  int index = 0;
  for (std::vector<int>::const_iterator col = cols.begin(); col != cols.end(); ++col) {
    for (std::vector<int>::const_iterator row = rows.begin(); row != rows.end(); ++row) {
      dM_submatrix.row(index) = dM.block(*row + *col * M_rows, q_start, 1, q_subvector_size);
      index++;
    }
  }
  return dM_submatrix;
}

template<int QSubvectorSize, typename Derived, std::size_t NRows, std::size_t NCols>
DLLEXPORT typename GetSubMatrixGradientArray<QSubvectorSize, Derived, NRows, NCols>::type
getSubMatrixGradient(const Eigen::MatrixBase<Derived>& dM,
  const std::array<int, NRows>& rows,
  const std::array<int, NCols>& cols,
  typename Derived::Index M_rows, int q_start = 0, typename Derived::Index q_subvector_size = QSubvectorSize) {
  if (q_subvector_size == Eigen::Dynamic) {
    q_subvector_size = dM.cols() - q_start;
  }
  Eigen::Matrix<typename Derived::Scalar, NRows * NCols, Derived::ColsAtCompileTime> dM_submatrix(NRows * NCols, q_subvector_size);
  int index = 0;
  for (typename std::array<int, NCols>::const_iterator col = cols.begin(); col != cols.end(); ++col) {
    for (typename std::array<int, NRows>::const_iterator row = rows.begin(); row != rows.end(); ++row) {
      dM_submatrix.row(index++) = dM.template block<1, QSubvectorSize> ((*row) + (*col) * M_rows, q_start, 1, q_subvector_size);
    }
  }
  return dM_submatrix;
};

template<int QSubvectorSize, typename Derived>
DLLEXPORT typename GetSubMatrixGradientSingleElement<QSubvectorSize, Derived>::type
getSubMatrixGradient(const Eigen::MatrixBase<Derived>& dM, int row, int col, typename Derived::Index M_rows,
    typename Derived::Index q_start = 0, typename Derived::Index q_subvector_size = QSubvectorSize) {
  if (q_subvector_size == Eigen::Dynamic) {
    q_subvector_size = dM.cols() - q_start;
  }
  return dM.template block<1, QSubvectorSize>(row + col * M_rows, q_start, 1, q_subvector_size);
};

template<typename DerivedA, typename DerivedB>
DLLEXPORT void setSubMatrixGradient(Eigen::MatrixBase<DerivedA>& dM, const Eigen::MatrixBase<DerivedB>& dM_submatrix,
    const std::vector<int>& rows, const std::vector<int>& cols, typename DerivedA::Index M_rows, typename DerivedA::Index q_start = 0, typename DerivedA::Index q_subvector_size = -1) {
  if (q_subvector_size < 0) {
    q_subvector_size = dM.cols() - q_start;
  }
  int index = 0;
  for (typename std::vector<int>::const_iterator col = cols.begin(); col != cols.end(); ++col) {
    for (typename std::vector<int>::const_iterator row = rows.begin(); row != rows.end(); ++row) {
      dM.block((*row) + (*col) * M_rows, q_start, 1, q_subvector_size) = dM_submatrix.row(index++);
    }
  }
};

template<int QSubvectorSize, typename DerivedA, typename DerivedB, std::size_t NRows, std::size_t NCols>
DLLEXPORT void setSubMatrixGradient(Eigen::MatrixBase<DerivedA>& dM, const Eigen::MatrixBase<DerivedB>& dM_submatrix,
    const std::array<int, NRows>& rows, const std::array<int, NCols>& cols, typename DerivedA::Index M_rows, typename DerivedA::Index q_start = 0, typename DerivedA::Index q_subvector_size = QSubvectorSize) {
  if (q_subvector_size == Eigen::Dynamic) {
    q_subvector_size = dM.cols() - q_start;
  }
  int index = 0;
  for (typename std::array<int, NCols>::const_iterator col = cols.begin(); col != cols.end(); ++col) {
    for (typename std::array<int, NRows>::const_iterator row = rows.begin(); row != rows.end(); ++row) {
      dM.template block<1, QSubvectorSize> ((*row) + (*col) * M_rows, q_start, 1, q_subvector_size) = dM_submatrix.row(index++);
    }
  }
};

template<int QSubvectorSize, typename DerivedDM, typename DerivedDMSub>
DLLEXPORT void setSubMatrixGradient(Eigen::MatrixBase<DerivedDM>& dM, const Eigen::MatrixBase<DerivedDMSub>& dM_submatrix, int row, int col, typename DerivedDM::Index M_rows,
    typename DerivedDM::Index q_start = 0, typename DerivedDM::Index q_subvector_size = QSubvectorSize) {
  if (q_subvector_size == Eigen::Dynamic) {
    q_subvector_size = dM.cols() - q_start;
  }
  dM.template block<1, QSubvectorSize>(row + col * M_rows, q_start, 1, q_subvector_size) = dM_submatrix;
};

#endif /* DRAKEGRADIENTUTIL_H_ */
