#ifndef DRAKE_DRAKECORE_H
#define DRAKE_DRAKECORE_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Path.h"

namespace Drake {

  /** @defgroup vector_concept Vector<ScalarType> Concept
   * @ingroup concepts
   * @{
   * @brief Describes a (potentially structured) data which can be operated on as a (finite-dimensional) column vector
   *
   * @nbsp
   *
   * | Valid Expressions (which must be implemented) |  |
   * ------------------|-------------------------------------------------------------|
   * | RowsAtCompileTime  | defined as a static constant int (or enum).  Can be Eigen::Dynamic. |
   * | size_t size()      | only required if RowsAtCompileTime==Eigen::Dynamic |
   * | template<Derived> Vector(const Eigen::MatrixBase<Derived>&)  | constructor taking an Eigen object |
   * | template<Derived> Vector& operator=(const Eigen::MatrixBase<Derived>&)   | assignment operator from an Eigen object |
   * | Eigen::Matrix<ScalarType,RowsAtCompileTime,1> toEigen(const Vector<ScalarType>&) | non-member namespace method which converts to the Eigen type |
   *
   * @}
   */

  /** EigenVector<Rows>::type<ScalarType>
   * @brief provides an alias for Eigen::Matrix<ScalarType,Rows,1> which is templated on only a single argument (the ScalarType)
   * @concept{vector_concept}
   */
  template <int Rows>
  struct EigenVector {
    template <typename ScalarType> using type = Eigen::Matrix<ScalarType,Rows,1>;
  };

  /** NullVector<ScalarType>
   * @brief provides the empty vector (templated by ScalarType)
   * @concept{vector_concept}
   */
  template <typename ScalarType> using NullVector = Eigen::Matrix<ScalarType,0,1>;

  template <typename ScalarType, int Rows> Eigen::Matrix<ScalarType,Rows,1> toEigen(const Eigen::Matrix<ScalarType,Rows,1>& vec) { return vec; }

  /**
   * @brief whether or not the given type is an Eigen column vector
   */
  template <typename StateVector>
  struct is_eigen_vector : public std::false_type {};

  template <typename Scalar, int Rows, int Options, int MaxRows>
  struct is_eigen_vector<Eigen::Matrix<Scalar, Rows, 1, Options, MaxRows, 1>> : public std::true_type {};

  /** getRandomVector()
   * @brief Returns a random vector of the desired type using Eigen::Random()
   * @concept{vector_concept}
   */
  template <template <typename> class VecType>
  VecType<double> getRandomVector(void) {
    static_assert(VecType<double>::RowsAtCompileTime>=0,"still need to handle dynamic case");
    VecType<double> rand(Eigen::Matrix<double,VecType<double>::RowsAtCompileTime,1>::Random());
    return rand;
  }

  namespace internal {
    template<typename VecType, class Enable = void>
    struct SizeDispatch {
      static std::size_t eval(const VecType &vec) { return VecType::RowsAtCompileTime; }
    };

    template<typename VecType>
    struct SizeDispatch<VecType, typename std::enable_if<VecType::RowsAtCompileTime == Eigen::Dynamic>::type > {
      static std::size_t eval(const VecType &vec) { return vec.size(); }
    };
  }
  /** size()
   * @brief Evaluate the size of a Vector
   * @concept{vector_concept}
   *
   * @retval RowsAtCompileTime or the result of size() for dynamically sized vectors
   */
  template <typename VecType> std::size_t size(const VecType& vec) { return internal::SizeDispatch<VecType>::eval(vec); }

  /** getCoordinateName()
   * @brief Implements default coordinate names for a generic vector.  Overload this to name your coordinates.
   * @concept{vector_concept}
   */
  template <typename Vector>
  std::string getCoordinateName(const Vector& vec, unsigned int index) { return "x"+std::to_string(index); }

  /** operator<<()
   * @brief Implements the ostream operator for your vector types.
   * @concept{vector_concept}
   */
/*
  template <typename Vector>  // warning: this will match almost anything!
  std::ostream& operator<<(std::ostream& os, const Vector& vec)
  {
    for (int i=0; i<=size(vec); i++)
      os << getCoordinateName(vec,i) << " = " << vec(i) << std::endl;
    return os;
  }
*/
  /** CombinedVector<ScalarType,Vector1,Vector2>
   *
   * @brief produces a new vector type which is the columnwise composition of vector1 and vector2
   */
  template <typename ScalarType, template <typename> class Vector1, template <typename> class Vector2>
  class CombinedVector {
  public:
    CombinedVector() {};  // allow use of default constructors for vec1 and vec2, also
    CombinedVector(const Vector1<ScalarType>& first, const Vector2<ScalarType>& second) : vec1(first), vec2(second) {};

    template <typename Derived> CombinedVector(const Eigen::MatrixBase<Derived>& x) : vec1(x.topRows(Vector1<ScalarType>::RowsAtCompileTime)), vec2(x.bottomRows(Vector2<ScalarType>::RowsAtCompileTime)) {
      static_assert(RowsAtCompileTime != Eigen::Dynamic, "Cannot determine sizes of subvectors because sizes are not known at compile time."); // TODO: could handle cases where only one of the subvectors has dynamic size
    };

    template <typename Derived1, typename Derived2> CombinedVector(const Eigen::MatrixBase<Derived1>& x1, const Eigen::MatrixBase<Derived2>& x2) : vec1(x1), vec2(x2) {};

    template <typename Derived>
    CombinedVector& operator=(const Eigen::MatrixBase<Derived>& x) {
      vec1 = x.topRows(Drake::size(vec1));
      vec2 = x.bottomRows(Drake::size(vec2));
      return *this;
    }

    const Vector1<ScalarType>& first(void) const { return vec1; }
    const Vector2<ScalarType>& second(void) const { return vec2; }

    friend std::ostream& operator<<(std::ostream& os, const CombinedVector& x)
    {
      os << x.vec1 << std::endl;
      os << x.vec2 << std::endl;
      return os;
    }

    enum {
      RowsAtCompileTime = (Vector1<ScalarType>::RowsAtCompileTime == Eigen::Dynamic || Vector2<ScalarType>::RowsAtCompileTime == Eigen::Dynamic) ? Eigen::Dynamic : Vector1<ScalarType>::RowsAtCompileTime + Vector2<ScalarType>::RowsAtCompileTime
    };

    std::size_t size() const {
      return Drake::size(vec1) + Drake::size(vec2);
    }

    friend Eigen::Matrix<ScalarType,RowsAtCompileTime,1> toEigen(const CombinedVector<ScalarType,Vector1,Vector2>& vec) {
      Eigen::Matrix<ScalarType, RowsAtCompileTime, 1> x(vec.size());
      x << toEigen(vec.vec1), toEigen(vec.vec2);
      return x;
    }

  private:
    Vector1<ScalarType> vec1;
    Vector2<ScalarType> vec2;
  };

  /**
 * @brief whether or not the given type is a CombinedVector
 */
  template <typename StateVector>
  struct is_combined_vector : public std::false_type {};

  template <typename Scalar, template <typename> class Vector1, template <typename> class Vector2>
  struct is_combined_vector<CombinedVector<Scalar, Vector1, Vector2> > : public std::true_type {};

  namespace FunctionalForm {
    // these are meant as tags which can be used to inform algorithms about underlying structure in a function (or system)
    // e.g., linear, affine, polynomial, analytic, differentiable, continuous, measurable, and -- lastly -- arbritrary
    struct Arbitrary {};
    struct Polynomial : public Arbitrary {};
    struct Affine : public Polynomial {};
    struct Linear : public Affine {};
  }

}


#endif //DRAKE_DRAKECORE_H
