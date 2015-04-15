#include "PiecewisePolynomial.h"
#include <stdexcept>
#include <cassert>
#include <algorithm>

using namespace std;

PiecewisePolynomial::PiecewisePolynomial(std::vector<Polynomial> const& polynomials, std::vector<double> const& segment_times) :
    PiecewisePolynomialBase(segment_times),
    polynomials(polynomials)
{
  assert(segment_times.size() == (polynomials.size() + 1));

  for (int i = 0; i < polynomials.size(); i++) {
    if (segment_times[i + 1] < segment_times[i])
      throw runtime_error("times must be increasing");
  }
}

int PiecewisePolynomial::getSegmentIndex(double t) {
  // clip to min/max times
  if (t < getStartTime())
    t = getStartTime();

  if (t > getEndTime())
    t = getEndTime();

  int segment_index = 0;
  while (t >= getEndTime(segment_index) && segment_index < getNumberOfSegments() - 1)
    segment_index++;

  return segment_index;
}

PiecewisePolynomial PiecewisePolynomial::derivative(int derivative_order) const {
  vector<Polynomial> derivative_polynomials;
  derivative_polynomials.reserve(polynomials.size());
  for (auto it = polynomials.begin(); it != polynomials.end(); ++it) {
    derivative_polynomials.push_back(it->derivative(derivative_order));
  }
  return PiecewisePolynomial(derivative_polynomials, segment_times);
}

PiecewisePolynomial PiecewisePolynomial::integral(double integration_constant) const {
  vector<Polynomial> integral_polynomials;
  integral_polynomials.reserve(polynomials.size());
  for (auto it = polynomials.begin(); it != polynomials.end(); ++it) {
    integral_polynomials.push_back(it->integral(integration_constant));
  }
  return PiecewisePolynomial(integral_polynomials, segment_times);
}

double PiecewisePolynomial::value(double t) {
  int segment_index = getSegmentIndex(t);
  return polynomials[segment_index].value(t - getStartTime(segment_index));
}

const Polynomial& PiecewisePolynomial::getPolynomial(int segment_index) const {
  segmentNumberRangeCheck(segment_index);
  return polynomials[segment_index];
}

int PiecewisePolynomial::getSegmentPolynomialOrder(int segment_index) const {
  segmentNumberRangeCheck(segment_index);
  return polynomials[segment_index].getOrder();
}

PiecewisePolynomial& PiecewisePolynomial::operator+=(const PiecewisePolynomial& other) {
  if (!segmentTimesEqual(other, 1e-10))
    throw runtime_error("Addition not yet implemented when segment times are not equal");
  for (int i = 0; i < polynomials.size(); i++)
    polynomials[i] += other.polynomials[i];
  return *this;
}

PiecewisePolynomial& PiecewisePolynomial::operator*=(const PiecewisePolynomial& other) {
  if (!segmentTimesEqual(other, 1e-10))
    throw runtime_error("Multiplication not yet implemented when segment times are not equal");
  for (int i = 0; i < polynomials.size(); i++)
    polynomials[i] *= other.polynomials[i];
  return *this;
}

const PiecewisePolynomial PiecewisePolynomial::operator+(const PiecewisePolynomial &other) const {
  PiecewisePolynomial ret = *this;
  ret += other;
  return ret;
}

const PiecewisePolynomial PiecewisePolynomial::operator*(const PiecewisePolynomial &other) const {
  PiecewisePolynomial ret = *this;
  ret *= other;
  return ret;
}