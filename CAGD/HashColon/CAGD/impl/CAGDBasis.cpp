#ifndef EIGEN_INITIALIZE_MATRICES_BY_ZERO
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#endif

#include <vector>
#include <array>
#include <HashColon/Core/Real.hpp>
#include <HashColon/CAGD/CAGDBasis.hpp>

using namespace std;

/* SurfaceBasis */
namespace HashColon::CAGD
{	
	Eigen::MatrixXR& SurfaceBasis::PreCompute(
		const std::vector<std::array<Real, 2>>& uv,
		bool noThrow)
	{
		preComputedValue = make_shared<Eigen::MatrixXR>(Value(uv, noThrow));
		return *preComputedValue;
	}
	
	const Eigen::MatrixXR& SurfaceBasis::PreComputedValue() const 
	{
		if (!preComputedValue) throw Exception("Pre-computed value not computed yet");
		else return *preComputedValue;
	}
}

namespace HashColon::CAGD::_helper
{
	Eigen::VectorXR FlattenMatrixColwise(Eigen::MatrixXR A)
	{
		using namespace Eigen;
		return VectorXR(Map<VectorXR>(A.data(), A.cols() * A.rows()));
	}

	Eigen::VectorXR FlattenMatrixRowwise(Eigen::MatrixXR A)
	{
		using namespace Eigen;
		Eigen::MatrixXR B = A.transpose();
		return VectorXR(Map<VectorXR>(B.data(), B.cols() * B.rows()));
	}
}