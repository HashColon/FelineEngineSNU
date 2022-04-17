#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <HashColon/Core/Real.hpp>
#include <HashColon/CAGD/Geometry.hpp>
#include <HashColon/CAGD/GeometryFunctions.hpp>

using namespace std;
using namespace HashColon;
using namespace HashColon::CAGD;

int main()
{
	Eigen::Matrix<Real, 1, 16> controlpoints;
	controlpoints <<
		1, 3, 3, 1, 3, 1, 1, 3, 3, 2, 2, 3, 1, 3, 3, 1;

	Eigen::Matrix<Real, 1, -1> re = ComputeCoonsPatch<1>(4, 4, controlpoints);

	Eigen::Matrix<Real, 1, -1> cTest(4);
	cTest << 0, 3, 3, 0;

	Curve<1, 3> bc(make_shared<BezierBasis<3>>(), cTest);

	Real epsilon = 1e-12;
	Real step = 0.01;
	for (Real u = 0; u <= 1.0 + epsilon; u += step)
	{
		cout << bc.PointOnCurve(u) << "\t";
	}
	
	

			
}