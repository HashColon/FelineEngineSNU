////#include <HashColon/Core/Exception.hpp>
////#include <HashColon/CAGD/GeometryFunctions.hpp>

namespace HashColon::CAGD
{
	//template<>
	//PointArrT<1> ComputeBilinearPatch<1>(size_t uOrder, size_t vOrder, const PointArrT<1>& controlPoints)
	//{
	//	using namespace std;
	//	Real uN = uOrder - 1.0;
	//	Real vN = vOrder - 1.0;

	//	// set return value
	//	PointArrT<1> bilinear = controlPoints;
	//	// set index map
	//	vector<vector<size_t>> indices = SurfaceIndexHelper::UIdx({ {uOrder, vOrder} });
	//	// set four corners
	//	array<size_t, 4> cornerIdx
	//	{ {
	//		indices.front().front(),
	//		indices.front().back(),
	//		indices.back().front(),
	//		indices.back().back()
	//	} };

	//	for (size_t i = 0; i < uOrder; i++)
	//	{
	//		for (size_t j = 0; j < vOrder; j++)
	//		{
	//			// compute CP
	//			Real tmp =
	//				((Real)((uN - i) * (vN - j)) / (Real)(uN * vN)) * bilinear(cornerIdx[0]) +
	//				((Real)((uN - i) * j) / (Real)(uN * vN)) * bilinear(cornerIdx[1]) +
	//				((Real)(i * (vN - j)) / (Real)(uN * vN)) * bilinear(cornerIdx[2]) +
	//				((Real)(i * j) / (Real)(uN * vN)) * bilinear(cornerIdx[3]);
	//			bilinear(indices[i][j]) = tmp;
	//		}
	//	}
	//	return bilinear;
	//}
}