#ifndef HASHCOLON_CAGD_GEOMETRYFUNCTIONS_IMPL
#define HASHCOLON_CAGD_GEOMETRYFUNCTIONS_IMPL

#include <HashColon/Core/Exception.hpp>
#include <HashColon/CAGD/GeometryFunctions.hpp>

namespace HashColon::CAGD
{
	HASHCOLON_NAMED_EXCEPTION_DEFINITION(GeometryFunction);

	template <size_t Dim>
	PointArrT<Dim> ComputeRuledPatch(size_t uOrder, size_t vOrder, const PointArrT<Dim>& controlPoints, bool isUdirection)
	{
		using namespace std;		
		// set return value
		PointArrT<Dim> ruled = controlPoints;

		// set infos for computations
		vector<vector<size_t>> indices;
		PointArrT<Dim> s(Dim, 1);	// start curve
		PointArrT<Dim> e(Dim, 1);	// end curve
		size_t order;

		// get control point index map 
		if (isUdirection)
		{
			indices = SurfaceIndexHelper::UIdx({ {uOrder, vOrder} });
			s.resize(Dim, uOrder); e.resize(Dim, uOrder);
			order = uOrder;
		}
		else
		{
			indices = SurfaceIndexHelper::VIdx({ {uOrder, vOrder} });
			s.resize(Dim, vOrder); e.resize(Dim, vOrder);
			order = vOrder;
		}

		size_t N = indices.size();
		Real n = (Real)N;

		// extract start/end curve
		for (size_t i = 0; i < order; i++)
			s.col(i) = controlPoints.col(indices.front()[i]);		
		for (size_t i = 0; i < order; i++)
			e.col(i) = controlPoints.col(indices.back()[i]);
		

		// fill out 2nd row to n-1 row
		for (size_t i = 1; i < N - 1; i++)
		{
			PointArrT<Dim> c = ((n - 1 - i) / (n - 1)) * s + ((i) / (n - 1)) * e;
			SetControlPointsByIndex<Dim>(ruled, c, indices[i]);
		}
		/*shared_ptr<PointArrT<Dim>> re = make_shared<PointArrT<Dim>>();
		(*re) = ruled;
		return re;*/
		return ruled;
	}

	template <size_t Dim>
	PointArrT<Dim> ComputeBilinearPatch(size_t uOrder, size_t vOrder, const PointArrT<Dim>& controlPoints)
	{
		using namespace std;
		Real uN = uOrder - 1.0;
		Real vN = vOrder - 1.0;

		// set return value
		PointArrT<Dim> bilinear = controlPoints;
		// set index map
		vector<vector<size_t>> indices = SurfaceIndexHelper::UIdx({ {uOrder, vOrder} });
		// set four corners
		array<size_t, 4> cornerIdx
			{ {
				indices.front().front(),
				indices.front().back(),				
				indices.back().front(),
				indices.back().back()
			} };

		for (size_t i = 0; i < uOrder; i++)
		{
			for (size_t j = 0; j < vOrder; j++)
			{
				// compute CP
				PointT<Dim> tmp =
					((Real)((uN - i) * (vN - j)) / (Real)(uN * vN)) * bilinear.col(cornerIdx[0]) +
					((Real)((uN - i) * j) / (Real)(uN * vN)) * bilinear.col(cornerIdx[1]) +
					((Real)(i * (vN - j)) / (Real)(uN * vN)) * bilinear.col(cornerIdx[2]) +
					((Real)(i * j) / (Real)(uN * vN)) * bilinear.col(cornerIdx[3]);
				bilinear.col(indices[i][j]) = tmp;
			}
		}
		/*shared_ptr<PointArrT<Dim>> re = make_shared<PointArrT<Dim>>();
		(*re) = bilinear;
		return re;*/
		return bilinear;
	}
			
	template <size_t Dim>
	PointArrT<Dim> ComputeCoonsPatch(size_t uOrder, size_t vOrder, const PointArrT<Dim>& controlPoints)
	{		
		// Build u-ruled surface		
		PointArrT<Dim> uRuled = ComputeRuledPatch<Dim>(uOrder, vOrder, controlPoints, true /*u*/);
		PointArrT<Dim> vRuled = ComputeRuledPatch<Dim>(uOrder, vOrder, controlPoints, false /*v*/);
		PointArrT<Dim> bilinear = ComputeBilinearPatch<Dim>(uOrder, vOrder, controlPoints);
	
		return  uRuled + vRuled - bilinear;
	}

	/*
	* Build Coons Patch
	* input curves should be sequencial(ccw)
	*/
	template <size_t Dim, size_t uDegree, size_t vDegree>
	std::shared_ptr<QuadSurface<Dim, uDegree, vDegree>> ComputeCoonsPatch(
		const Curve<Dim, uDegree>& _c1, const Curve<Dim, vDegree>& _c2,
		const Curve<Dim, uDegree>& _c3, const Curve<Dim, vDegree>& _c4,
		const std::shared_ptr<QuadSurfaceBasis<uDegree, vDegree>> basis)
	{	
		using namespace std;
		
		Curve<Dim, uDegree> c1 = _c1;
		Curve<Dim, vDegree> c2 = _c2;
		Curve<Dim, uDegree> c3 = _c3;
		Curve<Dim, vDegree> c4 = _c4;
		
		// Check input curves 
		// Basis type/order of c1 & c3 should match
		if( typeid(c1.Basis()) != typeid(c3.Basis())) throw GeometryFunctionException("Basis type of u-direction curves do not match.");
		if (c1.Basis()->Order() != c3.Basis()->Order()) throw GeometryFunctionException("Order of u-direction curves do not match.");

		// Basis type/order of c2 & c4 should match
		if (typeid(c2.Basis()) != typeid(c4.Basis())) throw GeometryFunctionException("Basis type of v-direction curves do not match.");
		if (c2.Basis()->Order() != c4.Basis()->Order()) throw GeometryFunctionException("Order of v-direction curves do not match.");

		// End points of the curves should match		
		const size_t uBegin = 0, uEnd = c1.Basis()->Order() - 1;
		const size_t vBegin = 0, vEnd = c2.Basis()->Order() - 1;	

		// c1[end] -> c2[begin]
		if (c1.CP.col(uEnd) != c2.CP.col(vBegin)) 
		{
			// if c2 is in opposite direction( c1[end] == c2[end] ), flip it / else fuck it.
			if (c1.CP.col(uEnd) == c2.CP.col(vEnd)) { c2 = c2.Flip(); }
			else if (c1.CP.col(uBegin) == c2.CP.col(vBegin)) { c1 = c1.Flip(); }
			else if (c1.CP.col(uBegin) == c2.CP.col(vEnd)) { c1 = c1.Flip(); c2 = c2.Flip(); }
			else 
				throw GeometryFunctionException("End point of c1 does not match with start point of c2");
		}		
		// c2[end] -> c3[begin]
		if (c2.CP.col(vEnd) != c3.CP.col(uBegin))
		{
			// if c3 is in opposite direction( c2[end] == c3[end] ), flip it / else fuck it.
			if (c2.CP.col(vEnd) == c3.CP.col(uEnd)) { c3 = c3.Flip(); }
			else
				throw GeometryFunctionException("End point of c2 does not match with start point of c3");			
		}		
		// c3[end] -> c4[begin]
		if (c3.CP.col(uEnd) != c4.CP.col(vBegin))
		{
			// if c4 is in opposite direction( c3[end] == c4[end] ), flip it / else fuck it.
			if (c3.CP.col(uEnd) == c4.CP.col(vEnd)) { c4 = c4.Flip(); }
			else throw GeometryFunctionException("End point of c3 does not match with start point of c4");			
		}		
		// c4[end] -> c1[begin]
		if (c4.CP.col(vEnd) != c1.CP.col(uBegin)) {			
			// we do not double check the direction of c1, since c1 was where we started.
			throw GeometryFunctionException("End point of c4 does not match with start point of c1");
		}
			

		// set control points 		
		PointArrT<Dim> cplist(Dim, c1.Basis()->Order() * c2.Basis()->Order());
		auto bidx = SurfaceIndexHelper::BIdx({ {c1.Basis()->Order(), c2.Basis()->Order()} });
		SetControlPointsByIndex<Dim>(cplist, c1.CP, bidx[0]);
		SetControlPointsByIndex<Dim>(cplist, c2.CP, bidx[1]);
		SetControlPointsByIndex<Dim>(cplist, c3.CP, bidx[2]);
		SetControlPointsByIndex<Dim>(cplist, c4.CP, bidx[3]);
		
		shared_ptr<QuadSurfaceBasis<uDegree, vDegree>> newbasis;
		if (basis) newbasis = basis;
		else newbasis = make_shared<QuadSurfaceBasis<uDegree, vDegree>>(c1.Basis(), c2.Basis());
				
		return make_shared<QuadSurface<Dim, uDegree, vDegree>>(
			newbasis, ComputeCoonsPatch<Dim>(c1.Basis()->Order(), c2.Basis()->Order(), cplist));
	}

	template <size_t Dim>
	PointArrT<Dim> ComputeTriLinearPatch(size_t order, const PointArrT<Dim>& controlPoints)
	{	
		const size_t& n = order;
		PointArrT<Dim> lin = controlPoints;

		// set end points
		PointT<Dim> U = controlPoints.col((n * (n + 1)) / 2 - 1);
		PointT<Dim> V = controlPoints.col(n - 1);
		PointT<Dim> W = controlPoints.col(0);		

		size_t idx = 0;
		for (size_t i = 0; i < n; i++)
		{
			for (size_t j = 0; j < n - i; j++)
			{
				Real u = (Real)i / (n - 1.0);
				Real v = (Real)j / (n - 1.0);
				Real w = (n - i - j - 1.0) / (n - 1.0);				

				lin.col(idx) = u * U + v * V + w * W;
				idx++;
			}
		}
		return lin;
	}

	template <size_t Dim>
	PointArrT<Dim> ComputeTriRuledPatch(size_t order, const PointArrT<Dim>& controlPoints, size_t direction)
	{
		using namespace std;		
		PointArrT<Dim> rule = controlPoints;

		// Set curves for ruling
		// Set previous curve (reverse sequence)
		vector<size_t> prevIdx = TriSurfaceIndexHelper::BIdx(order)[(direction + 2) % 3];
		PointArrT<Dim> prev(Dim, order);
		for (size_t i = 0; i < order; i++) { prev.col(order - 1 - i) = controlPoints.col(prevIdx[i]); }
		// Set next curve
		vector<size_t> nextIdx = TriSurfaceIndexHelper::BIdx(order)[(direction + 1) % 3];
		PointArrT<Dim> next(Dim, order);
		for (size_t i = 0; i < order; i++) { next.col(i) = controlPoints.col(nextIdx[i]); }

		vector<vector<size_t>> ruleIdx =
			direction == TriSurfaceIndexHelper::cu ? TriSurfaceIndexHelper::UIdx(order) :
			direction == TriSurfaceIndexHelper::cv ? TriSurfaceIndexHelper::VIdx(order) :
			direction == TriSurfaceIndexHelper::cw ? TriSurfaceIndexHelper::WIdx(order) :
			throw CAGD_GeometryException("Invalid direction name given.");		

		for (size_t i = 0; i < ruleIdx.size(); i++)
		{
			for (size_t j = 0; j < ruleIdx[i].size() - 1; j++)
			{
				Real n = ruleIdx[i].size() - 1.0;
				rule.col(ruleIdx[i][j]) = ((n - j) / n)* prev.col(i) + (j / n) * next.col(i);
			}
		}
		return rule;
	}

	template <size_t Dim>
	PointArrT<Dim> ComputeTriCoonsPatch(size_t order, const PointArrT<Dim>& controlPoints)
	{	
		using namespace std;		
		auto uRuled = ComputeTriRuledPatch<Dim>(order, controlPoints, 0/*u*/);		
		auto vRuled = ComputeTriRuledPatch<Dim>(order, controlPoints, 1/*v*/);
		auto wRuled = ComputeTriRuledPatch<Dim>(order, controlPoints, 2/*w*/);
		auto linear = ComputeTriLinearPatch<Dim>(order, controlPoints);

		return uRuled + vRuled + wRuled - linear - linear;
	}

	/*
	* Build Coons Patch
	* input curves should be sequencial(ccw)
	*/
	template <size_t Dim, size_t Degree>
	std::shared_ptr<TriSurface<Dim, Degree>> ComputeTriCoonsPatch(
		const Curve<Dim, Degree>& _cu, const Curve<Dim, Degree>& _cv, const Curve<Dim, Degree>& _cw,
		const std::shared_ptr<TriSurfaceBasis<Degree>> basis)
	{
		using namespace std;

		Curve<Dim, Degree> cu = _cu;
		Curve<Dim, Degree> cv = _cv;
		Curve<Dim, Degree> cw = _cw;

		// Check input curves 
		// Basis type/order of c1 & c3 should match
		if ((typeid(cu.Basis()) != typeid(cv.Basis())) || (typeid(cv.Basis()) != typeid(cw.Basis())))
			throw GeometryFunctionException("Basis type of the boundary curves do not match.");
		if ((cu.Basis()->Order() != cv.Basis()->Order()) || (cu.Basis()->Order() != cv.Basis()->Order()))
			throw GeometryFunctionException("Order of the boundary curves do not match.");
				
		size_t order = cu.Basis()->Order();

		// End points of the curves should match		
		const size_t begin = 0, end = order - 1;		

		// cu[end] -> cv[begin]
		if (cu.CP.col(end) != cv.CP.col(begin))
		{
			// if cv is in opposite direction( c1[end] == c2[end] ), flip it / else fuck it.
			if (cu.CP.col(end) == cv.CP.col(end)) { cv = cv.Flip(); }
			else if (cu.CP.col(begin) == cv.CP.col(begin)) { cu = cu.Flip(); }
			else if (cu.CP.col(begin) == cv.CP.col(end)) { cu = cu.Flip(); cv = cv.Flip(); }
			else throw GeometryFunctionException("End point of cu does not match with start point of cv");
		}
		// cv[end] -> cw[begin]
		if (cv.CP.col(end) != cw.CP.col(begin))
		{
			// if cw is in opposite direction( cv[end] == cw[end] ), flip it / else fuck it.
			if (cv.CP.col(end) == cw.CP.col(end)) { cw = cw.Flip(); }
			else throw GeometryFunctionException("End point of cv does not match with start point of cw");			
		}		
		// cw[end] -> cu[begin]
		if (cw.CP.col(end) != cu.CP.col(begin))
			throw GeometryFunctionException("End point of cw does not match with start point of cv");		

		// set control points 		
		PointArrT<Dim> cplist(Dim, order * (order + 1) / 2);
		auto bidx = TriSurfaceIndexHelper::BIdx(order);
		SetControlPointsByIndex<Dim>(cplist, cu.CP, bidx[0]);
		SetControlPointsByIndex<Dim>(cplist, cv.CP, bidx[1]);
		SetControlPointsByIndex<Dim>(cplist, cw.CP, bidx[2]);
		
		return make_shared<TriSurface<Dim, Degree>>(
			basis, ComputeTriCoonsPatch<Dim>(basis->Order(), cplist));
	}
}

#endif