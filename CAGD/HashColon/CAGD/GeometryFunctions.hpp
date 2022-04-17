#ifndef HASHCOLON_CAGD_GEOMETRYFUNCTIONS_HPP
#define HASHCOLON_CAGD_GEOMETRYFUNCTIONS_HPP

#include <HashColon/CAGD/CAGDBasis.hpp>
#include <HashColon/CAGD/Geometry.hpp>

namespace HashColon::CAGD
{
	/* Compute control points of ruled patch:
	*   - uOrder / vOrder : order of basis
	*   - controlPoints : control points of quad surface 
	*    (ignore other control points except boundary points on u/v direction)
	*/ 
	template <size_t Dim>
	PointArrT<Dim> ComputeRuledPatch(size_t uOrder, size_t vOrder, const PointArrT<Dim>& controlPoints, bool isUdirection);

	/* Compute control points of bilinear patch:
	*   - uOrder / vOrder : order of basis
	*   - controlPoints : control points of quad surface
	*    (ignore other control points except 4 boundary vertex points)
	*/
	template <size_t Dim>
	PointArrT<Dim> ComputeBilinearPatch(size_t uOrder, size_t vOrder, const PointArrT<Dim>& controlPoints);

	/* Compute control points of coons patch:
	*   - uOrder / vOrder : order of basis
	*   - controlPoints : control points of quad surface
	*    (ignore off-boundary points)
	*/
	template <size_t Dim>
	PointArrT<Dim> ComputeCoonsPatch(size_t uOrder, size_t vOrder, const PointArrT<Dim>& controlPoints);
	
	/* Compute control points of coons patch:
	*   - c1 ~ c4: boundary curves in right-hand seq order 
	*   - basis: surface basis provided. if nullptr, create a new basis	
	*/
	template <size_t Dim, size_t uDegree, size_t vDegree>
	std::shared_ptr<QuadSurface<Dim, uDegree, vDegree>> ComputeCoonsPatch(
		const Curve<Dim, uDegree>& c1, const Curve<Dim, vDegree>& c2,
		const Curve<Dim, uDegree>& c3, const Curve<Dim, vDegree>& c4,
		const std::shared_ptr<QuadSurfaceBasis<uDegree, vDegree>> basis = nullptr);

	/* Compute control points of flat triangle patch:
	*   - uOrder / vOrder : order of basis
	*   - controlPoints : control points of quad surface
	*    (ignore other control points except 3 boundary vertex points)
	*/
	template <size_t Dim>
	PointArrT<Dim> ComputeTriLinearPatch(size_t order, const PointArrT<Dim>& controlPoints);

	/* Compute control points of rule tri-patch:
	*   - uOrder / vOrder : order of basis
	*   - controlPoints : control points of quad surface
	*    (ignore other control points except on-boundary points)
	*/
	template <size_t Dim>
	PointArrT<Dim> ComputeTriRuledPatch(size_t order, const PointArrT<Dim>& controlPoints, size_t direction);

	/* Compute control points of coons patch:
	*   - uOrder / vOrder : order of basis
	*   - controlPoints : control points of quad surface
	*    (ignore off-boundary points)
	*/
	template <size_t Dim>
	PointArrT<Dim> ComputeTriCoonsPatch(size_t order, const PointArrT<Dim>& controlPoints);

	/* Compute control points of coons patch:
	*   - c1 ~ c3: boundary curves in right-hand seq order
	*   - basis: surface basis provided. if nullptr, create a new basis
	*/
	template <size_t Dim, size_t Degree>
	std::shared_ptr<TriSurface<Dim, Degree>> ComputeTriCoonsPatch(
		const Curve<Dim, Degree>& cu, const Curve<Dim, Degree>& cv, const Curve<Dim, Degree>& cw,
		const std::shared_ptr<TriSurfaceBasis<Degree>> basis);
}

#endif	

#include <HashColon/CAGD/impl/GeometryFunctions_Impl.hpp>