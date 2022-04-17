#ifndef HASHCOLON_CAGD_GEOMETRY_HPP
#define HASHCOLON_CAGD_GEOMETRY_HPP

#include <vector>
#include <array>
#include <map>
#include <Eigen/Eigen>
#include <HashColon/Core/Exception.hpp>
#include <HashColon/Core/Real.hpp>
#include <HashColon/CAGD/CAGDBasis.hpp>

/* Types */
namespace HashColon::CAGD
{
	// point as a column vector
	template <size_t Dim>
	using PointT = Eigen::Matrix<HashColon::Real, Dim, 1>;

	/* Stores points as a list of columns
	* ex)	| x x x ... x |
	*		| y y y ... y |
	*		| z z z ... z |	*/
	template <size_t Dim>
	using PointArrT = Eigen::Matrix<HashColon::Real, Dim, -1>;
}

/* GeometryBase */
namespace HashColon::CAGD
{
	/*
	* GeometryBase: Generic base class for all geometry which can be defined with feature points
	*	ex)	Curves/Surfaces: control points
	*		Square: 4 points, Min-max box: 2 points, Circle: 1 point, etc....
	*	Feature points are stored as a columnwise matrix form. Each column refers to a feature point
	*/
	template <size_t Dim>
	class GeometryBase
	{
	public:
		PointArrT<Dim> FeaturePoints;

		//EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	protected:
		GeometryBase(
			const PointArrT<Dim>& featurePoints = PointArrT<Dim>::Zero(Dim, 1)
		) : FeaturePoints(featurePoints)
		{};

		GeometryBase(const GeometryBase& rhs) 
			: FeaturePoints(rhs.FeaturePoints)
		{};

		GeometryBase& operator=(const GeometryBase& rhs) { FeaturePoints = rhs.FeaturePoints; };
	};
}

 /* Curves */
namespace HashColon::CAGD
{
	template <size_t Dim, size_t Degree>
	class Curve : public GeometryBase<Dim>
	{
	protected:
		std::shared_ptr<CAGDBasisBase<Degree>> basis;
	public:
		HASHCOLON_CLASS_EXCEPTION_DEFINITION(Curve);		

		Curve(const std::shared_ptr<CAGDBasisBase<Degree>>& iBasis);

		Curve(
			const std::shared_ptr<CAGDBasisBase<Degree>>& iBasis, 
			const PointArrT<Dim>& featurePoints);

		Curve(const Curve& rhs)
			: GeometryBase<Dim>(rhs.FeaturePoints), basis(rhs.basis)
		{};

		Curve& operator=(const Curve& rhs) { this->FeaturePoints = rhs.FeaturePoints; basis = rhs.basis; return *this; };

		const std::shared_ptr<CAGDBasisBase<Degree>> Basis() const { return basis; };

		// alias for feature points. We call them "Control points" here.
		PointArrT<Dim>& ControlPoints = this->FeaturePoints;
		PointArrT<Dim>& CP = this->FeaturePoints;

		PointT<Dim> PointOnCurve(Real u) const;
		PointArrT<Dim> PointOnCurve_PreComputed() const;

		// Flip the sequence of the control points
		Curve<Dim, Degree> Flip();

		// Subdivision at parameter u
		//Curve<Dim, Degree> Subdivide(Real u, bool isForePart);

		// Clipping from parameter u1 to u2
		//Curve<Dim, Degree> Clip(Real u1, Real u2);

		// Degree elevation/reduction 
		/*template <size_t toDegree>
		Curve<Dim, toDegree> ElevateDegree();*/
		/*template <size_t toDegree>
		Curve<Dim, toDegree> ReduceDegree();*/
	};
}

/* Surfaces */
namespace HashColon::CAGD
{
	/*
	* Surface: Generic square surface class 
	* 
	* - The class follows the generic rule of feature points from the parent class
	* - The class does NOT use the common matrix form uCv^T. Instead, uses [Cij][uivi] form.
	*	This is because, it is difficult to handle matrix of points C, which should be treated as 3D tensor.
	*	The example of computing point on surface from order 4 basis is shown below.
	* 
	*	ex) uCv^T
	*	             | C00 C01 C02 C03 | |v0|
	*	[u0 u1 u2 u3]| C10 C11 C12 C13 | |v1|   where ui, vj are basis function of u, v directions, 
	*	             | C20 C21 C22 C23 | |v2|   Cij are control points.
	*	             | C30 C31 C32 C33 | |v3|
	* 
	*	ex) [Cij][uivi]
	*   | C00x C01x C02x C03x C10x C11x C12x C13x ... C32x C33x || u0 * v0 |
	*   | C00y C01y C02y C03y C10x C11x C12x C13x ... C32x C33x || u0 * v1 |
	*                      ...16 points(col) ...                 | u0 * v2 |
	*                                                            | u0 * v3 |
	*                                                            | u1 * v0 |
	*                                                            |   ...   |
	*                                                            | 16 rows |
	*                                                            |   ...   |
	*                                                            | u3 * v2 |
	*                                                            | u3 * v3 |
	* All control points are stored as the form described in [Cij] form
	*/
	template <size_t Dim>
	class Surface : public GeometryBase<Dim>
	{
	protected:
		std::shared_ptr<SurfaceBasis> basis;		

	public:
		HASHCOLON_CLASS_EXCEPTION_DEFINITION(Surface);

		Surface(const std::shared_ptr<SurfaceBasis>& basis);
		Surface(const std::shared_ptr<SurfaceBasis>& basis, const PointArrT<Dim>& featurePoints);
		Surface(const Surface& rhs);

		const std::shared_ptr<SurfaceBasis> Basis() { return basis; };

		// alias for feature points. We call them "Control points" here.
		PointArrT<Dim>& ControlPoints = this->FeaturePoints;
		PointArrT<Dim>& CP = this->FeaturePoints;

		/* point on surface */
		PointT<Dim> PointOnSurface(Real u, Real v) const;		

		/* Compute multiple point on Surface */
		PointArrT<Dim> MeshOnSurface(std::vector<std::array<Real, 2>> uv, bool noThrow = true) const;
		PointArrT<Dim> MeshOnSurface_PreComputed() const;		
	};

	template <size_t Dim, size_t uDegree, size_t vDegree>
	class QuadSurface : public Surface<Dim>
	{	
	public:		
		HASHCOLON_CLASS_EXCEPTION_DEFINITION(Surface);

		QuadSurface(std::shared_ptr<QuadSurfaceBasis<uDegree, vDegree>> ibasis);
		QuadSurface(std::shared_ptr<QuadSurfaceBasis<uDegree, vDegree>> ibasis,
			const PointArrT<Dim>& featurePoints);
		QuadSurface(const Surface<Dim>& rhs);
				
		/* Surface CP idxs: CP in surface is stored as following
		* (for Order 4x4 example)
		* | C00 C01 C02 C03 |	|  0  1  2  3 |
		* | C10 C11 C12 C13 | = |  4  5  6  7 |
		* | C20 C21 C22 C23 |   |  8  9 10 11 |
		* | C30 C31 C32 C33 |	| 12 13 14 15 |
		*
		* Index/Alias of each boundary curves are as follows
		* 0/u1:  0  4  8 12
		* 1/v1: 12 13 14 15
		* 2/u2: 15 11  7  3
		* 3/v2:  3  2  1  0
		*/

		enum BoundaryCurveAlias { u1 = 0, v1 = 1, u2 = 2, v2 = 3 };
	};
	
	/* TriSurface: Generic triangular surface class
	* 
	* Control points in this class are indexed as followings.
	* ex) for(order 4)
	*       u       
	*    0  1  2  3 
	*    4  5  6
	* v  7  8    w
	*    9
	* 
	* point u: [3, 0, (0)]: 9
	* point v: [0, 3, (0)]: 3
	* point w: [0, 0, (3)]: 0	* 
	* 
	* all index mapping is based on right-hand law;
	* for example, boundary index mapping goes like 
	*   cu: 3 2 1 0
	*   cv: 0 4 7 9
	*   cw: 9 8 6 3	
	* for example, u-row index mapping goes like
	*   row0: 3 2 1 0 (w->u->v direction)
	*   row1: 6 5 4
	*   row2: 8 7
	*   row3: 9
	* for example, v-row index mapping goes like
	*   row0: 0 4 7 9 (u->v->w direction)
	*   row1: 1 5 8
	*   row2: 2 6
	*   row3: 3
	* for example, w-row index mapping goes like
	*   row0: 9 8 6 3 (v->w->u direction)
	*   row1: 7 5 28
	*   row2: 4 1
	*   row3: 0
	*/
	template <size_t Dim, size_t Degree>
	class TriSurface : public Surface<Dim>
	{
	protected:
		std::shared_ptr<TriSurfaceBasis<Degree>> basis;
	public:
		HASHCOLON_CLASS_EXCEPTION_DEFINITION(TriSurface);

		TriSurface(const std::shared_ptr<TriSurfaceBasis<Degree>>& ibasis);
		TriSurface(const std::shared_ptr<TriSurfaceBasis<Degree>>& ibasis, const PointArrT<Dim>& featurePoints);
		TriSurface(const TriSurface<Dim, Degree>& rhs);

		const std::shared_ptr<TriSurfaceBasis<Degree>> Basis() { return basis; };

		// alias for feature points. We call them "Control points" here.
		PointArrT<Dim>& ControlPoints = this->FeaturePoints;
		PointArrT<Dim>& CP = this->FeaturePoints;

		enum BoundaryCurveAlias { cu = 0, cv = 1, cw = 2 };
		/* point on surface */
		//PointT<Dim> PointOnSurface(Real u, Real v) const;

		/* Compute multiple point on Surface */
		//PointArrT<Dim> MeshOnSurface(std::vector<std::array<Real, 2>> uv, bool noThrow = true) const;
		//PointArrT<Dim> MeshOnSurface_PreComputed() const;		
	};	
}

namespace HashColon::CAGD
{
	HASHCOLON_NAMED_EXCEPTION_DEFINITION(CAGD_Geometry);

	/* Flatten matrix to a vector (rowwise)
	* ex)	|  1  2  3  4 |
	*		|  5  6  7  8 |
	*		|  9 10 11 12 |
	*		| 13 14 15 16 | */
	Eigen::VectorXR FlattenMatrixColwise(Eigen::MatrixXR A);
	/* Flatten matrix to a vector (rowwise)
	* ex)	| 1 5  9 13 |
	*		| 2 6 10 14 |
	*		| 3 7 11 15 |
	*		| 4 8 12 16 | */
	Eigen::VectorXR FlattenMatrixRowwise(Eigen::MatrixXR A);

	/* Does the shitty & complecated building tasks
	* & save it for reuse */
	class SurfaceIndexHelper
	{
	private:
		static inline std::map<
			std::array<size_t, 2>,
			std::array<std::vector<size_t>, 4>
			> _bidxlist;

		static inline std::map<
			std::array<size_t, 2>,
			std::vector<std::vector<size_t>>
		> _uidxlist;

		static inline std::map<
			std::array<size_t, 2>,
			std::vector<std::vector<size_t>>
		> _vidxlist;


	public:
		enum BoundaryCurveAlias { u1 = 0, v1 = 1, u2 = 2, v2 = 3 };
		static std::array<std::vector<size_t>, 4> BIdx(std::array<size_t, 2> uvorders);
		static std::vector<std::vector<size_t>> UIdx(std::array<size_t, 2> uvorders);
		static std::vector<std::vector<size_t>> VIdx(std::array<size_t, 2> uvorders);
	};

	template <size_t Dim>
	PointArrT<Dim>& SetControlPointsByIndex(
		PointArrT<Dim>& SurfCP, const PointArrT<Dim>& CPs, 
		const std::vector<size_t>& indices);

	class TriSurfaceIndexHelper
	{
	private:
		static inline std::map<
			size_t,
			std::array<std::vector<size_t>, 3>
		> _bidxlist;

		static inline std::map<
			size_t,
			std::vector<std::vector<size_t>>
		> _uidxlist;

		static inline std::map<
			size_t,
			std::vector<std::vector<size_t>>
		> _vidxlist;

		static inline std::map<
			size_t,
			std::vector<std::vector<size_t>>
		> _widxlist;

	public:
		enum BoundaryCurveAlias { cu = 0, cv = 1, cw = 2 };

		static std::array<std::vector<size_t>, 3> BIdx(size_t order);
		static std::vector<std::vector<size_t>> UIdx(size_t order);
		static std::vector<std::vector<size_t>> VIdx(size_t order);
		static std::vector<std::vector<size_t>> WIdx(size_t order);
	};
}

#endif

#include <HashColon/CAGD/impl/Geometry_Impl.hpp>