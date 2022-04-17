#ifndef HASHCOLON_CAGD_GEOMETRY_IMPL
#define HASHCOLON_CAGD_GEOMETRY_IMPL

#include <HashColon/CAGD/Geometry.hpp>

/* Curve */
namespace HashColon::CAGD
{
	template <size_t Dim, size_t Degree>
	Curve<Dim, Degree>::Curve(
		const std::shared_ptr<CAGDBasisBase<Degree>>& iBasis)
		: GeometryBase<Dim>(PointArrT<Dim>::Zero(Dim, iBasis->Order())), basis(iBasis)
	{};

	template <size_t Dim, size_t Degree>
	Curve<Dim, Degree>::Curve(
		const std::shared_ptr<CAGDBasisBase<Degree>>& iBasis,
		const PointArrT<Dim>& featurePoints)
		: GeometryBase<Dim>(featurePoints), basis(iBasis)
	{
		// check validity of basis
		if (basis->Order() != featurePoints.cols())
			throw Exception(
				"Invalid number control points: " + std::to_string(featurePoints.cols())
				+ ", order: " + std::to_string(basis->Order()));
	}
	
	template <size_t Dim, size_t Degree>
	PointT<Dim> Curve<Dim, Degree>::PointOnCurve(Real u) const
	{
		return this->FeaturePoints * basis->Curve(u);
	}

	template <size_t Dim, size_t Degree>
	PointArrT<Dim> Curve<Dim, Degree>::PointOnCurve_PreComputed() const
	{
		return this->FeaturePoints * basis->PreComputedValue();
	}

	template <size_t Dim, size_t Degree>
	Curve<Dim, Degree> Curve<Dim, Degree>::Flip()
	{
		Curve<Dim, Degree> re = (*this);
		Eigen::MatrixXR flipMatrix = Eigen::MatrixXR::Zero(this->Basis()->Order(), this->Basis()->Order());
		for (size_t i = 0; i < this->Basis()->Order(); i++)
		{
			flipMatrix(i, this->Basis()->Order() - i - 1) = 1;
		}
		re.CP = this->CP * flipMatrix;
		return re;
	}
}

/* Surface */
namespace HashColon::CAGD
{
	template<size_t Dim>
	Surface<Dim>::Surface(
		const std::shared_ptr<SurfaceBasis>& ibasis)
		: GeometryBase<Dim>(PointArrT<Dim>::Zero(Dim, 1)), basis(ibasis)
	{};

	template<size_t Dim>
	Surface<Dim>::Surface(
		const std::shared_ptr<SurfaceBasis>& ibasis,
		const PointArrT<Dim>& featurePoints)
		: GeometryBase<Dim>(featurePoints), basis(ibasis)
	{};

	template<size_t Dim>
	Surface<Dim>::Surface(const Surface& rhs)
		: GeometryBase<Dim>(rhs.FeaturePoints), basis(rhs.basis)
	{};

	template<size_t Dim>
	PointT<Dim> Surface<Dim>::PointOnSurface(Real u, Real v) const
	{
		return CP * basis->Value({ {{u,v}} });
	}

	template<size_t Dim>
	PointArrT<Dim> Surface<Dim>::MeshOnSurface(std::vector<std::array<Real, 2>> uv, bool noThrow) const
	{
		return CP * basis->Value(uv, noThrow);
	}

	template<size_t Dim>
	PointArrT<Dim> Surface<Dim>::MeshOnSurface_PreComputed() const
	{			
		return CP * basis->PreComputedValue();
	}
}

/* QuadSurface */
namespace HashColon::CAGD
{
	template<size_t Dim, size_t uDegree, size_t vDegree>
	QuadSurface<Dim, uDegree, vDegree>::QuadSurface(
		std::shared_ptr<QuadSurfaceBasis<uDegree, vDegree>> ibasis)
		: Surface<Dim>(ibasis)
	{};

	template<size_t Dim, size_t uDegree, size_t vDegree>
	QuadSurface<Dim, uDegree, vDegree>::QuadSurface(
		std::shared_ptr<QuadSurfaceBasis<uDegree, vDegree>> ibasis,
		const PointArrT<Dim>& featurePoints)
		: Surface<Dim>(ibasis, featurePoints)
	{
		// check validity of basis
		if (ibasis->NumOfCP() != featurePoints.cols())
			throw Exception(
				"Invalid number control points: " + std::to_string(featurePoints.cols())
				+ ", u-order: " + std::to_string(ibasis->uOrder())
				+ ", v-order: " + std::to_string(ibasis->vOrder())
			);
	};

	template<size_t Dim, size_t uDegree, size_t vDegree>
	QuadSurface<Dim, uDegree, vDegree>::QuadSurface(const Surface<Dim>& rhs)
		: Surface<Dim>(rhs.basis, rhs.FeaturePoints)
	{};
}

/* TriSurface */
namespace HashColon::CAGD
{
	template <size_t Dim, size_t Degree>
	TriSurface<Dim, Degree>::TriSurface(
		const std::shared_ptr<TriSurfaceBasis<Degree>>& ibasis)
		: Surface<Dim>(ibasis, PointArrT<Dim>::Zero(Dim, ibasis->NumOfCP()))
	{};
	
	template <size_t Dim, size_t Degree>
	TriSurface<Dim, Degree>::TriSurface(
		const std::shared_ptr<TriSurfaceBasis<Degree>>& ibasis, const PointArrT<Dim>& featurePoints)
		: Surface<Dim>(ibasis, featurePoints)
	{
		// check validity of basis
		if (ibasis->NumOfCP() != featurePoints.cols())
			throw Exception(
				"Invalid number control points: " + std::to_string(featurePoints.cols())
				+ ", order: " + std::to_string(ibasis->Order())				
			);
	};

	template <size_t Dim, size_t Degree>
	TriSurface<Dim, Degree>::TriSurface(const TriSurface<Dim, Degree>& rhs)
		: Surface<Dim>(rhs.basis, rhs.FeaturePoints)
	{};

}

namespace HashColon::CAGD
{
	template <size_t Dim>
	PointArrT<Dim>& SetControlPointsByIndex(
		PointArrT<Dim>& SurfCP, const PointArrT<Dim>& CPs,
		const std::vector<size_t>& indices)
	{
		for (size_t i = 0; i < indices.size(); i++)		
			SurfCP.col(indices[i]) = CPs.col(i);
		return SurfCP;
	}
}


#endif
