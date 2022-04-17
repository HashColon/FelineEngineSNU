#ifndef HASHCOLON_CAGD_CAGDBASIS_HPP
#define HASHCOLON_CAGD_CAGDBASIS_HPP

#include <string>
#include <Eigen/Eigen>

#include <HashColon/Core/Real.hpp>
#include <HashColon/Core/Exception.hpp>

/* Curve basis (1D in domain space) */
namespace HashColon::CAGD
{
	template <size_t DegreeN>
	class CAGDBasisBase
	{
	protected:
		std::shared_ptr<Eigen::MatrixXR> preComputedValue;
		CAGDBasisBase() : preComputedValue(nullptr) {};
		HASHCOLON_CLASS_EXCEPTION_DEFINITION(CAGDBasisBase);
	public:		
		virtual size_t Order() const = 0;
		virtual size_t Degree() const = 0;

		// Basis value of index and parameter
		virtual Real Value(size_t i, Real u, bool noThrow) const = 0;
		
		/* Basis value matrix for computing points on curves
		*    POC = CP * Value(u[]) 
		*    - Row: basis indices from 0 to Order
		*    - Col: parameter u array 
		*/
		virtual Eigen::MatrixXR Value(			
			const std::vector<Real>& u, bool noThrow = true) const = 0;

		/* Compute pre-compute values 
		*   runs Value(u[], noThrow) and stores the result internally
		*/
		Eigen::MatrixXR& PreCompute(			
			const std::vector<Real>& u, bool noThrow = true);
		//Eigen::MatrixXR PreComputedValue() const { return *PreComputedValue; };
		const Eigen::MatrixXR& PreComputedValue() const;
	};

	template<size_t DegreeN>
	class BezierBasis : public CAGDBasisBase<DegreeN>
	{		
	public:
		BezierBasis() : CAGDBasisBase<DegreeN>() {};

		HASHCOLON_CLASS_EXCEPTION_DEFINITION(BezierBasis);

		virtual size_t Order() const override { return DegreeN + 1; };
		virtual size_t Degree() const override { return DegreeN; };

		virtual Real Value(size_t i, Real u, bool noThrow) const override;
		virtual Eigen::MatrixXR Value(			
			const std::vector<Real>& u, bool noThrow = true) const override;
	};
}

/* Surface basis (2D in domain space) */
namespace HashColon::CAGD
{
	class SurfaceBasis
	{
	protected:
		std::shared_ptr<Eigen::MatrixXR> preComputedValue;		
		SurfaceBasis() : preComputedValue(nullptr) {};		
	public:
		HASHCOLON_CLASS_EXCEPTION_DEFINITION(SurfaceBasis);

		// Basis value of index and parameter
		virtual Real Value(size_t i, Real u, size_t j, Real v, bool noThrow) const = 0;

		/* Basis value matrix for computing points on curves
		*    POC = CP * Value(u[])
		*    - Row: basis indices from 0 to Order (Defined by inherited class)
		*    - Col: parameter u/v array
		*/
		virtual Eigen::MatrixXR Value(			
			const std::vector<std::array<Real, 2>>& uv,
			bool nothrow = true) const = 0;		

		virtual size_t NumOfCP() const = 0;

		/* Compute pre-compute values
		*   runs Value(u[], noThrow) and stores the result internally
		*/
		Eigen::MatrixXR& PreCompute(			
			const std::vector<std::array<Real, 2>>& uv,
			bool noThrow = true);
		//Eigen::MatrixXR PreComputedValue() const { return *PreComputedValue; };
		const Eigen::MatrixXR& PreComputedValue() const;
	};

	template <size_t uDegree, size_t vDegree>
	class QuadSurfaceBasis : public SurfaceBasis
	{
	protected:
		std::shared_ptr<CAGDBasisBase<uDegree>> uBasis;
		std::shared_ptr<CAGDBasisBase<vDegree>> vBasis;
	public:
		QuadSurfaceBasis(
			const std::shared_ptr<CAGDBasisBase<uDegree>>& ubasis,
			const std::shared_ptr<CAGDBasisBase<vDegree>>& vbasis);

		HASHCOLON_CLASS_EXCEPTION_DEFINITION(QuadSurfaceBasis);

		size_t uOrder() const { return uBasis->Order(); };
		size_t vOrder() const { return vBasis->Order(); };
		size_t uDegree() const { return uBasis->Degree(); };
		size_t vDegree() const { return vBasis->Degree(); };

		virtual Real Value(size_t ui, Real u, size_t vi, Real v, bool noThrow = true) const override;
		virtual Eigen::MatrixXR Value(
			const std::vector<std::array<Real, 2>>& uv,
			bool noThrow = true) const override;
		virtual size_t NumOfCP() const override;
		
		
	};

	template <size_t DegreeN>
	class TriSurfaceBasis : public SurfaceBasis
	{
	protected:
		TriSurfaceBasis() : SurfaceBasis() {};
	public:
		virtual size_t Degree() const { return DegreeN; };
		virtual size_t Order() const = 0;		
	};
		
	template <size_t DegreeN>
	class TriBezierSurfaceBasis : public TriSurfaceBasis<DegreeN>
	{
	public:
		TriBezierSurfaceBasis() : TriSurfaceBasis<DegreeN>() {};

		HASHCOLON_CLASS_EXCEPTION_DEFINITION(TriBezierSurfaceBasis);				

		virtual size_t Order() const override { return DegreeN + 1; };
		virtual size_t Degree() const override { return DegreeN; };

		virtual Real Value(size_t i, Real u, size_t j, Real v, bool noThrow) const override;
		virtual Eigen::MatrixXR Value(
			const std::vector<std::array<Real, 2>>& uv,
			bool nothrow = true) const override;
		virtual size_t NumOfCP() const override;
	};
}

namespace HashColon::CAGD::_helper
{
	Eigen::VectorXR FlattenMatrixColwise(Eigen::MatrixXR A);
	Eigen::VectorXR FlattenMatrixRowwise(Eigen::MatrixXR A);
}

#endif

#include <HashColon/CAGD/impl/CAGDBasis_Impl.hpp>