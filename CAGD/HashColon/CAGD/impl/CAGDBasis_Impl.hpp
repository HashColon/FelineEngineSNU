#ifndef HASHCOLON_CAGD_CAGDBASIS_IMPL
#define HASHCOLON_CAGD_CAGDBASIS_IMPL

#include <memory>
#include <mutex>
#include <vector>
#include <cmath>
#include <Eigen/Eigen>
#include <HashColon/CAGD/CAGDBasis.hpp>

/* helpers */
namespace HashColon::CAGD
{
	class CAGDBasisPreCalculationHelper
	{
	private:
		static inline std::shared_ptr<CAGDBasisPreCalculationHelper> _instance = nullptr;
		static inline std::once_flag _onlyOne;
		static inline std::mutex _m;

		std::vector<size_t> _factorial;

		// Constructor & destructor	
		CAGDBasisPreCalculationHelper() : _factorial(1, 1) {};
		CAGDBasisPreCalculationHelper(const CAGDBasisPreCalculationHelper& rs) { _instance = rs._instance; };
		CAGDBasisPreCalculationHelper& operator=(const CAGDBasisPreCalculationHelper& rs)
		{
			if (this != &rs) { _instance = rs._instance; }
			return *this;
		};

	public:
		static CAGDBasisPreCalculationHelper& GetInstance()
		{
			using namespace std;
			call_once(_onlyOne,
				[](size_t idx) {
					_instance.reset(new CAGDBasisPreCalculationHelper());
				}, 0);
			return *_instance;
		};

		size_t Factorial(size_t k)
		{
			using namespace std;
			if (_factorial.size() <= k + 1)
			{
				lock_guard<mutex> _lg(_m);
				for (; _factorial.size() <= k + 1; )
				{
					size_t i = _factorial.size();
					_factorial.push_back(_factorial.back() * i);
				}
			}
			return _factorial[k];
		};

		size_t N_Choose_K(size_t n, size_t k)
		{
			return Factorial(n) / Factorial(k) / Factorial(n - k);			
		}
	};
};

/* CAGDBasisBase */
namespace HashColon::CAGD
{
	template<size_t D>
	Eigen::MatrixXR& CAGDBasisBase<D>::PreCompute(
		const std::vector<Real>& u, bool noThrow)
	{
		preComputedValue = std::make_shared<Eigen::MatrixXR>(Value(u, noThrow));
		return *preComputedValue;
	}

	template<size_t D>
	const Eigen::MatrixXR& CAGDBasisBase<D>::PreComputedValue() const 
	{
		if (!preComputedValue) throw Exception("Pre-computed value not computed yet");
		else return *preComputedValue; 
	}
}

/* BezierBasis */
namespace HashColon::CAGD
{
	template<size_t D>
	Real BezierBasis<D>::Value(size_t i, Real u, bool noThrow) const 
	{
		using namespace std;

		if ((i >= Order()) || (u < 0.0 && u > 1.0))
		{
			if (noThrow) return 0.0;
			else
				throw Exception(
					"Invalid input for Bezier basis: D=" + to_string(D) + ", i=" + to_string(i) + " u=" + to_string(u));
		}

		return
			(Real)CAGDBasisPreCalculationHelper::GetInstance().N_Choose_K(D, i) *
			pow(u, (Real)i) * pow(1.0 - u, (Real)(D - i));
	}

	template <size_t D>
	Eigen::MatrixXR BezierBasis<D>::Value(
		const std::vector<Real>& u, bool noThrow) const
	{
		Eigen::MatrixXR re(Order(), u.size());
		for (size_t i = 0; i < Order(); i++)
		{
			for (size_t j = 0; j < u.size(); j++)
				re(i, j) = Value(i, u[j], noThrow);
		}
		return re;
	}
}

/* SurfaceBasis is defined in cpp file */

/* QuadSurfaceBasis */
namespace HashColon::CAGD
{
	template <size_t uD, size_t vD>
	QuadSurfaceBasis<uD, vD>::QuadSurfaceBasis(
		const std::shared_ptr<CAGDBasisBase<uD>>& ubasis,
		const std::shared_ptr<CAGDBasisBase<vD>>& vbasis)
		: SurfaceBasis(), uBasis(ubasis), vBasis(vbasis)
	{}

	template <size_t uD, size_t vD>
	Real QuadSurfaceBasis<uD, vD>::Value(size_t i, Real u, size_t j, Real v, bool noThrow) const 
	{
		const CAGDBasisBase<uD>& ub = (*uBasis);
		const CAGDBasisBase<vD>& vb = (*vBasis);
		return ub.Value(i, u, noThrow) * vb.Value(j, v, noThrow);
	}

	template <size_t uD, size_t vD>
	Eigen::MatrixXR QuadSurfaceBasis<uD, vD>::Value(
		const std::vector<std::array<Real, 2>>& uv,
		bool noThrow) const 
	{
		using namespace Eigen;
		using namespace std;

		MatrixXR re(NumOfCP(), uv.size());

		for (size_t j = 0; j < uv.size(); j++)
		{
			vector<Real> u{ uv[j][0] };
			vector<Real> v{ uv[j][1] };

			// compute[ui * vj]s
			MatrixXR basisMesh = uBasis->Value(u, noThrow) * vBasis->Value(v, noThrow).transpose();

			// convert squared [ ui * vi ] to linearlized vector (left part of the below eq.)
			re.col(j) = _helper::FlattenMatrixRowwise(basisMesh);
		}		
		return re;
	}

	template <size_t uD, size_t vD>
	size_t QuadSurfaceBasis<uD, vD>::NumOfCP() const
	{
		return uBasis->Order() * vBasis->Order();
	}
}

/* TriBezierSurfaceBasis */
namespace HashColon::CAGD
{	
	template <size_t D>
	Real TriBezierSurfaceBasis<D>::Value(size_t i, Real u, size_t j, Real v, bool noThrow) const
	{	
		if (!noThrow) 
		{
			assert(u + v <= 1.0);
			if (u + v > 1.0) throw Exception("u + v exceeds 1.0.");
			assert(i + j <= Degree());
			if (i + j > Degree()) throw Exception("Sum of indices exceeds order");
		}

		Real re = 1.0;
		size_t k = Degree() - i - j;
		Real w = 1.0 - u - v;
				
		re *= (CAGDBasisPreCalculationHelper::GetInstance().Factorial(Degree()) 
			/ CAGDBasisPreCalculationHelper::GetInstance().Factorial(i) 
			/ CAGDBasisPreCalculationHelper::GetInstance().Factorial(j) 
			/ CAGDBasisPreCalculationHelper::GetInstance().Factorial(k));

		for (size_t p = 0; p < i; p++) re *= u;
		for (size_t p = 0; p < j; p++) re *= v;
		for (size_t p = 0; p < k; p++) re *= w;

		return re;
	}

	template <size_t D>
	Eigen::MatrixXR TriBezierSurfaceBasis<D>::Value(
		const std::vector<std::array<Real, 2>>& uv,
		bool noThrow) const
	{
		Eigen::MatrixXR re(NumOfCP(), uv.size());
		for (size_t c = 0; c < uv.size(); c++)
		{
			size_t cnt = 0;
			const Real& u = uv[c][0];
			const Real& v = uv[c][1];

			for (size_t i = 0; i < Order(); i++)
			{
				for (size_t j = 0; j < Order() - i; j++) {
					re(cnt, c) = Value(i, u, j, v, noThrow);
					cnt++;
				}
			}
		}
		return re;
	}

	template <size_t D>
	size_t TriBezierSurfaceBasis<D>::NumOfCP() const
	{
		return (D + 1) * (D + 2) / 2;
	}
}


#endif