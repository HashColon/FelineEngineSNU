// HashColon config
#include <HashColon/HashColon_config.h>
// std libraries
#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>
#include <queue>
#include <vector>
// dependant external libraries
#include <Eigen/Eigen>
#include <flann/flann.hpp>
// modified external libraries
#include <CLI11_modified/CLI11_extended.hpp>
#include <dkm_modified/dkm_parallel.hpp>
// HashColon libraries
#include <HashColon/Exception.hpp>
#include <HashColon/Log.hpp>
#include <HashColon/Real.hpp>
#include <HashColon/SingletonCLI.hpp>
#include <HashColon/Statistics.hpp>
// header file for this source file
#include <HashColon/Clustering.hpp>

using namespace std;
using namespace Eigen;
using namespace HashColon;
using namespace HashColon::LogUtils;
using namespace HashColon::Statistics;

using PointType = std::vector<HashColon::Real>;

// DBSCAN
namespace HashColon::Clustering
{
	void DBSCAN::Initialize(const string configFilePath)
	{
		CLI::App* cli = SingletonCLI::GetInstance().GetCLI("Clustering.DBSCAN");

		if (!configFilePath.empty())
		{
			SingletonCLI::GetInstance().AddConfigFile(configFilePath);
		}

		cli->add_option("--minPts", _cDefault.minPts, "minPts value. If a point is a core point, at least minPts number of points should be in range of epsilon.");
		cli->add_option("--DbscanEpsilon", _cDefault.DbscanEpsilon, "Epsilon value. Range value for cheking density");
	}

	void DBSCAN::TrainModel(
		const typename ClusteringBase<PointType>::DataListType& iTrainingData,
		typename ClusteringBase<PointType>::LabelsPtr oLabels,
		typename ClusteringBase<PointType>::ProbListPtr oProbabilities)
	{
		CommonLogger logger;
		{
			lock_guard<mutex> _lg(CommonLogger::_mutex);
			logger.Log({ { Tag::lvl, 3} }) << "DBSCAN: Started. " << endl;
		}

		// Assertion	
		// assert at least 1 data is given for clustering
		assert(_c.minPts > 0);
		assert(iTrainingData.size() >= _c.minPts);
		assert(!ClusteringBase<PointType>::isTrained);
		assert(oLabels != nullptr);
		if (ClusteringBase<PointType>::isTrained) throw Exception("Already trained.");
		if (oLabels == nullptr) throw Exception("DBSCAN needs cluster label output for input.");

		// compute distance matrix & neighbors
		vector<vector<size_t>> neighbors = GetNeighbors(iTrainingData);

		{
			lock_guard<mutex> _lg(CommonLogger::_mutex);
			logger.Log({ { Tag::lvl, 3 } }) << "DBSCAN: Neighbor computation finished. " << endl;

			stringstream ss;
			for (const auto& debugout1 : neighbors)
			{
				for (const auto& debugout2 : debugout1)
				{
					ss << debugout2 << "\t";
				}
				ss << "\n";
			}
			logger.Debug({ {__CODEINFO_TAGS__} }) << "\n" <<
				"Neighbor lists:\n" << ss.str() << endl;
		}

		// Initialize cluster state
		oLabels->clear();
		oLabels->resize(iTrainingData.size());

		// set meanings for the cluster idx
		// unclassified: 0, noise: 1, clustered: 2~		
		size_t clusterIdx = 2;

		for (size_t i = 0; i < iTrainingData.size(); i++)
		{
			// if the point is classfied already, continue;
			if (oLabels->at(i) != unclassified) continue;

			// if the point is core point: has more then minPts in radius epsilon,
			// run bfs
			if (neighbors[i].size() >= _c.minPts)
			{
				DbscanBfs(i, clusterIdx, neighbors, (*oLabels));
				clusterIdx++;
			}
			// else this point is noise
			else
			{
				oLabels->at(i) = noise;
			}
		}

		// as all classification is finished,
			// remove unclassified flags 
		for (size_t i = 0; i < iTrainingData.size(); i++)
		{
			assert(oLabels->at(i) > 0);
			oLabels->at(i)--;
		}

		{
			lock_guard<mutex> _lg(CommonLogger::_mutex);
			logger.Log({ { Tag::lvl, 3} }) << "DBSCAN: Finished. " << endl;
		}
		_numOfClusters = clusterIdx;
		ClusteringBase<PointType>::isTrained = true;
	}

	vector<vector<size_t>> DBSCAN::GetNeighbors(
		const typename ClusteringBase<vector<Real>>::DataListType& iTrainingData) const
	{
		vector<vector<size_t>> re;

		// use flann to build neighbor information
		size_t dataN = iTrainingData.size();
		size_t dataD = iTrainingData[0].size();
		vector<Real> flannDataset_; flannDataset_.resize(dataN * dataD);
		flann::Matrix<Real> flannDataset(flannDataset_.data(), dataN, dataD);

		// copy sample data to flann Matrix
		#pragma omp parallel for
		for (size_t i = 0; i < dataN; i++)
		{
			for (size_t j = 0; j < dataD; j++)
				flannDataset[i][j] = iTrainingData[i][j];
		}
		// build kd-tree using flann using 2-norm
		flann::Index<flann::L2<double>> kdtree(flannDataset, flann::KDTreeIndexParams());
		kdtree.buildIndex();

		// results from flann
		vector<vector<int>> idxs;
		vector<vector<double>> dists;

		kdtree.radiusSearch(
			flannDataset,
			idxs, dists,
			(float)_c.DbscanEpsilon,
			flann::SearchParams()
		);

		// we don't need distance info 
		dists.clear();

		// convert vectr<vector<int>> to vectr<vector<size_t>>
		re.resize(idxs.size());
		#pragma omp parallel for
		for (size_t i = 0; i < idxs.size(); i++)
		{
			re[i].resize(idxs[i].size());
			for (size_t j = 0; j < idxs[i].size(); i++)
			{
				re[i][j] = idxs[i][j];
			}
		}

		return re;
	}

	void DBSCAN::DbscanBfs(size_t initP, size_t clusterIdx,
		vector<vector<size_t>>& neighbors, vector<size_t>& labels) const
	{
		using namespace std;
		assert(neighbors[initP].size() >= _c.minPts);
		assert(labels[initP] == unclassified);
		queue<size_t> q; q.push(initP);
		labels[initP] = clusterIdx;
		do
		{
			// pop a point from queue, add it to current cluster
			size_t p = q.front(); q.pop();
			//labels[p] = clusterIdx;

			// if current point has sufficient neighbors,
			if (neighbors[p].size() >= _c.minPts)
			{
				// push neighbors to queue 
				for (auto& n : neighbors[p]) 
				{
					if(labels[n] != clusterIdx)
					{
						q.push(n);
						labels[n] = clusterIdx;
					}					
				}
			}
		} while (!q.empty());
	}
}

// Kmenas
namespace HashColon::Clustering
{
	void Kmeans::Initialize(const string configFilePath)
	{
		CLI::App* cli = SingletonCLI::GetInstance().GetCLI("Clustering.Kmeans");

		if (!configFilePath.empty())
		{
			SingletonCLI::GetInstance().AddConfigFile(configFilePath);
		}

		cli->add_option("--k", _cDefault.k, "K value for K-means clustering");
		cli->add_option("--kmeansEpsilon", _cDefault.kmeansEpsilon, "Difference criteria for K-means clustering");
		cli->add_option("--kmeansIteration", _cDefault.kmeansIteration, "Max iteration number for K-means clustering");
	}

	void Kmeans::TrainModel(
		const typename ClusteringBase<PointType>::DataListType& iTrainingData,
		typename ClusteringBase<PointType>::LabelsPtr oLabels,
		typename ClusteringBase<PointType>::ProbListPtr oProbabilities)
	{
		CommonLogger logger;
		logger.Log({ {Tag::lvl, 3} }) << "Kmeans: Started. " << endl;

		// K-means clustering using dkm
		dkm::clustering_parameters<Real> dkm_params((unsigned int)_c.k);
		dkm_params.set_max_iteration(_c.kmeansIteration);
		dkm_params.set_min_delta(_c.kmeansEpsilon);

		vector<vector<Real>> means;
		auto [dkm_means, dkm_labels] = dkm::kmeans_lloyd_parallel(iTrainingData, dkm_params);
		//auto [dkm_means, dkm_labels] = dkm::kmeans_lloyd(samples, dkm_params);
		logger.Log({ {Tag::lvl, 3} }) << "Kmeans: K-means clustering finished." << endl;

		// save clustering results in oLabels
		oLabels->clear();
		for (auto& label : dkm_labels)
			oLabels->push_back(label);

		// compute likelihood: Not Implemented

		logger.Log({ {Tag::lvl, 3} }) << "Kmeans: Finished." << endl;

		// Training finished.
		ClusteringBase<PointType>::isTrained = true;
	}

	size_t Kmeans::GetClusterOf(
		const PointType& iTestValue,
		typename ClusteringBase<PointType>::ProbPtr oProbabilities)
	{
		throw NotImplementedException;
	}


}

// Evaluation functions for clustering
namespace HashColon::Clustering
{
	// Get Min/Max/Mean/Median/Standard deviation from a given cluster
	SimpleStatisticsAnalysisResults DistanceAnalysis(const vector<Real>& DistancesInCluster)
	{
		SimpleStatisticsAnalysisResults re;
		re.NewAll();
		/*(*re.Min) = Min(DistancesInCluster);
		(*re.Max) = Max(DistancesInCluster);
		(*re.Mean) = Mean(DistancesInCluster);
		(*re.Median) = Median(DistancesInCluster);
		(*re.Variance) = Variance(DistancesInCluster);
		(*re.StdDeviation) = StdDeviation(DistancesInCluster);*/

		return re;
	}

	DistancesAnalysisResults DistanceAnalysis(const vector<vector<Real>>& DistancesInClusters)
	{
		DistancesAnalysisResults re;
		vector<Real> wholeDist;
		re.ClusterResults.resize(DistancesInClusters.size());
		for (size_t i = 0; i < DistancesInClusters.size(); i++)
		{
			re.ClusterResults[i] = DistanceAnalysis(DistancesInClusters[i]);
			copy(DistancesInClusters[i].begin(), DistancesInClusters[i].end(), back_inserter(wholeDist));
		}
		re.TotalResults = DistanceAnalysis(wholeDist);
		return re;
	}

	vector<vector<Real>> SortedDistanceGraph(
		const vector<size_t>& clusterResult, const MatrixXR& DistanceMatrix)
	{
		// DistanceMatrix should be square matrix & size of cluster result should be equal to size of distance matrix
		assert((long int)clusterResult.size() == DistanceMatrix.rows());
		assert(DistanceMatrix.rows() == DistanceMatrix.cols());
		const auto N = clusterResult.size();

		// count number of clusters
		size_t clusterNo = (*max_element(clusterResult.begin(), clusterResult.end())) + 1;

		// Distances grouped by each cluster
		vector<vector<Real>> DistancesInClusters;
		DistancesInClusters.resize(clusterNo);
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = i + 1; j < N; j++)
			{
				if (clusterResult[i] == clusterResult[j])
				{
					DistancesInClusters[clusterResult[i]].push_back(DistanceMatrix(i, j));
				}
			}
		}

		// sort results
		for (auto& distances : DistancesInClusters)
		{
			sort(distances.begin(), distances.end());
		}

		return DistancesInClusters;
	}

	vector<vector<Real>> SortedDistanceGraph(
		const vector<size_t>& clusterResult, const MatrixXR& DistanceMatrix,
		DistancesAnalysisResults& additionals)
	{
		// DistanceMatrix should be square matrix & size of cluster result should be equal to size of distance matrix
		assert((long int)clusterResult.size() == DistanceMatrix.rows());
		assert(DistanceMatrix.rows() == DistanceMatrix.cols());
		const auto N = clusterResult.size();

		// count number of clusters
		size_t clusterNo = (*max_element(clusterResult.begin(), clusterResult.end())) + 1;

		// Distances grouped by each cluster
		vector<vector<Real>> DistancesInClusters;
		DistancesInClusters.resize(clusterNo);
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = i + 1; j < N; j++)
			{
				if (clusterResult[i] == clusterResult[j])
				{
					DistancesInClusters[clusterResult[i]].push_back(DistanceMatrix(i, j));
				}
			}
		}

		// sort results
		for (auto& distances : DistancesInClusters)
		{
			sort(distances.begin(), distances.end());
		}

		// compute distance analysis
		additionals = DistanceAnalysis(DistancesInClusters);

		return DistancesInClusters;
	}

	vector<size_t> PseudoMedian(
		const vector<size_t>& clusterResult, const Eigen::MatrixXR& DistanceMatrix)
	{
		// DistanceMatrix should be square matrix & size of cluster result should be equal to size of distance matrix
		assert((long int)clusterResult.size() == DistanceMatrix.rows());
		assert(DistanceMatrix.rows() == DistanceMatrix.cols());
		const auto N = clusterResult.size();

		size_t clusterNo = (*max_element(clusterResult.begin(), clusterResult.end())) + 1;
		vector<size_t> medians(clusterNo);
		vector<Real> sum(clusterNo, numeric_limits<Real>::max());

		for (size_t i = 0; i < N; i++)
		{
			Real tmpsum = 0;
			for (size_t j = 0; j < N; j++)
			{
				if (i == j) continue;
				else if (clusterResult[i] == clusterResult[j])
				{
					tmpsum += DistanceMatrix((int)i, (int)j);
				}
			}
			if (tmpsum < sum[clusterResult[i]])
			{
				sum[clusterResult[i]] = tmpsum;
				medians[clusterResult[i]] = i;
			}
		}

		return medians;
	}
	
	vector<Real> PseudoDaviesBouldin(
		const vector<size_t>& clusterResult, const Eigen::MatrixXR& DistanceMatrix,
		vector<size_t>& pseudoMedian)
	{
		// DistanceMatrix should be square matrix & size of cluster result should be equal to size of distance matrix
		assert((long int)clusterResult.size() == DistanceMatrix.rows());
		assert(DistanceMatrix.rows() == DistanceMatrix.cols());
		const auto N = clusterResult.size();

		size_t clusterNo = (*max_element(clusterResult.begin(), clusterResult.end())) + 1;
		vector<Real> re(clusterNo, 0.0);
		vector<Real> cnt(clusterNo, 0.0);

		pseudoMedian = PseudoMedian(clusterResult, DistanceMatrix);

		for (size_t i = 0; i < N; i++)
		{
			Real d = DistanceMatrix(i, pseudoMedian[clusterResult[i]]);
			re[clusterResult[i]] += d * d;
			cnt[clusterResult[i]] += 1.0;
		}

		for (size_t i = 0; i < clusterNo; i++)
		{
			if (cnt[i] == 0)
				re[i] = 0;
			else
				re[i] = sqrt(re[i] / cnt[i]);
		}

		return re;
	}

	vector<Real> PseudoDaviesBouldin(
		const vector<size_t>& clusterResult, const Eigen::MatrixXR& DistanceMatrix)
	{
		vector<size_t> pseudoMedian;
		return PseudoDaviesBouldin(clusterResult, DistanceMatrix, pseudoMedian);
	}

	Real Silhouette(size_t itemIdx, const vector<size_t>& clusterResult, const Eigen::MatrixXR& DistanceMatrix)
	{
		// DistanceMatrix should be square matrix & size of cluster result should be equal to size of distance matrix
		assert((long int)clusterResult.size() == DistanceMatrix.rows());
		assert(DistanceMatrix.rows() == DistanceMatrix.cols());
		const auto N = clusterResult.size();

		// compute fitness to each clusters
		size_t clusterNo = (*max_element(clusterResult.begin(), clusterResult.end())) + 1;
		vector<Real> A(clusterNo, 0.0);
		vector<Real> cnt(clusterNo, 0.0);
		for (size_t i = 0; i < N; i++)
		{
			if (i == itemIdx) continue;
			A[clusterResult[i]] += DistanceMatrix(itemIdx, i);
			cnt[clusterResult[i]] += 1.0;
		}

		// if cnt of the cluster including itemIdx is equal or lesser than 1.0, return 0.0;
		if (cnt[clusterResult[itemIdx]] <= 1.0) return 0.0;

		// compute a(itemIdx) && b(itemIdx)
		Real a_item = 0.0;
		Real b_item = numeric_limits<Real>::max();
		for (size_t i = 0; i < A.size(); i++)
		{
			Real tmpVal = A[i] / (cnt[i]);
			if (i == clusterResult[itemIdx])
			{
				a_item = tmpVal;
			}
			else if (tmpVal < b_item)
			{
				b_item = tmpVal;
			}
		}

		Real div = a_item > b_item ? a_item : b_item;
		return (b_item - a_item) / div;
	}

	vector<Real> Silhouette(const vector<size_t>& clusterResult, const Eigen::MatrixXR& DistanceMatrix)
	{
		// DistanceMatrix should be square matrix & size of cluster result should be equal to size of distance matrix
		assert((long int)clusterResult.size() == DistanceMatrix.rows());
		assert(DistanceMatrix.rows() == DistanceMatrix.cols());
		const auto N = clusterResult.size();

		vector<Real> re(N);

#pragma omp parallel for
		for (size_t i = 0; i < N; i++)
		{
			re[i] = Silhouette(i, clusterResult, DistanceMatrix);
		}

		return re;
	}
}

namespace HashColon::Clustering::Evaluation
{
	void DoAssertion(const vector<size_t>& Labels, const Eigen::MatrixXR& DistanceMatrix)
	{
		assert(DistanceMatrix.rows() == DistanceMatrix.cols());
		assert(DistanceMatrix.rows() == (long int)Labels.size());

		if (DistanceMatrix.rows() != DistanceMatrix.cols())
			throw EvaluationException("Distance matrix is not square matrix.");
		if (DistanceMatrix.rows() != (long int)Labels.size())
			throw EvaluationException("Distance matrix size does not match with labels size");
	}

	size_t GetClusterNo(const vector<size_t>& Labels)
	{
		return (*max_element(Labels.begin(), Labels.end())) + 1;
	}

	vector<bool> EvalOptions::GetExcludeItemFlags(const size_t itemN)
	{
		vector<bool> re(itemN, false);
		for (size_t idx : ExcludeItemIndices) re[idx] = true;
		return re;		
	}

	vector<bool> EvalOptions::GetExcludeClusterFlags(const size_t clusterN)
	{
		vector<bool> re(clusterN, false);
		for (size_t idx : ExcludeClusterIndices) re[idx] = true;
		return re;
	}

	size_t GetMedian(const size_t clusterIdx, 
		const vector<size_t>& Labels, const MatrixXR& DistanceMatrix, 
		EvalOptions opt)
	{
		DoAssertion(Labels, DistanceMatrix);		
		const size_t N = Labels.size();
		vector<bool> excludeItemFlag = opt.GetExcludeItemFlags(N);

		size_t medianIdx;
		Real minSum = numeric_limits<Real>::max();

		for (size_t i = 0; i < N; i++)
		{
			if (Labels[i] != clusterIdx) continue;
			if (excludeItemFlag[i]) continue;

			Real tmpsum = 0;
			for (size_t j = 0; j < N; j++)
			{
				if (i == j) continue;
				if (excludeItemFlag[j]) continue;

				else if (clusterIdx == Labels[j])
				{
					tmpsum += DistanceMatrix((int)i, (int)j);
				}
			}
			if (tmpsum < minSum)
			{
				minSum = tmpsum;
				medianIdx = i;
			}
		}

		return medianIdx;
	}

	vector<size_t> GetMedians(
		const vector<size_t>& Labels, const Eigen::MatrixXR& DistanceMatrix,
		EvalOptions opt)
	{
		DoAssertion(Labels, DistanceMatrix);
		size_t clusterNo = GetClusterNo(Labels);
		vector<bool> excludeClusterFlag = opt.GetExcludeClusterFlags(clusterNo);

		vector<size_t> medians(clusterNo, numeric_limits<size_t>::max());

		#pragma omp parallel for
		for (size_t i = 0; i < clusterNo; i++)
		{
			if (excludeClusterFlag[i]) continue;
			medians[i] = GetMedian(i, Labels, DistanceMatrix, opt);
		}

		return medians;
	}

	Real GetSilhouette(size_t itemIdx,
		const vector<size_t>& Labels, const Eigen::MatrixXR& DistanceMatrix,
		EvalOptions opt)
	{
		// DistanceMatrix should be square matrix & size of cluster result should be equal to size of distance matrix
		DoAssertion(Labels, DistanceMatrix);
		size_t clusterNo = GetClusterNo(Labels);
		const auto N = Labels.size();	
		vector<bool> excludeItemFlag = opt.GetExcludeItemFlags(N);
		vector<bool> excludeClusterFlag = opt.GetExcludeClusterFlags(clusterNo);
		
		// compute fitness to each clusters		
		vector<Real> A(clusterNo, 0.0);
		vector<Real> cnt(clusterNo, 0.0);
		for (size_t i = 0; i < N; i++)
		{
			if (i == itemIdx) continue;
			if (excludeItemFlag[i]) continue;
			A[Labels[i]] += DistanceMatrix(itemIdx, i);
			cnt[Labels[i]] += 1.0;
		}

		// if cnt of the cluster including itemIdx is equal or lesser than 1.0, return 0.0;
		if (cnt[Labels[itemIdx]] <= 1.0) return 0.0;

		// compute a(itemIdx) && b(itemIdx)
		Real a_item = 0.0;
		Real b_item = numeric_limits<Real>::max();
		for (size_t i = 0; i < A.size(); i++)
		{
			Real tmpVal = A[i] / (cnt[i]);

			if (excludeClusterFlag[i]) continue;
			else if (i == Labels[itemIdx]) a_item = tmpVal;			
			else if (tmpVal < b_item) b_item = tmpVal;
		}

		Real div = a_item > b_item ? a_item : b_item;
		return (b_item - a_item) / div;
	}

	vector<Real> GetSilhouettes(
		const vector<size_t>& Labels, const Eigen::MatrixXR& DistanceMatrix,
		EvalOptions opt)
	{
		// DistanceMatrix should be square matrix & size of cluster result should be equal to size of distance matrix
		DoAssertion(Labels, DistanceMatrix);
		size_t clusterNo = GetClusterNo(Labels);
		const auto N = Labels.size();
		vector<bool> excludeItemFlag = opt.GetExcludeItemFlags(N);
		vector<bool> excludeClusterFlag = opt.GetExcludeClusterFlags(clusterNo);

		vector<Real> re(N);
		#pragma omp parallel for
		for (size_t i = 0; i < N; i++)
		{
			if (excludeItemFlag[i]) continue;
			else if (excludeClusterFlag[Labels[i]]) continue;
			else re[i] = GetSilhouette(i, Labels, DistanceMatrix, opt);
		}

		return re;
	}

	Real GetSilhouetteCoefficient(
		const vector<size_t>& Labels, const MatrixXR& DistanceMatrix,
		const shared_ptr<vector<Real>>& Silhouettes,
		EvalOptions opt)
	{
		// DistanceMatrix should be square matrix & size of cluster result should be equal to size of distance matrix
		DoAssertion(Labels, DistanceMatrix);
		size_t clusterNo = GetClusterNo(Labels);
		const auto N = Labels.size();
		vector<bool> excludeItemFlag = opt.GetExcludeItemFlags(N);
		vector<bool> excludeClusterFlag = opt.GetExcludeClusterFlags(clusterNo);

		vector<Real> _silhouettes;
		vector<Real> _meanSilhouettes(clusterNo, 0.0);
		vector<size_t> _clusterSize(clusterNo, 0);
		Real resSC = 2.0;

		// if Silhouettes are given, use it
		if (Silhouettes && Silhouettes->size() == N) _silhouettes = (*Silhouettes);
		// else compute it
		else
		{
			_silhouettes = GetSilhouettes(Labels, DistanceMatrix, opt);
			if (Silhouettes) (*Silhouettes) = _silhouettes;
		}

		// compute mean silhouettes
		for (size_t i = 0; i < N; i++)
		{
			if (excludeItemFlag[i]) continue;
			if (excludeClusterFlag[Labels[i]]) continue;

			_meanSilhouettes[Labels[i]] += _silhouettes[i];
			_clusterSize[Labels[i]]++;
		}		
		for (size_t c = 0; c < clusterNo; c++)
		{
			if (excludeClusterFlag[c]) continue;
			_meanSilhouettes[c] /= (Real)_clusterSize[c];
		}

		// get minimum mean silhouette
		for (size_t c = 0; c < clusterNo; c++)
		{
			if (excludeClusterFlag[c]) continue;
			else if (_meanSilhouettes[c] < resSC) resSC = _meanSilhouettes[c];
		}

		return resSC;
	}

	Real DaviesBouldinIndex_usingMedian(
		const vector<size_t>& Labels, const MatrixXR& DistanceMatrix,
		const shared_ptr<vector<size_t>> Medians, EvalOptions opt, const Real nearZero)
	{
		DoAssertion(Labels, DistanceMatrix);
		size_t clusterNo = GetClusterNo(Labels);
		const auto N = Labels.size();
		vector<bool> excludeItemFlag = opt.GetExcludeItemFlags(N);
		vector<bool> excludeClusterFlag = opt.GetExcludeClusterFlags(clusterNo);

		Real DBI = 0;

		// get medians. if median is given, use the given. else, compute medians
		vector<size_t> _medians = (Medians && Medians->size() == clusterNo) ?
			(*Medians) :
			GetMedians(Labels, DistanceMatrix, opt);

		// compute S_i
		vector<Real> S_i(clusterNo, 0.0);
		vector<size_t> cnt_i(clusterNo, 0);
		for (size_t i = 0; i < N; i++)
		{
			const size_t& c = Labels[i];
			if (excludeItemFlag[i]) continue;
			if (excludeClusterFlag[c]) continue;
			
			S_i[c] += DistanceMatrix(_medians[c], i);
			cnt_i[c]++;
		}
		for (size_t c = 0; c < clusterNo; c++)
		{
			if (excludeClusterFlag[c]) continue;
			S_i[c] /= (Real)cnt_i[c];
		}
			

		// compute M_ij
		vector<vector<Real>> M_ij(clusterNo, vector<Real>(clusterNo, 0.0));
		for (size_t i = 0; i < clusterNo; i++)
		{
			if (excludeClusterFlag[i]) continue;
			for (size_t j = i+1; j < clusterNo; j++)
			{
				if (excludeClusterFlag[j]) continue;
				auto d = DistanceMatrix(_medians[i], _medians[j]);
				M_ij[i][j] = M_ij[j][i] = d <= 0.0 ? nearZero : d;
			}
		}

		// compute R_ij, D_i
		for (size_t i = 0; i < clusterNo; i++)
		{
			if (excludeClusterFlag[i]) continue;
			Real D_i = 0;
			for (size_t j = 0; j < clusterNo; j++)
			{
				if (excludeClusterFlag[j]) continue;
				else if (i == j) continue;
				Real R_ij = (S_i[i] + S_i[j]) / M_ij[i][j];

				if (R_ij > D_i) D_i = R_ij;
			}
			DBI += D_i;
		}
		DBI /= (Real)clusterNo;
		return DBI;
	}
}