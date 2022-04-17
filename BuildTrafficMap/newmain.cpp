#include <HashColon/HashColon_config.h>
#include <mutex>
#include <string>
#include <vector>
#include <HashColon/Clustering.hpp>
#include <HashColon/Exception.hpp>
#include <HashColon/Log.hpp>
#include <HashColon/Real.hpp>
#include <HashColon/SingletonCLI.hpp>
#include <HashColon/ThreadPool.hpp>
#include <HashColon/Feline/FelineJsonIO.hpp>
#include <HashColon/Feline/GeoValues.hpp>

using namespace std;
using namespace HashColon;
using namespace HashColon::Feline;
using namespace HashColon::LogUtils;

struct _Params
{
	int numOfTreads = 1;

	vector<string> inputTrajFilepath;

	// flags
	bool PreClustering_byDeptArri;
	bool PreClustering_byShipType;
};

_Params _c;

size_t GetClusterNum(const vector<size_t>& labels) { return (*max_element(labels.begin(), labels.end())) + 1; }

void PreClustering_byShipType(
	const vector<AisTrajectory<>>& iTrajClusters,
	vector<size_t>& oLabels)
{	
	// resize oLabels
	oLabels.clear(); oLabels.resize(iTrajClusters.size());

	set<int> shiptypes;	
	for (const auto& traj : iTrajClusters) shiptypes.insert(traj.staticInfo.shiptypecode);
	
	map<int, size_t> shiptypes_map; 
	size_t cnt = 0;
	for (const auto& shiptype : shiptypes)
	{
		shiptypes_map.insert({ shiptype, cnt });
		cnt++;
	}

	for (size_t i = 0; i < iTrajClusters.size(); i++) 
		oLabels[i] = shiptypes_map[iTrajClusters[i].staticInfo.shiptypecode];		
}

void PreClustering_byShipType(
	const vector<AisTrajectory<>>& iTrajClusters,
	vector<size_t>& iPreLabels,
	vector<size_t>& oLabels)
{	
	assert(iTrajClusters.size() == iPreLabels.size());
	if (iTrajClusters.size() != iPreLabels.size()) throw;

	// resize oLabels
	oLabels.clear(); oLabels.resize(iTrajClusters.size());

	auto toStrCode = [](int shipcode, size_t prelabel) 
	{ 
		return to_string(shipcode) + "_" + to_string(prelabel); 
	};

	set<string> shiptypes;
	for (size_t i = 0; i < iTrajClusters.size(); i++){
		shiptypes.insert(
			toStrCode(iTrajClusters[i].staticInfo.shiptypecode, iPreLabels[i])
		);
	}

	map<string, size_t> shiptypes_map;
	size_t cnt = 0;
	for (const auto& shiptype : shiptypes)
	{
		shiptypes_map.insert({ shiptype, cnt });
		cnt++;
	}

	for (size_t i = 0; i < iTrajClusters.size(); i++) {
		oLabels[i] = shiptypes_map[
			toStrCode(iTrajClusters[i].staticInfo.shiptypecode, iPreLabels[i])
		];
	}
}

void PreClustering_byDeptArri(
	const vector<AisTrajectory<>>& iTrajClusters,
	vector<size_t>& oLabels,
	Position Dept, Position Arri,
	Real radius, size_t minPts,
	bool enableReverse = false)
{
	// assertions
	assert(iTrajClusters.size() == oLabels.size());
	if (iTrajClusters.size() != oLabels.size()) throw;

	//run DBSCAN to cluster with dept/arri
}

void PreClustering_byDeptArri(
	const vector<AisTrajectory<>>& iTrajClusters,
	vector<size_t>& iPreLabels,
	vector<size_t>& oLabels,
	Position Dept, Position Arri,
	Real radius, size_t minPts,
	bool enableReverse = false)
{
	// assertions
	assert(iTrajClusters.size() == iPreLabels.size());
	if (iTrajClusters.size() != iPreLabels.size()) throw;
	const size_t& N = iTrajClusters.size();

	size_t C = GetClusterNum(iPreLabels);
	vector<vector<size_t>> tmpLabels(C, vector<size_t>());

	// define distance 
	struct DeptArriDist : Clustering::DistanceMeasureBase<AisTrajectory<>>
	{

	};

	#pragma omp parallel for
	for (size_t c = 0; c < C; c++)
	{
		vector<AisTrajectory<>> tmpCluster;		
		for (size_t i = 0; i < N; i++) {
			if (iPreLabels[i] == c) tmpCluster.push_back(iTrajClusters[i]);
		}

		
	}

}

enum class TrajOutputFileType { FelineJson, GeoJson };

void ExportByClusters(
	const vector<AisTrajectory<>>& iTrajClusters,
	vector<size_t>& labels, 
	const string& OutputDirPath, TrajOutputFileType type)
{
	assert(iTrajClusters.size() == labels.size());
	if (iTrajClusters.size() != labels.size()) throw;
	const size_t& N = iTrajClusters.size();

	size_t C = GetClusterNum(labels);

	#pragma omp parallel for
	for (size_t c = 0; c < C; c++)	
	{
		vector<AisTrajectory<>> tmpCluster;
		for (size_t i = 0; i < N; i++) {
			if (labels[i] == c) tmpCluster.push_back(iTrajClusters[i]);
		}

		const string filepath = OutputDirPath + "/Cluster_" + to_string(c) + ".json";
		if (type == TrajOutputFileType::FelineJson)
			IO::WriteFelineJsonFile(filepath, tmpCluster);
		else
			IO::WriteGeoJsonFile(filepath, tmpCluster);
	}	
}




void run() 
{
	CommonLogger logger;
	logger.Log({ {Tag::lvl, 1} }) << "Maritime Traffic Map building started." << endl;

	// Read Trajectories
	vector<string> felineJsonFilenames = Fs::GetFilesFromPaths(_c.inputTrajFilepath, ".*\\.json");
	vector<AisTrajectory<>> inputTrajs;
	vector<vector<AisTrajectory<>>> rawInputTrajs(felineJsonFilenames.size(), vector<AisTrajectory<>>());
	#pragma omp parallel for
	for (size_t i = 0; i < felineJsonFilenames.size(); i++)
	{
		rawInputTrajs[i] = IO::ReadFelineJsonFile<AisTrajectory<>>(felineJsonFilenames[i]);
	}
	for (size_t i = 0; i < felineJsonFilenames.size(); i++)
		inputTrajs.insert(end(inputTrajs), begin(rawInputTrajs[i]), end(rawInputTrajs[i]));
	rawInputTrajs.clear();

	// trajectory clusters to extract representative routes
	vector<size_t> trajClusterLabels;

	Clustering::DistanceBasedDBSCAN test;
	//test


	
}	