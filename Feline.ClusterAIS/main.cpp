#include <HashColon/HashColon_config.h>
#include <any>
#include <cmath>
#include <functional>
#include <limits> 
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <Eigen/Eigen>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/writer.h>
#include <HashColon/Clustering.hpp>
#include <HashColon/Exception.hpp>
#include <HashColon/Helper.hpp>
#include <HashColon/Log.hpp>
#include <HashColon/Real.hpp>
#include <HashColon/SingletonCLI.hpp>
#include <HashColon/Feline/FelineJsonIO.hpp>
#include <HashColon/Feline/TrajectoryClustering.hpp>
#include <HashColon/Feline/XtdTrajectoryClustering.hpp>

using namespace std;
using namespace Eigen;
using namespace HashColon;
using namespace HashColon::Fs;
using namespace HashColon::Clustering;
using namespace HashColon::Feline;

using Tag = HashColon::LogUtils::Tag;
namespace Json = rapidjson;

//#define DEBUG

HASHCOLON_NAMED_EXCEPTION_DEFINITION(main);

struct _Params
{
	std::string experimentName = "";
	std::string experimentComments = "";
	std::string experimentDatetime = "";	// not input: set in preprocess

	bool enableNoiseCluster = false;	// not input: set in clustering methods
	bool enableEvaluation = true;
	bool enableBuildReportFiles = true;
	bool enableScreenResults = false;
	bool enableExportDistMatrix = true;
	bool enablePreComputedDistMat = false;

	std::string inputDir = "";
	std::string outputDir = "";
	std::string predefinedLabelPath = "";
	
	std::string trajectoryTypeName = "";
	std::string clusteringMethodName = "";
	std::string measuringMethodName = "";	
	Real reptModeRadius = 0.0;

	bool enableUniformSampling = false;
	std::size_t UniformSamplingNumber = 0;
};

struct TrajOriginList
{
	std::string filename = "";
	int routeNum = -1;
	int clusterTag = -1;
};

static _Params _c;
static Json::MemoryPoolAllocator<> alloc;

template <typename TrajT>
struct Method
{
	typename DistanceMeasureBase<TrajT>::Ptr measure(string name) const;
	typename DistanceBasedClustering<TrajT>::Ptr clustering(
		string name, typename DistanceMeasureBase<TrajT>::Ptr measureFunc) const;
};

// similarity/distance methods, clustering methods for XYList
template <>
struct Method<XYList>
{
	using TrajT = XYList;

	// similarity/distance methods
	typename DistanceMeasureBase<TrajT>::Ptr measure(string name) const
	{
		using namespace HashColon::Feline::TrajectoryClustering;

		if (name == "Euclidean") { return make_shared<Euclidean>(); }
		if (name == "Hausdorff") { return make_shared<Hausdorff>(); }
		if (name == "LCSS") { return make_shared<LCSS>(); }
		if (name == "Merge") { return make_shared<Merge>(); }
		if (name == "ProjectedPCA") { return make_shared<ProjectedPCA>(); }
		if (name == "ModifiedHausdorff") { return make_shared<ModifiedHausdorff>(); }
		if (name == "DynamicTimeWarping") { return make_shared<DynamicTimeWarping>(); }

		throw mainException(
			"Invalid distance/similarity measuring method name[" + name + "]. Check option --measuringMethodName.",
			__CODEINFO__);
	}

	// clustering methods
	typename DistanceBasedClustering<TrajT>::Ptr clustering(
		string name, typename DistanceMeasureBase<TrajT>::Ptr measureFunc) const
	{
		using namespace HashColon::Feline::TrajectoryClustering;

		if (name == "NJW") { _c.enableNoiseCluster = false; return make_shared<NJW<TrajT>>(measureFunc); }
		if (name == "DistanceBasedDBSCAN") { _c.enableNoiseCluster = true; return make_shared<DistanceBasedDBSCAN<TrajT>>(measureFunc); }

		throw mainException(
			"Invalid clustering method name[" + name + "]. Check option --clusteringMethodName.",
			__CODEINFO__);
	}
};

// similarity/distance methods, clustering methods for XYXtdList
template <>
struct Method<XYXtdList>
{
	using TrajT = XYXtdList;

	// similarity/distance methods
	typename DistanceMeasureBase<TrajT>::Ptr measure(string name) const
	{
		using namespace HashColon::Feline::XtdTrajectoryClustering;

		if (name == "DtwXtd") { return make_shared<DtwXtd>(); }
		if (name == "DtwXtd_JS") { return make_shared<DtwXtd_usingJSDivergence>(); }
		if (name == "DtwXtd_EMD") { return make_shared<DtwXtd_usingWasserstein>(); }
		if (name == "DtwXtd_Blended") { return make_shared<DtwXtd_BlendedDistance>(); }

		throw mainException(
			"Invalid distance/similarity measuring method name[" + name + "]. Check option --measuringMethodName.",
			__CODEINFO__);
	}

	// clustering methods
	DistanceBasedClustering<TrajT>::Ptr clustering(
		string name, typename DistanceMeasureBase<TrajT>::Ptr measureFunc) const
	{
		using namespace HashColon::Feline::XtdTrajectoryClustering;

		if (name == "NJW") {
			_c.enableNoiseCluster = false; 
			return make_shared<NJW<TrajT>>(measureFunc); 
		}
		else if (name == "DistanceBasedDBSCAN") {
			_c.enableNoiseCluster = true; 
			return make_shared<DistanceBasedDBSCAN<TrajT>>(measureFunc); 
		}
		else 
			throw mainException(
				"Invalid clustering method name[" + name + "]. Check option --clusteringMethodName.",
				__CODEINFO__);
	}
};

vector<XYList> run_uniformSampling(const vector<XYList>& routes)
{
	vector<XYList> re(routes.size());
	#pragma omp parallel for
	for (size_t i = 0; i < routes.size(); i++)
	{
		re[i] = routes[i].GetUniformLengthSampled(_c.UniformSamplingNumber);		
	}
	return re;
}

vector<XYXtdList> run_uniformSampling(const vector<XYXtdList>& routes)
{
	vector<XYXtdList> re(routes.size());
	#pragma omp parallel for
	for (size_t i = 0; i < routes.size(); i++)
	{
		re[i] = routes[i].GetUniformLengthSampled(_c.UniformSamplingNumber);		
	}
	return re;
}

void run_additionals(vector<XYList>& routes,
	typename DistanceBasedClustering<XYList>::Ptr clusteringMethod)
{
	// if measurement method is projectd PCA, run PCA
	if (_c.measuringMethodName == "ProjectedPCA")
	{
		using namespace HashColon::Feline::TrajectoryClustering;

		shared_ptr<DistanceBasedClustering<XYList>> a;


		dynamic_pointer_cast<ProjectedPCA>(
			clusteringMethod->GetDistanceFunc()
			)->RunPCA(routes);
	}
}

void run_additionals(vector<XYXtdList>& routes,
	typename DistanceBasedClustering<XYXtdList>::Ptr clusteringMethod)
{

}

namespace BuildResult
{
	void Inputs_ParamsJson(const string& InputsDir, const string& experimentName)
	{
		string Path_Params = InputsDir + "/" + experimentName + "_Params.json";

		Json::Document doc;
		doc.Parse(
			SingletonCLI::GetInstance().GetCLI()->config_to_str().c_str()
		);
		ofstream ofs(Path_Params);
		Json::OStreamWrapper osw(ofs);
		Json::PrettyWriter<Json::OStreamWrapper> writer(osw);
		doc.Accept(writer);
	}

	void Inputs_FelineJson(const string& InputsDir, const string& inDir)
	{
		string Dir_Inputs_FelineJson = InputsDir + "/FelineJson";
		BuildDirectoryStructure(Dir_Inputs_FelineJson);		
		filesystem::copy(inDir, Dir_Inputs_FelineJson, 
			filesystem::copy_options::recursive | 
			filesystem::copy_options::update_existing);
	}

	void Inputs_GeoJson(const string& InputsDir)
	{
		string Dir_Inputs_GeoJson = InputsDir + "/GeoJson";
		BuildDirectoryStructure(Dir_Inputs_GeoJson);		
	}

	Json::Document Report(
		const string& ReportDir, const _Params& c,
		const vector<TrajOriginList>& TrajList,
		const Eigen::MatrixXR& DistMatrix, 
		const vector<vector<size_t>>& ClusterResults,				
		const vector<Real>& DaviesBouldinCluster,
		const vector<vector<Real>>& Silhouettes)
	{
		Json::Document doc;		
		doc.Parse("{}");

		// Add Experiment informations
		Json::Value ValExpName; ValExpName.SetString(c.experimentName.c_str(), alloc);
		doc.AddMember("ExperimentName", ValExpName, alloc);
		Json::Value ValExpComments; ValExpComments.SetString(c.experimentComments.c_str(), alloc);
		doc.AddMember("ExperimentComments", ValExpComments, alloc);
		Json::Value ValExpDatetime; ValExpDatetime.SetString(c.experimentDatetime.c_str(), alloc);
		doc.AddMember("ExperimentDatetime", ValExpDatetime, alloc);

		// Add Input parameters
		Json::Document InputParams_doc;
		InputParams_doc.Parse(
			SingletonCLI::GetInstance().GetCLI()->config_to_str().c_str()
		);		
		// remove duplicate/unused/needless items
		if (InputParams_doc.HasMember("Log")) InputParams_doc.RemoveMember("Log");
		if (InputParams_doc.HasMember("experimentName")) InputParams_doc.RemoveMember("experimentName");
		if (InputParams_doc.HasMember("experimentComments")) InputParams_doc.RemoveMember("experimentComments");		
		if (InputParams_doc.HasMember("enableEvaluation")) InputParams_doc.RemoveMember("enableEvaluation");
		if (InputParams_doc.HasMember("enableBuildReportFiles")) InputParams_doc.RemoveMember("enableBuildReportFiles");
		if (InputParams_doc.HasMember("enableScreenResults")) InputParams_doc.RemoveMember("enableScreenResults");
		if (InputParams_doc.HasMember("Clustering")) 
		{
			vector<string> delList;
			for (auto& item : InputParams_doc["Clustering"].GetObject())
				delList.push_back(string(item.name.GetString()));
			string methodname = c.clusteringMethodName + "_" + c.trajectoryTypeName;
			delList.erase(remove(delList.begin(), delList.end(), methodname), delList.end());
			for(auto& itemName : delList)
				InputParams_doc["Clustering"].RemoveMember(itemName.c_str());
		}
		// remove unused parameters from feline parameters
		if (InputParams_doc.HasMember("Feline"))
		{			
			string measureType;
			if (string(InputParams_doc["trajectoryTypeName"].GetString()) == "XYList")
			{
				// if traj type is "XYList", erase "XtdTrajectoryDistanceMeasure" objs
				InputParams_doc["Feline"].RemoveMember("XtdTrajectoryDistanceMeasure");
				measureType = "TrajectoryDistanceMeasure";

			
			}
			else if (string(InputParams_doc["trajectoryTypeName"].GetString()) == "XYXtdList")
			{
				// if traj type is "XYList", erase "XtdTrajectoryDistanceMeasure" objs
				InputParams_doc["Feline"].RemoveMember("TrajectoryDistanceMeasure");
				measureType = "XtdTrajectoryDistanceMeasure";				
			}
			// delete methods except used method
			vector<string> delList;
			for (auto& item : InputParams_doc["Feline"][measureType.c_str()].GetObject())
				delList.push_back(string(item.name.GetString()));
			delList.erase(remove(delList.begin(), delList.end(), "Enable_ReversedSequence"), delList.end());
			delList.erase(remove(delList.begin(), delList.end(), c.measuringMethodName), delList.end());
			for (auto& itemName : delList)
				InputParams_doc["Feline"][measureType.c_str()].RemoveMember(itemName.c_str());
		}
		
		doc.AddMember("InputParams", InputParams_doc, alloc);

		Json::Value ValIndexList(Json::kArrayType);
		for (const auto& Traj : TrajList)
		{
			Json::Value item(Json::kArrayType);
			Json::Value fname; fname.SetString(Traj.filename.c_str(), alloc);
			item.PushBack(fname, alloc);
			item.PushBack(Json::Value(Traj.routeNum), alloc);
			ValIndexList.PushBack(item, alloc);
		}
		doc.AddMember("IndexList", ValIndexList, alloc);

		Json::Value ValDistMatrix(Json::kArrayType);
		for (int i = 0; i < DistMatrix.rows(); i++)
		{
			Json::Value ValRow(Json::kArrayType);
			for (int j = 0; j < DistMatrix.cols(); j++)
			{
				ValRow.PushBack(Json::Value(DistMatrix(i, j)), alloc);
			}
			ValDistMatrix.PushBack(ValRow, alloc);
		}
		doc.AddMember("DistanceMatrix", ValDistMatrix, alloc);

		doc.AddMember("IsNoiseClusterDefined", Json::Value(c.enableNoiseCluster), alloc);

		Json::Value ValClusterNames(Json::kArrayType);
		for (size_t i = 0; i < ClusterResults.size(); i++)
		{
			if (c.enableNoiseCluster)
				if (i == 0)
					ValClusterNames.PushBack(Json::Value("Noise"), alloc);
				else
				{
					string clusterName = "Cluster_" + to_string(i - 1);
					Json::Value ValClusterName; ValClusterName.SetString(clusterName.c_str(), alloc);
					ValClusterNames.PushBack(ValClusterName, alloc);
				}
			else
			{
				string clusterName = "Cluster_" + to_string(i);
				Json::Value ValClusterName; ValClusterName.SetString(clusterName.c_str(), alloc);
				ValClusterNames.PushBack(ValClusterName, alloc);
			}
		}
		doc.AddMember("ClusterNames", ValClusterNames, alloc);

		Json::Value ValClusteringResults(Json::kArrayType);
		for (const auto& cluster : ClusterResults)
		{
			Json::Value ValCluster(Json::kArrayType);
			for (const auto& item : cluster)
			{
				ValCluster.PushBack(Json::Value(item), alloc);
			}
			ValClusteringResults.PushBack(ValCluster, alloc);
		}
		doc.AddMember("ClusteringResults", ValClusteringResults, alloc);	

		doc.AddMember("Evaluation", Json::Value(Json::kObjectType), alloc);
		Json::Value ValDaviesBouldin(Json::kArrayType);
		for (const auto& db : DaviesBouldinCluster)
		{
			ValDaviesBouldin.PushBack(Json::Value(db), alloc);
		}
		doc["Evaluation"].AddMember("DaviesBouldin", ValDaviesBouldin, alloc);
		
		Json::Value ValSilhouettes(Json::kArrayType);
		for (const auto& SilC : Silhouettes)
		{
			Json::Value Sil_byCluster(Json::kArrayType);
			for (const auto& item : SilC)
			{
				Sil_byCluster.PushBack(Json::Value(item), alloc);
			}
			ValSilhouettes.PushBack(Sil_byCluster, alloc);
		}
		doc["Evaluation"].AddMember("Silhouettes", ValSilhouettes, alloc);

		string Path_ReportJson = ReportDir + "/" + c.experimentName + "_Report.json";
		ofstream ofs(Path_ReportJson);
		Json::OStreamWrapper osw(ofs);
		Json::PrettyWriter<Json::OStreamWrapper> writer(osw);
		bool jsonWriterSuccess = doc.Accept(writer);
		if (!jsonWriterSuccess)
		{
			throw mainException("Error while writing report json file. Invalid value nan/inf found in results.");
		}

		return doc;
	}

	template <typename TList>
	Json::Document Results(
		const string& ResultGeoJsonDir, 
		const vector<vector<TList>>& Clusters,
		const vector<TList>& Representatives,
		const string& experimentName, const bool hasNoiseCluster)
	{
		Json::Value re(Json::kArrayType);
		Json::Document doc;		
		doc.SetObject();		

		for (size_t i = 0; i < Clusters.size(); i++)
		{	
			if (hasNoiseCluster)
			{
				if (i == 0)
				{
					string filename = experimentName + "_Noise.json";
					string filepath = ResultGeoJsonDir + "/" + filename;
					Json::Document re = IO::WriteGeoJsonFile(filepath, Clusters[i], alloc);
					Json::Value ValFilename; ValFilename.SetString(filename.c_str(), alloc);
					doc.AddMember(ValFilename, re.GetObject(), alloc);
				}
				else
				{
					string filename1 = experimentName + "_Cluster_" + to_string(i - 1) + ".json";
					string filepath1 = ResultGeoJsonDir + "/" + filename1;
					Json::Document re1 = IO::WriteGeoJsonFile(filepath1, Clusters[i], alloc);
					Json::Value ValFilename1; ValFilename1.SetString(filename1.c_str(), alloc);
					doc.AddMember(ValFilename1, re1.GetObject(), alloc);

					string filename2 = experimentName + "_Representative_Cluster_" + to_string(i - 1) + ".json";
					string filepath2 = ResultGeoJsonDir + "/" + filename2;
					vector<TList> representative; representative.push_back(Representatives[i]);
					Json::Document re2 = IO::WriteGeoJsonFile(filepath2, representative, alloc);
					Json::Value ValFilename2; ValFilename2.SetString(filename2.c_str(), alloc);
					doc.AddMember(ValFilename2, re2.GetObject(), alloc);
				}
			}
			else
			{
				string filename1 = experimentName + "_Cluster_" + to_string(i) + ".json";
				string filepath1 = ResultGeoJsonDir + "/" + filename1;
				Json::Document re1 = IO::WriteGeoJsonFile(filepath1, Clusters[i], alloc);
				Json::Value ValFilename1; ValFilename1.SetString(filename1.c_str(), alloc);
				doc.AddMember(ValFilename1, re1.GetObject(), alloc);

				string filename2 = experimentName + "_Representative_Cluster_" + to_string(i - 1) + ".json";
				string filepath2 = ResultGeoJsonDir + "/" + filename2;
				vector<TList> representative; representative.push_back(Representatives[i]);
				Json::Document re2 = IO::WriteGeoJsonFile(filepath2, representative, alloc);
				Json::Value ValFilename2; ValFilename2.SetString(filename2.c_str(), alloc);
				doc.AddMember(ValFilename2, re2.GetObject(), alloc);
			}
		}
		return doc;
	}
		
	Json::Document ExportDistMatrix(const string& ResultDir, const Eigen::MatrixXR& DistMatrix)
	{
		Json::Document doc;
		doc.Parse("{}");

		Json::Value ValDistMatrix(Json::kArrayType);
		for (int i = 0; i < DistMatrix.rows(); i++)
		{
			Json::Value ValRow(Json::kArrayType);
			for (int j = 0; j < DistMatrix.cols(); j++)
			{
				ValRow.PushBack(Json::Value(DistMatrix(i, j)), alloc);
			}
			ValDistMatrix.PushBack(ValRow, alloc);
		}
		doc.AddMember("DistanceMatrix", ValDistMatrix, alloc);

		string Path_DistMatrixJson = ResultDir + "/DistMatrix.json";
		ofstream ofs(Path_DistMatrixJson);
		Json::OStreamWrapper osw(ofs);
		Json::PrettyWriter<Json::OStreamWrapper> writer(osw);
		bool jsonWriterSuccess = doc.Accept(writer);
		if (!jsonWriterSuccess)
		{
			throw mainException("Error while writing report json file. Invalid value nan/inf found in results.");
		}

		return doc;
	}

	void Preprocess(const _Params& c)
	{
		// stamp experiment datetime
		auto now = chrono::system_clock::now();
		auto now_time_t = chrono::system_clock::to_time_t(now);
		stringstream ss; 
		ss << put_time(localtime(&now_time_t), "%Y-%m-%d %H:%M:%S");
		_c.experimentDatetime = ss.str();

		// build root dir
		string Dir_root = _c.outputDir + "/" + _c.experimentName;
		if (!BuildDirectoryStructure(Dir_root))
		{
			throw mainException(
				"Output directory does not exist or experiment directory cannot be built["
				+ Dir_root + "]. Check option --outputDir, --experimentName.",
				__CODEINFO__);
		}

		// clear root directory 
		//RemoveAllInDirectory(Dir_root);
		
		// build inputs
		string Dir_Inputs = Dir_root + "/Inputs";
		BuildDirectoryStructure(Dir_Inputs);
		RemoveAllInDirectory(Dir_Inputs);
		Inputs_FelineJson(Dir_Inputs, _c.inputDir);
		Inputs_GeoJson(Dir_Inputs);
		Inputs_ParamsJson(Dir_Inputs, _c.experimentName);

		// intercept Log files from commonlogger
		string Dir_Log = Dir_root + "/Log";
		BuildDirectoryStructure(Dir_Log);
		//RemoveAllInDirectory(Dir_Log);

		string Path_logFile = Dir_Log + "/Log.log";
		shared_ptr<ofstream> logFile = make_shared<ofstream>(Path_logFile);
		CommonLogger::GetDefaultParams().LogFile = logFile;

		string Path_errlogFile = Dir_Log + "/ErrorLog.log";
		shared_ptr<ofstream> errlogFile = make_shared<ofstream>(Path_errlogFile);
		CommonLogger::GetDefaultParams().ErrorFile = errlogFile;
	}

	template <typename TList>
	void PrintResult(const _Params& c,
		const vector<TList>& routes, const size_t& numOfClusters,
		const vector<TrajOriginList>& TrajList,
		const Eigen::MatrixXR& DistMatrix, 
		const vector<size_t>& ClusterResults,
		const vector<TList>& RepresentativeTrajs,
		const shared_ptr<vector<Real>>& DaviesBouldin,
		const shared_ptr<vector<vector<Real>>>& Silhouettes)
	{
		string Dir_root = _c.outputDir + "/" + _c.experimentName;
		Json::Document doc;		
		doc.SetObject();

		// Report 
		if (_c.enableBuildReportFiles)
		{
			// build Report directory
			string Dir_Report = Dir_root + "/Report";
			BuildDirectoryStructure(Dir_Report);
			
			// Gather clustering results by clusters
			vector<vector<size_t>> clusterResults_byCluster;
			clusterResults_byCluster.resize(numOfClusters);
			for (size_t i = 0; i < ClusterResults.size(); i++)
				clusterResults_byCluster[ClusterResults[i]].push_back(i);				
			
			// Print out report
			Json::Document DocReport = Report(
				Dir_Report, c, TrajList, DistMatrix, 
				clusterResults_byCluster,
				(*DaviesBouldin), (*Silhouettes) );						
			if (c.enableScreenResults)
			{
				string name = c.experimentName + "_Report.json";
				Json::Value ValName; ValName.SetString(name.c_str(), alloc);				
				doc.AddMember(ValName, DocReport, alloc);
			}
		}

		// Results		
		const string Dir_Results = Dir_root + "/Results";
		const string Dir_Results_GeoJson = Dir_Results + "/GeoJson";
		BuildDirectoryStructure(Dir_Results);
		//RemoveAllInDirectory(Dir_Results_GeoJson);
		BuildDirectoryStructure(Dir_Results_GeoJson);
		

		// DistMatrix output
		if (_c.enableExportDistMatrix)
			Json::Document DocDistMatrix = ExportDistMatrix(Dir_Results, DistMatrix);		

		// Gather clustering results by clusters
		vector<vector<TList>> clusteredTrajs;
		clusteredTrajs.resize(numOfClusters);
		for (size_t i = 0; i < ClusterResults.size(); i++)
			clusteredTrajs[ClusterResults[i]].push_back(routes[i]);			
		Json::Document DocResults = Results(
			Dir_Results_GeoJson, clusteredTrajs, RepresentativeTrajs,
			c.experimentName, c.enableNoiseCluster);		
		if (c.enableScreenResults)
		{
			doc.AddMember("Clusters", DocResults, alloc);
		}

		// if screen output enabled, print
		if (c.enableScreenResults)
		{			
			cout << endl;
			Json::OStreamWrapper osw(cout);
			Json::PrettyWriter<Json::OStreamWrapper> writer(osw);
			doc.Accept(writer);			
			cout << endl;
		}
	}
}

bool ReadDistMatrix(MatrixXR& DistMatrix)
{	
	const string filepath = _c.outputDir + "/" + _c.experimentName + "/Results/DistMatrix.json";
	
	// try reading file
	ifstream ifs(filepath);
	if (ifs.fail()) return false;

	Json::IStreamWrapper isw(ifs);
	Json::Document doc;
	Json::ParseResult parseSucceed = doc.ParseStream(isw);
	if (parseSucceed)
	{
		// check validity
		if (!doc.IsObject()) return false;
		if (!doc.HasMember("DistanceMatrix")) return false;
		Json::Value& ValDistMatrix = doc["DistanceMatrix"];
		if (!ValDistMatrix.IsArray()) return false;

		// resize distmatrix
		size_t N = ValDistMatrix.GetArray().Size();
		if (N > 0)
			DistMatrix.resize(N, N);
		else
			return false;
		
		// loop and read
		size_t i = 0, j = 0;
		for (auto& row : ValDistMatrix.GetArray())
		{			
			if (!row.IsArray()) return false;
			j = 0;
			for (auto& item : row.GetArray())
			{
				if (!item.IsDouble()) return false;
				DistMatrix(i, j) = item.GetDouble();
				j++;
			}
			i++;
		}

		return true;
	}
	else {		
		return false;
	}
		
}

template <typename TList>
void ExportPseudoMode(
	const std::vector<TList>& inputTraj,
	const std::vector<size_t>& clusterResult, 
	const Eigen::MatrixXR& DistanceMatrix, 
	const Real radius)
{
	CommonLogger logger;	

	size_t N = (size_t)DistanceMatrix.rows();
	size_t C = (*max_element(clusterResult.begin(), clusterResult.end())) + 1;
	vector<size_t> modeCnt(C, 0);
	vector<vector<size_t>> modeIndices(C, vector<size_t>());	
	
	for (size_t i = 0; i < N; i++)
	{
		size_t cnt = 0;
		size_t clusterIdx = clusterResult[i];
		if (_c.enableNoiseCluster && clusterIdx == 0) continue;
		for (size_t j = 0; j < N; j++)
		{
			if (i == j) continue;		
			double test = DistanceMatrix(i, j);
			if (DistanceMatrix(i, j) <= radius) 
			{
				cnt++;
			}
		}	
		
		if (cnt > modeCnt[clusterIdx])
		{
			modeCnt[clusterIdx] = cnt;
			modeIndices[clusterIdx].clear();
			modeIndices[clusterIdx].push_back(i);			
		}
		else if (cnt == modeCnt[clusterIdx])
		{
			modeIndices[clusterIdx].push_back(i);
		}
	}

	string Dir_root = _c.outputDir + "/" + _c.experimentName;
	string Dir_Results = Dir_root + "/Results";
	string Dir_Results_GeoJson = Dir_Results + "/GeoJson";
	BuildDirectoryStructure(Dir_Results);
	RemoveAllInDirectory(Dir_Results_GeoJson);
	BuildDirectoryStructure(Dir_Results_GeoJson);

	const size_t noiseIdxFactor = _c.enableNoiseCluster ? 1 : 0;
	const string dirpath = _c.outputDir + "/" + _c.experimentName + "/Results/GeoJson";		

	for (size_t i = noiseIdxFactor; i < C; i++)
	{
		vector<TList> outputTrajs;
		outputTrajs.clear();
		for (const auto& idx : modeIndices[i]) {
			outputTrajs.push_back(inputTraj[idx]);
		}		
		logger.Log() << "PseudoMode: " << "cluster " << (i - noiseIdxFactor)
			<< "max count : " << modeCnt[i] << " / number of modes : " << modeIndices[i].size() << endl;

		string filename = _c.experimentName + "_Mode_Cluster_" + to_string(i - noiseIdxFactor) + ".json";
		string filepath = dirpath + "/" + filename;		
		IO::WriteGeoJsonFile(filepath, outputTrajs);		
	}
}

template <typename TList>
void run_core(
	const vector<string>& inFilePaths, string& outDir)
{
	CommonLogger logger;
	vector<vector<TList>> routes_raw; routes_raw.resize(inFilePaths.size());

	// parse json to routes	
	logger.Log({ {Tag::lvl, 1} }) << "Starting the process." << endl;
#pragma omp parallel for 
	for (size_t i = 0; i < inFilePaths.size(); i++)
	{
		routes_raw[i] = IO::ReadFelineJsonFile<TList>(inFilePaths[i]);
		// if report option is on, export raw routes to geojson
		if (_c.enableBuildReportFiles)
		{
			string filename = String::Split(inFilePaths[i], "/").back();
			string filepath = _c.outputDir + "/" + _c.experimentName + "/Inputs/GeoJson/" + filename;
			IO::WriteGeoJsonFile(filepath, routes_raw[i]);
		}
	}

	// merge routes into a array
	vector<TList> routes_input;
	vector<TrajOriginList> trajOrigins;
	for (size_t i = 0; i < routes_raw.size(); i++)
	{
		for (size_t j = 0; j < routes_raw[i].size(); j++)
		{
			routes_input.push_back(routes_raw[i][j]);
			trajOrigins.push_back({ inFilePaths[i], (int)j });
		}
	}
	routes_raw.clear();
	logger.Log({ {Tag::lvl, 1} }) << "Reading input trajectories finished. " << routes_input.size() << " trajectories found." << endl;

	// check if input routes is not empty
	if (routes_input.size() <= 0) {
		logger.Log({ {Tag::lvl, 1} }) << "No trajectory is given. Clustering Finished." << endl;
		logger.Log({ {Tag::lvl, 1} }) << "Process finished." << endl;
		return;
	}

	// check if uniform sampling enabled.
	vector<TList> routes;
	if (_c.enableUniformSampling)
	{
		routes = run_uniformSampling(routes_input);
		logger.Log({ {Tag::lvl, 1} }) << "Uniform sampling done. (Sample #: " << _c.UniformSamplingNumber << ")" << endl;
	}
	else
		routes = routes_input;

	// set clustering method.
	Method<TList> method;
	typename DistanceBasedClustering<TList>::Ptr clusteringMethod =
		method.clustering(_c.clusteringMethodName,
			method.measure(_c.measuringMethodName));

	// check if additional process is needed
	run_additionals(routes, clusteringMethod);

	// run clustering
	shared_ptr<vector<size_t>> labels = make_shared<vector<size_t>>();
	MatrixXR DistMatrix;
	if (_c.enablePreComputedDistMat)
	{
		bool readDistMatSucceeded = ReadDistMatrix(DistMatrix);
		// if reading failed, build dist matrix
		if (!readDistMatSucceeded) 
		{
			logger.Log({ {Tag::lvl, 1} }) << "Pre-computed distance matrix failed. Start recomputeing distance matrix." << endl;
			DistMatrix = clusteringMethod->ComputeDistanceMatrix(routes, true);			
		}
		else
		{
			logger.Log({ {Tag::lvl, 1} }) << "Pre-computed distance matrix loaded." << endl;
		}		
	}
	else
	{
		logger.Log({ {Tag::lvl, 1} }) << "Start computeing distance matrix." << endl;
		DistMatrix = clusteringMethod->ComputeDistanceMatrix(routes, true);
	}

	clusteringMethod->TrainModel(
		DistMatrix,
		clusteringMethod->GetDistanceFunc()->GetMeasureType() == HashColon::Clustering::distance,
		labels);
	size_t numOfClusters = clusteringMethod->GetNumOfClusters();
	//clusteringMethod->TrainModel(routes, labels);

	logger.Log({ {Tag::lvl, 1} }) << "Clustering Trajectories finished. Total " << numOfClusters << " clusters." << endl;

	// check if input routes is not empty
	if (numOfClusters <= 0) {
		logger.Log({ {Tag::lvl, 1} }) << "No clusters are generated or all trajectories are assigned as noise." << endl;
		logger.Log({ {Tag::lvl, 1} }) << "Process finished." << endl;
		return;
	}

	// compute pseudo-median trajectory of clusters
	vector<size_t> reptTrajIdx = PseudoMedian((*labels), DistMatrix);
	vector<TList> reptTrajList; reptTrajList.resize(reptTrajIdx.size());
	for (size_t i = 0; i < reptTrajIdx.size(); i++)
	{
		reptTrajList[i] = routes_input[reptTrajIdx[i]];
	}
	logger.Log({ {Tag::lvl, 1} }) << "Representative trajectory extracted." << endl;

	// export mode trajs
	ExportPseudoMode(routes_input, (*labels), DistMatrix, _c.reptModeRadius);


	// compute evaluations	
	shared_ptr<vector<Real>> DaviesBouldin;
	shared_ptr<vector<vector<Real>>> SilhouettesByClusters;
	if (_c.enableEvaluation)
	{
		// compute Davies-Bouldin
		DaviesBouldin = make_shared<vector<Real>>(PseudoDaviesBouldin((*labels), DistMatrix, reptTrajIdx));
		// compute Silhouettes
		SilhouettesByClusters = make_shared<vector<vector<Real>>>(numOfClusters);
		vector<Real> Silhouettes = Silhouette((*labels), DistMatrix);
		for (size_t i = 0; i < Silhouettes.size(); i++)
		{
			SilhouettesByClusters->at(labels->at(i)).push_back(Silhouettes[i]);
		}
		logger.Log({ {Tag::lvl, 1} }) << "Computing cluster evaluation values finished." << endl;
	}

	// Print out results
	BuildResult::PrintResult(_c, routes_input, numOfClusters,
		trajOrigins, DistMatrix, (*labels), reptTrajList,
		DaviesBouldin, SilhouettesByClusters);
	if (_c.enableBuildReportFiles)
		logger.Log({ {Tag::lvl, 1} }) << "Printing outputs & building reports finished." << endl;
	else
		logger.Log({ {Tag::lvl, 1} }) << "Printing outputs finished." << endl;

	logger.Log({ {Tag::lvl, 1} }) << "Process finished." << endl;
}

void run()
{
	BuildResult::Preprocess(_c);

	CommonLogger logger;
	logger.Log({ {Tag::lvl, 1} }) << "Start!" << endl;

	//string configstr = SingletonCLI::GetInstance().GetCLI()->config_to_str();
	logger.Log({ {Tag::lvl, 1} }) << "Options: "
		<< SingletonCLI::GetInstance().GetCLI()->config_to_str() << endl;

	string dir = _c.inputDir;
	string outdir = _c.outputDir + "/" + _c.clusteringMethodName + "/" + _c.measuringMethodName;
	vector<string> filepaths = GetFilesInDirectory(dir);

	if (_c.trajectoryTypeName == "XYXtdList")
		run_core<XYXtdList>(filepaths, outdir);
	else if (_c.trajectoryTypeName == "XYList")
		run_core<XYList>(filepaths, outdir);
	else
		mainException("Invalid trajectory type name: " + _c.trajectoryTypeName + ". Check parameter --trajectoryTypeName.");
}

static void Initialize()
{
	// register config files	
	SingletonCLI::Initialize();
	CommonLogger::Initialize("./CommonLogger.json");

	GeoDistance::Initialize();

	NJW<XYList>::Initialize();
	DistanceBasedDBSCAN<XYList>::Initialize();
	TrajectoryClustering::Initialize_All_TrajectoryDistanceMeasure();

	NJW<XYXtdList>::Initialize();
	DistanceBasedDBSCAN<XYXtdList>::Initialize();
	XtdTrajectoryClustering::Initialize_All_XtdTrajectoryDistanceMeasure();

	CLI::App* cli = SingletonCLI::GetInstance().GetCLI();

	cli->add_option("--experimentName", _c.experimentName);
	cli->add_option("--experimentComments", _c.experimentComments);

	cli->add_option("--enableEvaluation", _c.enableEvaluation);
	cli->add_option("--enableBuildReportFiles", _c.enableBuildReportFiles);
	cli->add_option("--enableScreenResults", _c.enableScreenResults);
	cli->add_option("--enableExportDistMatrix", _c.enableExportDistMatrix);
	cli->add_option("--enablePreComputedDistMat", _c.enablePreComputedDistMat);	

	cli->add_option("--inputDir", _c.inputDir);
	cli->add_option("--outputDir", _c.outputDir);
	cli->add_option("--predefinedLabelPath", _c.predefinedLabelPath);

	cli->add_option("--clusteringMethodName", _c.clusteringMethodName);
	cli->add_option("--measuringMethodName", _c.measuringMethodName);
	cli->add_option("--trajectoryTypeName", _c.trajectoryTypeName);	
	cli->add_option("--reptModeRadius", _c.reptModeRadius);	

	cli->add_option("--enableUniformSampling", _c.enableUniformSampling);
	cli->add_option("--UniformSamplingNumber", _c.UniformSamplingNumber);	

	// Parse additional options from command line.
	SingletonCLI::GetInstance().GetCLI()->set_help_all_flag("--help-all");

	cli->callback(run);
}

int main(int argc, char** argv)
{
	Initialize();

	try
	{
		SingletonCLI::GetInstance().Parse(
			argc, argv,
			{
				//"./config/CommonLogger.json",
				"./TrajectoryClustering.json"
			});
	}
	catch (const CLI::Error& e)
	{
		SingletonCLI::GetInstance().GetCLI()->exit(e);
	}
	catch (const HashColon::Exception& e)
	{
		CommonLogger logger;
		logger.Error({ { Tag::file, e.file()}, { Tag::func, e.func() }, { Tag::line, e.line() } }) << e.what() << std::endl;		
	}

}