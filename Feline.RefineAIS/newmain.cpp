#include <HashColon/HashColon_config.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <HashColon/Exception.hpp>
#include <HashColon/Helper.hpp>
#include <HashColon/Log.hpp>
#include <HashColon/SingletonCLI.hpp>
#include <HashColon/Table.hpp>
#include <HashColon/ThreadPool.hpp>
#include <HashColon/Feline/AisPreprocess.hpp>
#include <HashColon/Feline/FelineJsonIO.hpp>
#include <HashColon/Feline/GeoData.hpp>
#include <HashColon/Feline/RouteSimplification.hpp>
#include <HashColon/Feline/XtdEstimation.hpp>

using namespace std;
namespace Json = rapidjson;
using namespace HashColon;
using namespace HashColon::Feline;
using namespace HashColon::LogUtils;

HASHCOLON_NAMED_EXCEPTION_DEFINITION(main);

struct _Params
{	
	int numOfThreads;

	bool runCsvRefinement;
	struct {		
		bool saveToFiles;
		string importDirectory_dynamic;		
		bool hasStaticFile;
		string importDirectory_static;
		string exportDirectory;
	} CsvRefinement;

	bool runTrajectoryExtraction;
	struct {		
		AisDataType dataType;
		bool saveToFiles;
		string importDirectory;
		string exportDirectory;
	} TrajectoryExtraction;

	bool runTrajectorySimplification;
	struct {		
		Real epsilon;
		bool saveToFiles;
		string importDirectory;
		string exportDirectory;
	} TrajectorySimplification;

	bool runXtdGeneration;
	struct {
		bool saveToFiles;
		string importDirectory;
		string exportDirectory;
	} XtdGeneration;
} _c;

vector<Table> run_CsvRefinement()
{
	// set params
	AisCsvReader::_Params params = AisCsvReader::GetDefaultParams();
	{
		params.threadCnt = _c.numOfThreads;
		params.aisPathAndDirectory = { _c.CsvRefinement.importDirectory_dynamic };
		if (_c.CsvRefinement.hasStaticFile)
		{
			params.staticPathAndDirectory = { _c.CsvRefinement.importDirectory_static };
		}
	}
	
	// run csv refinement
	AisCsvReader reader(params);	
	auto refinedCSV = reader.ReadAndRefineAisFiles_withLabel();
	// save to files
	if (_c.CsvRefinement.saveToFiles)
	{				
		#pragma omp parallel for	
		for(size_t i = 0 ; i < refinedCSV.size(); i++)		
		{
			auto it = refinedCSV.begin();
			advance(it, i);
			const auto& key = it->first;
			const auto& csv = it->second;
			string savefilepath = _c.CsvRefinement.exportDirectory + "/" + key + ".csv";
			csv.WriteToFile(savefilepath);
		}
	}	
	return reader.RemoveLabels(refinedCSV);
}

vector<Table> import_RefinedCsv(const string dir)
{
	// set params
	AisCsvReader::_Params params = AisCsvReader::GetDefaultParams();
	{
		params.threadCnt = _c.numOfThreads;
		params.aisPathAndDirectory = { dir };		
	}
	AisCsvReader reader(params);
	return reader.ReadAisFiles(nullptr);
}

vector<AisTrajectory<>> run_TrajectoryExtraction(vector<Table>& refinedCSV)
{
	// trajectory extraction
	shared_ptr<AisTrajectoryExtraction> extractor;
	if (_c.TrajectoryExtraction.dataType == ExactEarth)
	{
		extractor = make_shared<AisTrajectoryExtraction_ExactEarth>();
	}
	else if (_c.TrajectoryExtraction.dataType == KMOF)
	{
		extractor = make_shared<AisTrajectoryExtraction_KMOF>();
	}
	else
	{
		throw mainException(
			"run_ExtractAisTrajectories: undefined data type code "
			+ to_string(_c.TrajectoryExtraction.dataType) + ".");
	}
	vector<AisTrajectory<>> aisTrajectories = extractor->ConvertAll_fromRefinedAisCSV(refinedCSV);

	// save to files
	if (_c.TrajectoryExtraction.saveToFiles)
	{
		// lamda for key generation
		auto genkey = [](const AisTrajectory<>& traj) -> string 
		{
			return "imo" + to_string(traj.staticInfo.imo) + "_mmsi" + to_string(traj.staticInfo.mmsi); 
		};

		// initialize key count
		map<string, size_t> keycnt;
		for (auto& traj : aisTrajectories)		
			keycnt[genkey(traj)] = 0;
		
		// save to file
		#pragma omp parallel for
		for (size_t i = 0; i < aisTrajectories.size(); i++)
		{
			const AisTrajectory<>& traj = aisTrajectories[i];
			string key = genkey(traj);
			string filename = _c.TrajectoryExtraction.exportDirectory + "/" + key + "_" + to_string(keycnt[key]++) + ".json";
			IO::WriteFelineJsonFile<AisTrajectory<>>(filename, { traj });
		}
	}
	return aisTrajectories;
}

vector<AisTrajectory<>> import_Trajectories(const string dir)
{
	return IO::ReadFelineJsonFile<AisTrajectory<>>(dir);
}

vector<AisTrajectory<>>& run_TrajectorySimplification(vector<AisTrajectory<>>& trajs)
{
	// Trajectory Simplification with grounding evasion
	// define lambda task for rdp
	auto ReadAndRunRDP = [&](int Tidx, size_t idx)
	{		
		// prepare xtd estimation
		XTDEstimation xtdEst(Tidx);
		RasterGeoData::Ptr gData = SingletonRasterGeoData::GetInstance(Tidx).GetData("Bathymetry");

		ConditionedDouglasPeucker<XYVVaT, StaticType>::CondFunc cond = [&xtdEst, &gData](
			const vector<XYVVaT>& waypoints, size_t s, size_t e, size_t p, const StaticType& params) -> bool
		{
			Position tmps = (Position)waypoints[s];
			Position tmpe = (Position)waypoints[e];
			Position tmpp = (Position)waypoints[p];
			Real dist = abs(GeoDistance::CrossTrackDistance(tmpp, tmps, tmpe));

			if (dist > _c.TrajectorySimplification.epsilon) return true;

			auto p1 = gData->IndexAtPosition(tmps.dat, RasterGeoData::idxo_Closest);
			auto p2 = gData->IndexAtPosition(tmpe.dat, RasterGeoData::idxo_Closest);
			array<size_t, 2> pmin{ min(p1[0], p2[0]), min(p1[1], p2[1]) };
			array<size_t, 2> pmax{ max(p1[0], p2[0]), max(p1[1], p2[1]) };

			for (size_t i = pmin[0]; i <= pmax[0]; i++)
				for (size_t j = pmin[1]; j <= pmax[1]; j++)
				{
					if (gData->ValueAt(0, array<size_t,2>{ {i,j} }) > -(xtdEst.GetParams().DraughtMargin_D_coast + params.Dim.T))
					{
						return true;
					}
				}
			return false;
		};
		ConditionedDouglasPeucker<XYVVaT, StaticType> rdp(cond, trajs[idx].staticInfo);
		rdp.SimplifyRoute(trajs[idx].trajectory);
	};

	ThreadPool tp(_c.numOfThreads);
	for (size_t i = 0; trajs.size(); i++)
	{
		tp.Push(ReadAndRunRDP, i);
	}
	tp.Wait();

	// save to files
	if (_c.TrajectorySimplification.saveToFiles)
	{
		// lamda for key generation
		auto genkey = [](const AisTrajectory<>& traj) -> string
		{
			return "imo" + to_string(traj.staticInfo.imo) + "_mmsi" + to_string(traj.staticInfo.mmsi);
		};

		// initialize key count
		map<string, size_t> keycnt;
		for (auto& traj : trajs)
			keycnt[genkey(traj)] = 0;

		// save to file
		#pragma omp parallel for
		for (size_t i = 0; i < trajs.size(); i++)
		{
			const AisTrajectory<>& traj = trajs[i];
			string key = genkey(traj);
			string filename = _c.TrajectorySimplification.exportDirectory + "/" + key + "_" + to_string(keycnt[key]++) + ".json";
			IO::WriteFelineJsonFile<AisTrajectory<>>(filename, { traj });
		}
	}

	return trajs;
}

vector<AisTrajectory<XYVVaXtdTList>> run_GenerateXtd(vector<AisTrajectory<>> trajs)
{
	vector<AisTrajectory<XYVVaXtdTList>> re;	
	re.resize(trajs.size());

	#pragma omp parallel for
	for (size_t i = 0; i < trajs.size(); i++)
	{
		auto& traj = trajs[i];
		auto& reItem = re[i];

		XTDEstimation xtdEst(i);
		XYXtdList xtdRes = xtdEst.EstimateXTD(
			traj.trajectory.ToXYList(),
			traj.staticInfo.Dim.L, traj.staticInfo.Dim.B, traj.staticInfo.Dim.T);

		reItem.trajectory.resize(traj.trajectory.size());
		for (size_t j = 0; j < trajs[i].trajectory.size(); j++)
		{
			reItem.trajectory[j] = {
				traj.trajectory[j].Pos, traj.trajectory[j].Vel,
				xtdRes[j].Xtd, traj.trajectory[j].TP };
		}
	}

	// save to files
	if (_c.XtdGeneration.saveToFiles)
	{
		// lamda for key generation
		auto genkey = [](const AisTrajectory<XYVVaXtdTList>& traj) -> string
		{
			return "imo" + to_string(traj.staticInfo.imo) + "_mmsi" + to_string(traj.staticInfo.mmsi);
		};

		// initialize key count
		map<string, size_t> keycnt;
		for (auto& traj : re)
			keycnt[genkey(traj)] = 0;

		// save to file
		#pragma omp parallel for
		for (size_t i = 0; i < re.size(); i++)
		{
			const AisTrajectory<XYVVaXtdTList>& traj = re[i];
			string key = genkey(traj);
			string filename = _c.XtdGeneration.exportDirectory + "/" + key + "_" + to_string(keycnt[key]++) + ".json";
			IO::WriteFelineJsonFile<AisTrajectory<XYVVaXtdTList>>(filename, { traj });
		}
	}
	return re;
}

void run()
{
	// csv refinement
	vector<Table> refinedCSV;
	if (_c.runCsvRefinement)
	{
		refinedCSV = run_CsvRefinement();
	}

	// trajectory extraction
	vector<AisTrajectory<>> extractedTrajs;
	if (_c.runTrajectoryExtraction)
	{
		if (refinedCSV.size() <= 0)
		{
			refinedCSV = import_RefinedCsv(_c.TrajectoryExtraction.importDirectory);
		}
		extractedTrajs = run_TrajectoryExtraction(refinedCSV);
	}
	refinedCSV.clear();

	// trajectory simplification
	if (_c.runTrajectorySimplification)
	{
		if (extractedTrajs.size() <= 0)
		{
			extractedTrajs = import_Trajectories(_c.TrajectorySimplification.importDirectory);
		}
		run_TrajectorySimplification(extractedTrajs);
	}

	if (_c.runXtdGeneration)
	{
		if (extractedTrajs.size() <= 0)
		{
			extractedTrajs = import_Trajectories(_c.XtdGeneration.importDirectory);
		}
		run_TrajectorySimplification(extractedTrajs);
	}
}

void Initialize(int argc, char** argv)
{
	// Initialize core classes
	SingletonCLI::Initialize();	
	CommonLogger::Initialize("CommonLogger.json");

	// Initialize Feline Classes
	AisCsvReader::Initialize();
	AisTrajectoryExtraction::Initialize();
	AisTrajectoryExtraction_KMOF::Initialize();
	AisTrajectoryExtraction_ExactEarth::Initialize();	
	SingletonRasterGeoData::Initialize();	
	XTDEstimation::Initialize();

	// Initialize get options 
	CLI::App* cli = SingletonCLI::GetInstance().GetCLI();
	cli->add_option("--numOfThreads", _c.numOfThreads);
	cli->add_option("--runCsvRefinement", _c.runCsvRefinement);
	CLI::App* cli_csvref = SingletonCLI::GetInstance().GetCLI("CsvRefinement");
	{
		cli_csvref->add_option("--saveToFiles", _c.CsvRefinement.saveToFiles);
		cli_csvref->add_option("--importDirectory_dynamic", _c.CsvRefinement.importDirectory_dynamic);
		cli_csvref->add_option("--hasStaticFile", _c.CsvRefinement.hasStaticFile);
		cli_csvref->add_option("--importDirectory_static", _c.CsvRefinement.importDirectory_static);
		cli_csvref->add_option("--exportDirectory", _c.CsvRefinement.exportDirectory);
	}
	cli->add_option("--runTrajectoryExtraction", _c.runTrajectoryExtraction);
	CLI::App* cli_trajext = SingletonCLI::GetInstance().GetCLI("TrajectoryExtraction");
	{
		cli_trajext->add_option_function<string>("--AisDataType",
			[](const string& val)
			{
				if (val == "ExactEarth")
					_c.TrajectoryExtraction.dataType = ExactEarth;
				else if (val == "KMOF")
					_c.TrajectoryExtraction.dataType = KMOF;
				else
					throw mainException("Invalid option --AisDataType given. Should be either ExactEarth/KMOF.");				
			});
		cli_trajext->add_option("--saveToFiles", _c.TrajectoryExtraction.saveToFiles);
		cli_trajext->add_option("--importDirectory", _c.TrajectoryExtraction.importDirectory);				
		cli_trajext->add_option("--exportDirectory", _c.TrajectoryExtraction.exportDirectory);
	}
	cli->add_option("--runTrajectorySimplification", _c.runTrajectorySimplification);
	CLI::App* cli_trajabs = SingletonCLI::GetInstance().GetCLI("TrajectorySimplification");
	{
		cli_trajabs->add_option("--Epsilon", _c.TrajectorySimplification.epsilon);
		cli_trajabs->add_option("--saveToFiles", _c.TrajectorySimplification.saveToFiles);
		cli_trajabs->add_option("--importDirectory", _c.TrajectorySimplification.importDirectory);
		cli_trajabs->add_option("--exportDirectory", _c.TrajectorySimplification.exportDirectory);
	}
	cli->add_option("--runXtdGeneration", _c.runXtdGeneration);
	CLI::App* cli_xtdgen = SingletonCLI::GetInstance().GetCLI("XtdGeneration");
	{		
		cli_xtdgen->add_option("--saveToFiles", _c.XtdGeneration.saveToFiles);
		cli_xtdgen->add_option("--importDirectory", _c.XtdGeneration.importDirectory);
		cli_xtdgen->add_option("--exportDirectory", _c.XtdGeneration.exportDirectory);
	}
	
	// Parse additional options from command line.
	SingletonCLI::GetInstance().GetCLI()->set_help_all_flag("--help-all");

	cli->callback(run);
}

void Finalize()
{

}

int main(int argc, char** argv)
{
	Initialize(argc, argv);

	try
	{
		SingletonCLI::GetInstance().Parse(
			argc, argv,
			{
				//"./config/CommonLogger.json",
				"./RefineAIS.json"
			});
	}
	catch (const CLI::Error& e)
	{
		SingletonCLI::GetInstance().GetCLI()->exit(e);
	}
	catch (const HashColon::Exception& e)
	{
		CommonLogger logger;
		logger.Error({ { Tag::file, e.file()}, { Tag::func, e.func() }, { Tag::line, e.line() } }) << e.what() << endl;
		//string _file = e.file();
		//string _func = e.func();
		//int _line = e.line();
		//string _what = e.what();
		//logger.Error({ { Tag::file, _file }, { Tag::func, _func }, { Tag::line, _line } }) << _what << std::endl;
	}
	catch (const std::exception& e)
	{
		CommonLogger logger;
		logger.Error() << e.what() << endl;
	}
	catch (...)
	{
		CommonLogger logger;
		logger.Error() << "unknown error" << endl;
	}

	Finalize();
}
