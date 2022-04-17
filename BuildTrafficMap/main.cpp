/*
* Input: Vessel/Time sorted AIS data in CSV format.
*		 Assume that AIS data from a vessel does not exist in different CSV files.
* Output: Maritime Traffic map in GeoJson(for visualize) (or .shp file)
* 
*/
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
#include <HashColon/Feline/AisPreprocess.hpp>
#include <HashColon/Feline/FelineJsonIO.hpp>
#include <HashColon/Feline/GeoData.hpp>
#include <HashColon/Feline/RouteSimplification.hpp>
#include <HashColon/Feline/XtdEstimation.hpp>
#include <HashColon/Feline/GeoValues.hpp>

using namespace std;
namespace Json = rapidjson;
using namespace HashColon;
using namespace HashColon::Feline;
using namespace HashColon::LogUtils;

struct _Params
{	
	vector<string> felineJsonPaths;
	int numOfThreads = 1;
};

_Params _c;

void run()
{
	CommonLogger logger;
	logger.Log({ {Tag::lvl, 1} }) << "Maritime Traffic Map building started." << endl;

	// Gather file names
	vector<string> felineJsonFilenames = Fs::GetFilesFromPaths(_c.felineJsonPaths, ".*\\.json");

	// define mutex
	mutex _mutex_stdout;
	mutex _mutex_result;	
	atomic<size_t> progressCnt{ 0 };
	vector<AisTrajectory<>> resultTrajs;

	// define lambda task for rdp
	auto ReadAndRunRDP = [&](int Tidx, size_t idx)
	{
		// read feline json file
		vector<AisTrajectory<>> trajs = IO::ReadFelineJsonFile<AisTrajectory<>>(felineJsonFilenames[idx]);

		// prepare xtd estimation
		XTDEstimation xtdEst(Tidx);
		RasterGeoData::Ptr gData = SingletonRasterGeoData::GetInstance(Tidx).GetData("Bathymetry");
				
		ConditionedDouglasPeucker<XYVVaT, StaticType>::CondFunc cond = [&xtdEst, &gData](
			const vector<XYVVaT>& waypoints, size_t s, size_t e, size_t p, const StaticType& params) -> bool
		{
			Position tmps = (Position)waypoints[s];
			Position tmpe = (Position)waypoints[e];
			Position tmpp = (Position)waypoints[p];
			Real dist = GeoDistance::CrossTrackDistance(tmpp, tmps, tmpe);
			Degree angle = (s > 0 && (s < waypoints.size() - 1))
				? GeoDistance::Angle(tmps, (Position)waypoints[s + 1]) - GeoDistance::Angle((Position)waypoints[s - 1], tmps)
				: 0.0;

			XTD xtd = xtdEst.EstimateXTD(tmps, (Position)waypoints[s + 1], params.Dim.L, params.Dim.B, params.Dim.T, angle);

			if ((dist > 0 && dist > xtd.xtdPortside) && (dist < 0 && dist < xtd.xtdStarboard))
			{
				return true;
			}

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

		for (AisTrajectory traj : trajs)
		{
			ConditionedDouglasPeucker<XYVVaT, StaticType> rdp(cond, traj.staticInfo);
			rdp.SimplifyRoute(traj.trajectory);			
		}
		{
			lock_guard<mutex> _lg(_mutex_result);
			resultTrajs.insert(resultTrajs.end(), trajs.begin(), trajs.end());
		}
		{
			lock_guard<mutex> _lg(_mutex_stdout);
			stringstream tempss;
			tempss << "Trajectory simplification: " << progressCnt << "/" << felineJsonFilenames.size() 
				<< " " << Percentage(++progressCnt, (int)felineJsonFilenames.size());
			logger.Message << Flashl(tempss.str());
		}
	};

	// run rdp as async task
	ThreadPool tp(_c.numOfThreads);
	for (size_t i = 0; i < felineJsonFilenames.size(); i++)
	{
		tp.Push(ReadAndRunRDP, i);		
	}
	tp.Wait();
	logger.Log({ {Tag::lvl, 1} }) << "Trajectory simplification finished." << endl;

	

	// todo:  clustering....
	

}

void Initialize(int argc, char** argv)
{
	//// Initialize helper classes
	//SingletonCLI::Initialize();
	//CommonLogger::Initialize("CommonLogger.json");

	//// Initialize Feline Classes
	//SingletonRasterGeoData::Initialize("./RefineAIS.json");
	//XTDEstimation::Initialize();

	//// Initialize get options 
	//CLI::App* cli = SingletonCLI::GetInstance().GetCLI();
	//cli->add_option("--inputDir", _c.inputDir);
	//cli->add_option("--outputDir", _c.outputDir);
	//cli->add_option("--shiplistPath", _c.shiplistPath);
	//cli->add_option("--staticDataPath", _c.staticDataPath);
	//cli->add_option("--NumOfThreads", _c.NumOfThreads);
	//cli->add_option("--RouteSimplificationEpsilon", _c.RouteSimplificationEpsilon);

	//// Parse additional options from command line.
	//SingletonCLI::GetInstance().GetCLI()->set_help_all_flag("--help-all");

	//cli->callback(run);
}

int main(int argc, char** argv)
{
	//Initialize(argc, argv);

	//try
	//{
	//	SingletonCLI::GetInstance().Parse(
	//		argc, argv,
	//		{
	//			//"./config/CommonLogger.json",
	//			"./RefineAIS.json"
	//		});
	//}
	//catch (const CLI::Error& e)
	//{
	//	SingletonCLI::GetInstance().GetCLI()->exit(e);
	//}
	//catch (const HashColon::Exception& e)
	//{
	//	CommonLogger logger;
	//	logger.Error({ { Tag::file, e.file()}, { Tag::func, e.func() }, { Tag::line, e.line() } }) << e.what() << endl;
	//	//string _file = e.file();
	//	//string _func = e.func();
	//	//int _line = e.line();
	//	//string _what = e.what();
	//	//logger.Error({ { Tag::file, _file }, { Tag::func, _func }, { Tag::line, _line } }) << _what << std::endl;
	//}
	//catch (const std::exception& e)
	//{
	//	CommonLogger logger;
	//	logger.Error() << e.what() << endl;
	//}
	//catch (...)
	//{
	//	CommonLogger logger;
	//	logger.Error() << "unknown error" << endl;
	//}

	//Finalize();
}