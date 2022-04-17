//#include <HashColon/HashColon_config.h>
//#include <any>
//#include <cmath>
//#include <functional>
//#include <limits> 
//#include <memory>
//#include <string>
//#include <unordered_map>
//#include <vector>
//#include <Eigen/Eigen>
//#include <rapidjson/document.h>
//#include <rapidjson/istreamwrapper.h>
//#include <rapidjson/ostreamwrapper.h>
//#include <rapidjson/prettywriter.h>
//#include <rapidjson/writer.h>
//#include <HashColon/Clustering.hpp>
//#include <HashColon/Exception.hpp>
//#include <HashColon/Helper.hpp>
//#include <HashColon/Log.hpp>
//#include <HashColon/Real.hpp>
//#include <HashColon/SingletonCLI.hpp>
//#include <HashColon/Feline/FelineJsonIO.hpp>
//#include <HashColon/Feline/TrajectoryClustering.hpp>
//#include <HashColon/Feline/XtdTrajectoryClustering.hpp>
//
//using namespace std;
//using namespace Eigen;
//using namespace HashColon;
//using namespace HashColon::Fs;
//using namespace HashColon::Clustering;
//using namespace HashColon::Feline;
//
//using Tag = HashColon::LogUtils::Tag;
//namespace Json = rapidjson;
//
//void run()
//{
//	XYList testlist;
//	testlist.push_back({ 127.480433, 33.661609 });
//	testlist.push_back({ 127.480433, 34.661609 });
//	testlist.push_back({ 128.480433, 34.661609 });
//	testlist.push_back({ 128.480433, 33.661609 });
//
//	XY testA{ 127, 33 };
//	XY testB = testA.MoveTo(100000, 0);
//
//	cout << testB.longitude << ", " << testB.latitude << endl;
//	
//
//	//XYList unisampled = testlist.GetUniformLengthSampled(31);
//
//	//vector<XYList> output;
//	//output.push_back(testlist);
//	//output.push_back(unisampled);
//		
//	//IO::WriteGeoJsonFile("/home/cadit/WTK/testresult.json", output, true);	
//}
//
//static void Initialize()
//{
//	// register config files	
//	SingletonCLI::Initialize();
//	CommonLogger::Initialize("./CommonLogger.json");
//
//	GeoDistance::Initialize();
//
//	NJW<XYList>::Initialize();
//	DistanceBasedDBSCAN<XYList>::Initialize();
//	TrajectoryClustering::Initialize_All_TrajectoryDistanceMeasure();
//
//	NJW<XYXtdList>::Initialize();
//	DistanceBasedDBSCAN<XYXtdList>::Initialize();
//	XtdTrajectoryClustering::Initialize_All_XtdTrajectoryDistanceMeasure();
//
//	CLI::App* cli = SingletonCLI::GetInstance().GetCLI();
//	
//	// Parse additional options from command line.
//	SingletonCLI::GetInstance().GetCLI()->set_help_all_flag("--help-all");
//
//	cli->callback(run);
//}
//
//int main(int argc, char** argv)
//{
//	Initialize();
//
//	try
//	{
//		SingletonCLI::GetInstance().Parse(
//			argc, argv,
//			{
//				//"./config/CommonLogger.json",
//				"./TrajectoryClustering.json"
//			});
//	}
//	catch (const CLI::Error& e)
//	{
//		SingletonCLI::GetInstance().GetCLI()->exit(e);
//	}
//	catch (const HashColon::Exception& e)
//	{
//		CommonLogger logger;
//		logger.Error({ { Tag::file, e.file()}, { Tag::func, e.func() }, { Tag::line, e.line() } }) << e.what() << std::endl;
//	}
//
//}