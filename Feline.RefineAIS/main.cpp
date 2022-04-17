#include <HashColon/HashColon_config.h>
#include <algorithm>
#include <mutex>

#include <rapidjson/document.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/prettywriter.h>

#include <HashColon/Core/Exception.hpp>
#include <HashColon/Core/SingletonCLI.hpp>
#include <HashColon/Core/Log.hpp>
#include <HashColon/Core/ThreadPool.hpp>

#include <HashColon/Core/Helper.hpp>
#include <HashColon/Helper/CsvHelper.hpp>

#include <HashColon/Feline/GeoValues.hpp>
#include <HashColon/Feline/Data/GdalSupport.hpp>
#include <HashColon/Feline/TrajectoryXtdEstimation.hpp>
#include <HashColon/Feline/RouteSimplification.hpp>
#include <HashColon/Feline/FelineJsonIO.hpp>

#include <iostream>

using namespace std;
namespace Json = rapidjson;
using namespace HashColon;
using namespace HashColon::Helper;
using namespace HashColon::Feline;
using Tag = HashColon::LogUtils::Tag;

HASHCOLON_NAMED_EXCEPTION_DEFINITION(main);
mutex _mutex;

struct _Params
{
	string inputDir;
	string outputDir;
	string shiplistPath;
	string staticDataPath;
	size_t NumOfThreads;
	Real RouteSimplificationEpsilon;
} _c;

struct StaticDataType
{
	int mmsi;
	int imo;
	Real loa;
	Real beam;
	Real draught;
};

CSVRow SearchByMmsiImo(HashColonCSV staticdata, int mmsi, int imo)
{
	auto checkMmsiImo = [mmsi, imo](CSVRow row) { return row["mmsi"].Get<int>() == mmsi && row["imo"].Get<int>() == imo; };

	{
		lock_guard<mutex> lg(_mutex);
		auto re = find_if(staticdata.begin(), staticdata.end(), checkMmsiImo);
		if (re != staticdata.end()) return (*re);
		else
			throw mainException("mmsi/imo not found");
	}
}

void WriteJsonFile(const string filepath, StaticDataType sdata, vector<XYXtdList> routes)
{
	Json::Document doc;
	Json::Document::AllocatorType& alloc = doc.GetAllocator();

	doc.Parse(R"str(
	{
		"routes" : []
	}
	)str");

	for (size_t i = 0; i < routes.size(); i++)
	{
		// route value
		Json::Value ObjRoute(Json::kObjectType);		

		// set static value;
		Json::Value ObjStatic(Json::kObjectType);		
		ObjStatic.AddMember("mmsi", Json::Value(sdata.mmsi), alloc);
		ObjStatic.AddMember("imo", Json::Value(sdata.imo), alloc);
		ObjStatic.AddMember("loa", Json::Value(sdata.loa), alloc);
		ObjStatic.AddMember("beam", Json::Value(sdata.beam), alloc);
		ObjStatic.AddMember("draught", Json::Value(sdata.draught), alloc);
		ObjRoute.AddMember("static", ObjStatic, alloc);

		// set dynamic values
		Json::Value ObjDynamic; ObjDynamic.SetArray();
		for (size_t j = 0; j < routes[i].size(); j++)
		{
			Json::Value ObjWP(Json::kObjectType);
			ObjWP.AddMember("latitude", Json::Value(routes[i][j].Pos.latitude), alloc);
			ObjWP.AddMember("longitude", Json::Value(routes[i][j].Pos.longitude), alloc);
			ObjWP.AddMember("xtd_p", Json::Value(routes[i][j].Xtd.xtdPortside), alloc);
			ObjWP.AddMember("xtd_s", Json::Value(routes[i][j].Xtd.xtdStarboard), alloc);
			ObjDynamic.PushBack(ObjWP, alloc);
		}
		ObjRoute.AddMember("dynamic", ObjDynamic, alloc);

		// add route
		doc["routes"].PushBack(ObjRoute, alloc);
	}	

	ofstream ofs(filepath);
	Json::OStreamWrapper osw(ofs);
	Json::PrettyWriter<Json::OStreamWrapper> writer(osw);
	doc.Accept(writer);
}


void run()
{
	CommonLogger logger;
	logger.Log({ {Tag::lvl, 1} }) << "AIS Refinement started." << endl;

	// defines an item for AIS Refinement
	struct itemtype 
	{
		int mmsi;
		int imo;
		string filepath;
	};
	vector<itemtype> items; items.clear();

	// read shiplist file
	HashColonCSV shiplist; shiplist.ReadFromFile(_c.shiplistPath);

	// read static informatino file
	HashColonCSV staticData; staticData.ReadFromFile(_c.staticDataPath); 

	//// get file lists	
	//size_t debug = 1;
	//size_t debug2 = 0;
	//#pragma omp parallel for num_threads(20)
	//for(size_t sl_idx = 0 ; sl_idx < shiplist.size(); sl_idx++)
	////for (auto shiplistrow : shiplist)
	//{
	//	auto& shiplistrow = shiplist[sl_idx];

	//	int mmsi_int = shiplistrow["mmsi"].Get<int>();
	//	int imo_int = shiplistrow["imo"].Get<int>();
	//	string mmsi_str = shiplistrow["mmsi"].Get<string>();
	//	string imo_str = shiplistrow["imo"].Get<string>();

	//	vector<string> fileList = GetFilesInDirectory(_c.inputDir, {}, 
	//		"mmsi" + mmsi_str + "_imo" + imo_str + "_(\\d*)\\.json");

	//	for (auto afile : fileList)
	//	{
	//		lock_guard<mutex> _lg(_mutex);
	//		items.push_back({ mmsi_int, imo_int, afile });

	//		/*if (debug * 100 < items.size())
	//		{
	//			cout << debug * 100 << ", ";
	//			if (debug % 10 == 0) cout << endl;
	//			debug++;
	//		}*/
	//	}

	//	{
	//		lock_guard<mutex> _lg(_mutex);
	//		if (debug2 % 500 == 0) cout << debug2 << endl;
	//		debug2++;
	//	}
	//	
	//}

	vector<string> fileList = GetFilesInDirectory(_c.inputDir);
	//#pragma omp parallel for num_threads(20)
	for (size_t fileIdx = 0; fileIdx < fileList.size(); fileIdx++)
	{
		auto& afile = fileList[fileIdx];
		vector<string> extracted1 = Split(Split(afile, "/").back(), "msi_io.jsn");
		vector<string> extracted2;
		for (size_t tmpIdx = 0; tmpIdx < extracted1.size(); tmpIdx++)
			if (!extracted1[tmpIdx].empty())
				extracted2.push_back(extracted1[tmpIdx]);

		int mmsi_int = atoi(extracted2[0].c_str());
		int imo_int = atoi(extracted2[1].c_str());
		{
			lock_guard<mutex> _lg(_mutex);
			items.push_back({ mmsi_int, imo_int, afile });
		}
	}

	logger.Log({ { Tag::lvl, 3 } }) << items.size() << " files imported." << endl;

	// for each file in file list
	auto work = [&staticData](int Tidx, itemtype item)
	{
		using namespace HashColon::Feline;
		CommonLogger logger;
		{
			lock_guard<mutex> lg(_mutex);
			logger.Log({ {Tag::lvl, 1} }) << item.filepath << " started." << endl;
		}

		// read routes from a file
		vector<XYList> routes = IO::ReadFelineJsonFile<XYList>(item.filepath);
		vector<XYXtdList> routesWithXTD;

		// get static data
		CSVRow staticrow = SearchByMmsiImo(staticData, item.mmsi, item.imo);
		StaticDataType staticItem{
			item.mmsi, item.imo,
			staticrow["dimension_a"].Get<Real>() + staticrow["dimension_b"].Get<Real>(),
			staticrow["dimension_c"].Get<Real>() + staticrow["dimension_d"].Get<Real>(),
			staticrow["draught"].Get<Real>() / 10.0
		};

		// XTD estimation prep
		XTDEstimation xtdest(Tidx);

		// Route Simplification prep. 
		/*Data::RasterGeoData::Ptr gdata = make_shared<Data::RasterGeoData>(
			Data::RasterGeoData("/home/cadit/WTK/WTK_data/GeoData/gebco_2020_n45.0_s30.0_w120.0_e135.0.nc")
			);*/
		Data::RasterGeoData::Ptr gdata = Data::SingletonRasterGeoData::GetInstance(Tidx).GetData("Bathymetry");
		//gdata->BuildCache();
				
		auto pathValidattion = [&gdata, &xtdest](const Position& a, const Position& b, void* additionals) -> bool
		{
			auto p1 = gdata->IndexAtPosition(a.dat, Data::RasterGeoData::idxo_Closest);
			auto p2 = gdata->IndexAtPosition(b.dat, Data::RasterGeoData::idxo_Closest);

			array<size_t, 2> pmin{ min(p1[0], p2[0]), min(p1[1], p2[1]) };
			array<size_t, 2> pmax{ max(p1[0], p2[0]), max(p1[1], p2[1]) };
			StaticDataType adds = *((StaticDataType*)additionals);

			for (size_t i = pmin[0]; i <= pmax[0]; i++)
			{
				for (size_t j = pmin[1]; j <= pmax[1]; j++)
				{
					if (gdata->ValueAt(0, array<size_t, 2>{i, j}) > -(xtdest.GetParams().DraughtMargin_D_coast+adds.draught))
					{
						return false;
					}
				}
			}
			return true;
		};

		// simplify route & estimate xtds
		//for (XYList route : routes)
		for(size_t routeIdx = 0 ; routeIdx < routes.size(); routeIdx++)
		{			
			try 
			{
				auto& route = routes[routeIdx];
				SimplifyRoute<Position, SimplifyRouteMethods::RDP>(
					route, SimplifyRouteMethods::RDP{ _c.RouteSimplificationEpsilon }, pathValidattion, &staticItem);
				XYXtdList tmpXtdRoute = xtdest.EstimateXTD(route, staticItem.loa, staticItem.beam, staticItem.draught);
				routesWithXTD.push_back(tmpXtdRoute);
			}
			catch (const HashColon::Exception& e)
			{
				lock_guard<mutex> _lg(_mutex);
				CommonLogger logger;
				logger.Error({ { Tag::file, e.file()}, { Tag::func, e.func() }, { Tag::line, e.line() } }) << e.what() << endl;
			}
		}

		// write as json file (do not write this function in namespace IO)
		string filename = _c.outputDir + "/" + Split(item.filepath, "/").back();
		if (routesWithXTD.size() > 0)
		{			
			WriteJsonFile(filename, staticItem, routesWithXTD);
			{
				lock_guard<mutex> lg(_mutex);
				logger.Log({ {Tag::lvl, 1} }) << filename << " finished. (" <<  routesWithXTD.size() << " routes exported)" << endl;
			}
		}
		else
		{
			lock_guard<mutex> lg(_mutex);
			logger.Log({ {Tag::lvl, 1} }) << filename << " ignored. No routes." << endl;
		}
			
		
	};

	// Define a thread pool
	ThreadPool tp((int)_c.NumOfThreads);	
	for (auto& item : items) {		
		tp.Push(work, item); 		
		/*if (item.mmsi == 431005341 && item.imo == 9713002)
		{
			tp.Push(work, item);
			break;
		}*/
	}
	tp.Wait();	// wait for all threads to be finished.
	//work(0, items[0]);
	logger.Log({ {Tag::lvl, 1} }) << "AIS Refinement finished." << endl;
	
}

void Initialize(int argc, char** argv)
{
	// Initialize helper classes
	SingletonCLI::Initialize();
	CommonLogger::Initialize("CommonLogger.json");

	// Initialize Feline Classes
	Data::SingletonRasterGeoData::Initialize("./RefineAIS.json");
	XTDEstimation::Initialize();

	// Initialize get options 
	CLI::App* cli = SingletonCLI::GetInstance().GetCLI();
	cli->add_option("--inputDir", _c.inputDir);
	cli->add_option("--outputDir", _c.outputDir);
	cli->add_option("--shiplistPath", _c.shiplistPath);
	cli->add_option("--staticDataPath", _c.staticDataPath);
	cli->add_option("--NumOfThreads", _c.NumOfThreads);
	cli->add_option("--RouteSimplificationEpsilon", _c.RouteSimplificationEpsilon);
		
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