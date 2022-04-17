#define DISTANCE_TEST_EXE

#ifdef DISTANCE_TEST_EXE

#include <HashColon/HashColon_config.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <HashColon/Log.hpp>
#include <HashColon/Real.hpp>
#include <HashColon/Feline/GeoValues.hpp>
#include <HashColon/Feline/XtdTrajectoryClustering.hpp>

using namespace std;
using namespace HashColon;
using namespace HashColon::Feline;
using namespace HashColon::Feline::XtdTrajectoryClustering;

using Tag = HashColon::LogUtils::Tag;
namespace Json = rapidjson;

struct _Params 
{
	struct { 
		size_t range; 
		vector<Real> angles;
		string resultDir;
	} distExp;

	struct {
		vector<Real> distances;
		string resultDir;
	} angleExp;
	
	struct {
		Position pos;
		XTD xtd;
		Degree heading;
	} A_Info;
	
	struct {		
		XTD xtd;		
	} B_Info;

	_JSDivergenceOption JsOpt;	
	_WassersteinOption WsOpt;
} _c;

void DistExp_JS()
{
	CommonLogger logger;
	logger.Log() << "Starting [Distance varying experiment] with [JS metric]." << endl;

	vector<Real> re; re.resize(_c.distExp.range + 1);
	XYXtd A = { _c.A_Info.pos, _c.A_Info.xtd };
	Degree A_Dir = _c.A_Info.heading;

	for (size_t i = 0; i < _c.distExp.angles.size(); i++)
	{
		#pragma omp parallel for
		for (size_t d = 0; d <= _c.distExp.range; d++)
		{
			XYXtd B{ A.Pos.MoveTo((Real)d, 0), _c.B_Info.xtd };
			Degree B_Dir = (Real)_c.distExp.angles[i];

			re[d] = JSDivergenceDistance(A, A_Dir, B, B_Dir, _c.JsOpt);	
		}
		logger.Log() << "Exporting results of [Distance varying experiment] #" + to_string(i) + " with [JS metric]." << endl;
		string fulloutputpath = _c.distExp.resultDir + "/" + "distTest_JS_" + to_string(i) + ".txt";
		ofstream ofs(fulloutputpath);
		for (size_t d = 0; d <= _c.distExp.range; d++)
		{
			ofs << d << "\t" << re[d] << "\n";
		}
		ofs << flush;
	}
	logger.Log() << "Finished [Distance varying experiment] with [JS metric]." << endl;
}

void DistExp_EMD()
{
	CommonLogger logger;
	logger.Log() << "Starting [Distance varying experiment] with [EMD metric]." << endl;

	vector<Real> re; re.resize(_c.distExp.range + 1);
	XYXtd A = { _c.A_Info.pos, _c.A_Info.xtd };
	Degree A_Dir = _c.A_Info.heading;

	for (size_t i = 0; i < _c.distExp.angles.size(); i++)
	{
		Degree B_Dir = (Real)_c.distExp.angles[i];

		#pragma omp parallel for
		for (size_t d = 0; d <= _c.distExp.range; d++)
		{
			XYXtd B{ A.Pos.MoveTo((Real)d, 0), _c.B_Info.xtd };			

			re[d] = WassersteinDistance(A, A_Dir, B, B_Dir, _c.WsOpt);			
		}
		logger.Log() << "Exporting results of [Distance varying experiment] #" + to_string(i) + " with [EMD metric]." << endl;
		string fulloutputpath = _c.distExp.resultDir + "/" + "distTest_EMD_" + to_string(i) + ".txt";
		ofstream ofs(fulloutputpath);
		for (size_t d = 0; d <= _c.distExp.range; d++)
		{
			ofs << d << "\t" << re[d] << "\n";
		}
		ofs << flush;
	}
	logger.Log() << "Finished [Distance varying experiment] with [END metric]." << endl;
}

void AngleExp_JS()
{
	CommonLogger logger;
	logger.Log() << "Starting [Heading varying experiment] with [JS metric]." << endl;

	vector<Real> re; re.resize(360);
	XYXtd A = { _c.A_Info.pos, _c.A_Info.xtd };
	Degree A_Dir = _c.A_Info.heading;

	for (size_t i = 0; i < _c.angleExp.distances.size(); i++)
	{
		XYXtd B{ A.Pos.MoveTo((Real)_c.angleExp.distances[i], 0), _c.B_Info.xtd };

		#pragma omp parallel for
		for (size_t d = 0; d < 360; d++)
		{			
			Degree B_Dir = (Real)d;
			re[d] = JSDivergenceDistance(A, A_Dir, B, B_Dir, _c.JsOpt);
		}
		logger.Log() << "Exporting results of [Heading varying experiment] #" + to_string(i) + " with [JS metric]." << endl;
		string fulloutputpath = _c.angleExp.resultDir + "/" + "angleTest_JS_" + to_string(i) + ".txt";
		ofstream ofs(fulloutputpath);
		for (size_t d = 0; d < 360; d++)
		{
			ofs << d << "\t" << re[d] << "\n";
		}
		ofs << flush;
	}
	logger.Log() << "Finished [Heading varying experiment] with [JS metric]." << endl;
}

void AngleExp_EMD()
{
	CommonLogger logger;
	logger.Log() << "Starting [Heading varying experiment] with [EMD metric]." << endl;

	vector<Real> re; re.resize(360);
	XYXtd A = { _c.A_Info.pos, _c.A_Info.xtd };
	Degree A_Dir = _c.A_Info.heading;

	for (size_t i = 0; i < _c.angleExp.distances.size(); i++)
	{
		XYXtd B{ A.Pos.MoveTo((Real)_c.angleExp.distances[i], 0), _c.B_Info.xtd };

#pragma omp parallel for
		for (size_t d = 0; d < 360; d++)
		{
			Degree B_Dir = (Real)d;
			re[d] = WassersteinDistance(A, A_Dir, B, B_Dir, _c.WsOpt);
		}
		logger.Log() << "Exporting results of [Heading varying experiment] #" + to_string(i) + " with [EMD metric]." << endl;
		string fulloutputpath = _c.angleExp.resultDir + "/" + "angleTest_EMD_" + to_string(i) + ".txt";
		ofstream ofs(fulloutputpath);
		for (size_t d = 0; d < 360; d++)
		{
			ofs << d << "\t" << re[d] << "\n";
		}
		ofs << flush;
	}
	logger.Log() << "Finished [Heading varying experiment] with [EMD metric]." << endl;
}


void run()
{
	CommonLogger logger;
	logger.Log({ {Tag::lvl, 1} }) << "Options: "
		<< SingletonCLI::GetInstance().GetCLI()->config_to_str() << endl;
	DistExp_JS();
	AngleExp_JS();
	DistExp_EMD();
	AngleExp_EMD();
}

static void Initialize()
{
	// register config files	
	SingletonCLI::Initialize();
	CommonLogger::Initialize("./CommonLogger.json");

	GeoDistance::Initialize();	
	
	CLI::App* cli_distExp = SingletonCLI::GetInstance().GetCLI("distExp");
	cli_distExp->add_option("--range", _c.distExp.range);
	cli_distExp->add_option("--angles", _c.distExp.angles);
	cli_distExp->add_option("--resultDir", _c.distExp.resultDir);

	CLI::App* cli_angleExp = SingletonCLI::GetInstance().GetCLI("angleExp");
	cli_angleExp->add_option("--distances", _c.angleExp.distances);
	cli_angleExp->add_option("--resultDir", _c.angleExp.resultDir);

	CLI::App* cli_A = SingletonCLI::GetInstance().GetCLI("A_Info");
	cli_A->add_option("--lat", _c.A_Info.pos.latitude);
	cli_A->add_option("--lon", _c.A_Info.pos.longitude);
	cli_A->add_option("--xtdp", _c.A_Info.xtd.xtdPortside);
	cli_A->add_option("--xtds", _c.A_Info.xtd.xtdStarboard);
	cli_A->add_option("--headings", _c.A_Info.heading);

	CLI::App* cli_B = SingletonCLI::GetInstance().GetCLI("B_Info");
	cli_B->add_option("--xtdp", _c.B_Info.xtd.xtdPortside);
	cli_B->add_option("--xtds", _c.B_Info.xtd.xtdStarboard);
	
	CLI::App* cli_jsopt = SingletonCLI::GetInstance().GetCLI("JsOpt");
	cli_jsopt->add_option("--domainUnit", _c.JsOpt.domainUnit);
	cli_jsopt->add_option("--domainSize", _c.JsOpt.domainSize);
	cli_jsopt->add_option("--epsilon", _c.JsOpt.errorEpsilon);

	CLI::App* cli_wsopt = SingletonCLI::GetInstance().GetCLI("WsOpt");
	cli_wsopt->add_option("--domainUnit", _c.WsOpt.domainUnit);
	cli_wsopt->add_option("--domainSize", _c.WsOpt.domainSize);
	cli_wsopt->add_option("--epsilon", _c.WsOpt.errorEpsilon);

	// Parse additional options from command line.
	SingletonCLI::GetInstance().GetCLI()->set_help_all_flag("--help-all");

	CLI::App* cli = SingletonCLI::GetInstance().GetCLI();
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
				"DistanceTest.json"
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


#endif