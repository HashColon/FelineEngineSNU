//#define EMD_TEST_EXE

#ifdef EMD_TEST_EXE
#include <HashColon/HashColon_config.h>
#include <Wasserstein.hh>
#include <iostream>

int main()
{
	using namespace std;
	using namespace emd;
	using EMDParticle = EuclideanParticle2D<>;
	using EMD = EMDFloat64<EuclideanEvent2D, EuclideanDistance2D>;
	using PairwiseEMD = PairwiseEMD<EMD>;

	//PairwiseEMD emdObj(0.4, 1, true);
	EuclideanEvent2D<> ev0;
	//ev0.
	cout << "1:" << EMDParticle::dimension() << endl;
	cout << "2:" << EMDParticle::distance_name() << endl;
	cout << "3:" << EMDParticle::name() << endl;
	
	
	

	//cout << "3:" << EMDParticle::plain_distance( << endl;
	
	//emdObj.

}

#endif
