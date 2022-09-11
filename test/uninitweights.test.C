#include"../source/HopfieldNetwork.h"
#include<cstdlib>
int main()
{
	auto H = HopfieldNetwork(10, 0.1, 0.0);
	H.build_random_patterns();

	try{
		H.glauber_evolve(10);
	}
	catch(std::runtime_error)
	{
		return 255;
	}
}