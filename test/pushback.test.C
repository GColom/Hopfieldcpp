#include"../source/HopfieldNetwork.h"
#include<cstdlib>
int main()
{
	auto H = HopfieldNetwork(10, 0.0, 0.0);

	auto p = spin_pattern(5, +1);

	try{
		H.push_back_pattern(p);
	}
	catch(std::runtime_error)
	{
		return 255;
	}
}