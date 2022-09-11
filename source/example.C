#include"HopfieldNetwork.h"
#include<iostream>
#include<chrono>

int main()
{
	auto start = std::chrono::steady_clock::now();

	HopfieldNetwork H(10000, 0.0125, 0.0);
	std::cout<<"Network instantiated."<<std::endl;
	H.build_random_patterns();
/*
for (int step = 5; step > 0; --step)
{
	double cur_alpha = step*0.03;
	std::cout<<"Ramping alpha to "<<cur_alpha<<std::endl;
	H.set_alpha(cur_alpha);
	std::cout<<"Building weights..."<<std::endl;
	H.build_weights();
	std::cout<<"Init corrupted pattern..."<<std::endl;
	H.init_on_corrupted_pattern(0, 0.3);
	std::cout<<"Evolving system..."<<std::endl;
	H.metropolis_evolve(1000000);
	auto ov = H.max_overlap();
	std::cout<<"argmax, max"<<std::endl;
	std::cout<<ov.first<<", "<<ov.second<<std::endl;
}*/

for (int step = 5; step > 0; --step)
{
	double cur_T = step*0.2;
	std::cout<<"Setting T to "<<cur_T<<std::endl;
	H.set_temperature(cur_T);
	std::cout<<"Building weights..."<<std::endl;
	H.build_weights();
	std::cout<<"Init corrupted pattern..."<<std::endl;
	H.init_on_corrupted_pattern(0, 0.3);
	std::cout<<"Evolving system..."<<std::endl;
	H.glauber_evolve(1000000);
	auto ov = H.max_overlap();
	std::cout<<"argmax, max"<<std::endl;
	std::cout<<ov.first<<", "<<ov.second<<std::endl;
}
	auto finish = std::chrono::steady_clock::now();

	std::cout<<std::chrono::duration<double>(finish - start).count()<<std::endl;
}