#include"HopfieldNetwork.h"
#include<iostream>
#include<chrono>

int main()
{
	auto start = std::chrono::steady_clock::now();

	HopfieldNetwork H(10000, 0.05, 0.01);

	H.build_random_patterns();

	H.build_weights();

	H.init_on_corrupted_pattern(0, 0.2);

	H.glauber_evolve(10000);

	std::cout<<H.overlap_max()<<std::endl;
	std::cout<<H.overlap_argmax()<<std::endl;

	auto finish = std::chrono::steady_clock::now();

	std::cout<<std::chrono::duration<double>(finish - start).count()<<std::endl;
}