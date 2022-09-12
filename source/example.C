#include"HopfieldNetwork.h"
#include<iostream>
#include<fstream>
#include<chrono>
#include<numeric>
#include<string>
#include<iomanip>

static int precision = 4;

double mean(std::vector<double> data)
{
	if(data.size() <=1) throw(std::runtime_error("Cannot calculate the mean of less than two values."));

	return std::accumulate(data.begin(), data.end(), 0.0)/data.size();
}

double std_dev(std::vector<double> data)
{
	if(data.size() <=2) throw(std::runtime_error("Cannot calculate the standard deviation of less than three values."));
	double mean = std::accumulate(data.begin(), data.end(), 0.0)/data.size();

	double res = 0;
	std::for_each (data.begin(), data.end(), [&](const double d) {res += (d - mean) * (d - mean);});

	return std::sqrt(res/(data.size()-1));
}

int main(int argc, char ** argv)
{
	if(argc < 8) throw(std::runtime_error("Not enough parameters: calling signature is:\n ./example N_spins alpha T_lower_bound T_upper_bound T_steps measurements_per_step iterations_per_measurement output_filepath (optional)"));

	int N_spins 	= std::stoi(argv[1]);
	double alpha 	= std::stod(argv[2]);
	double low_T	= std::stod(argv[3]);
	double high_T	= std::stod(argv[4]);
	int T_steps		= std::stoi(argv[5]);
	int repetitions = std::stoi(argv[6]);
	int n_iter		= std::stoi(argv[7]);

	std::ofstream o;
	std::string outpath;

	if (argc == 9)
	{
		outpath = std::string(argv[8]);
		o.open(outpath);
	}

	auto start 		= std::chrono::steady_clock::now();

	HopfieldNetwork H(N_spins, alpha, 0.0);
	std::cout<<"Network instantiated."<<std::endl;
	std::cout<<"N = "<<H.get_N()<<std::endl;
	std::cout<<"alpha = "<<H.get_alpha()<<std::endl;
	std::cout<<"Temperature lower bound = "<<low_T<<std::endl;
	std::cout<<"Temperature upper bound = "<<high_T<<std::endl;
	std::cout<<"Temperature steps = "<<T_steps<<std::endl;
	std::cout<<"Measurements per step = "<<repetitions<<std::endl;

	if(repetitions > H.get_M()) throw(std::runtime_error("Too many measurements per step w.r.t. network size. Refer to manual for more information."));

	std::cout<<"Building Network patterns."<<std::endl;
	H.build_random_patterns();
	H.build_weights();

	std::vector<double> Ts (T_steps+1, 0.0);
	std::vector<std::vector<double>> m(T_steps+1, std::vector<double>(0,0.0));

	double delta_T = (high_T - low_T)/(T_steps);
for (int step = 0; step < T_steps+1; ++step)
{
	double cur_T = step*delta_T;
	Ts[step] = cur_T;
	std::cout<<"Setting T to "<<cur_T<<std::endl;
	H.set_temperature(cur_T);
	for(int rep = 0; rep < repetitions; ++rep)
	{
		std::cout<<"Init corrupted pattern "<<rep<<std::endl;
		H.init_on_corrupted_pattern(rep, 0.3);
		std::cout<<"Evolving system..."<<std::endl;
		H.glauber_evolve(n_iter);
		double ov = H.cur_overlap(rep);
		std::cout<<"Init pattern, Final overlap"<<std::endl;
		std::cout<<rep<<", "<<ov<<std::endl;
		m[step].push_back(ov);
	}

}

std::vector<double> averages;
std::vector<double> std_devs;

std::cout<<"T\t <m>\t sigma_m"<<std::endl;
std::cout<<"-------------------------"<<std::endl;

for (int i = 0; i < T_steps+1; ++i)
{
	averages.push_back(mean(m[i]));
	std_devs.push_back(std_dev(m[i]));
	std::cout<<std::setprecision(precision)<<std::setw(precision)<<Ts[i]<<"\t"<<averages[i]<<"\t"<<std_devs[i]<<std::endl;
	if (argc == 9)
		 {
		 	o<<std::setprecision(precision)<<std::setw(precision)<<Ts[i]<<"\t"<<averages[i]<<"\t"<<std_devs[i]<<std::endl;
		 }
}
	auto finish = std::chrono::steady_clock::now();

	if(argc == 9) std::cout<<"Output filepath specified, written results to: "<<outpath<<std::endl;
	else std::cout<<"No output filepath specified, written results to standard output only."<<std::endl;

	std::cout<<"Full run time: "<<std::chrono::duration<double>(finish - start).count()<<" s"<<std::endl;
}
