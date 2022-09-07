/**
 * @file HopfieldNetwork.h
 *
 *\author Giulio Colombini
 *
 * Hopfield Network with built-in parallel Monte Carlo Glauber evolution.
 *
 * This header contains the full implementation of a Hopfield Network,
 * along with some useful functions to work with it.
 */

#include<cmath>
#include<random>
#include<vector>
#include<algorithm>
#include<utility>
#include<stdexcept>

#ifndef N_PARALLEL_THREADS
//! The number of parallel threads used to run the program defaults to 1, but can be set at compile time depending on the
//! possibilities of the available machine by passing ```-DN_PARALLEL_THREADS=N```, where ```N``` is the number of desired threads.
	#define N_PARALLEL_THREADS 1
#endif

//! Signed char used to represent a single Ising spin which can only take the values +1 or -1 to optimize memory usage.
typedef signed char spin;

//! Shorthand to represent a pattern of spins to store in the Network.
typedef std::vector<spin> spin_pattern;


//! Common random engine for usage throughout the program.
static 	std::default_random_engine re;
//! Binary Random Variable.
static  std::uniform_int_distribution<short int> coin_toss(0,1);
//! Uniform real Random Variable in the interval [0,1].
static  std::uniform_real_distribution<double> rnd(0.0,1.0);

inline spin random_spin() 
//!
//! Returns a ```spin``` which has value +1 or -1 with equal probability 1/2.
//!
{return 2 * coin_toss(re) - 1;}

double overlap (spin * a, spin * b, int N)
//!
//! Returns the overlap between spin configuration \f$\sigma^a\f$ and \f$\sigma^b\f$, both of which must be of size N,
//! amounting to \f[ \mathrm{Overlap}(\sigma^a,\,\sigma^b) = \frac{1}{N}\sum_{i = 1}^{N} \sigma^a_i \sigma^b_i \textrm{.}\f].
//! When calculated between the current spin configuration and one of the memories, this quantity is often called
//! the *Mattis magnetisation* and indicated, e.g. for memory \f$\mu\f$, by \f$m^{\mu}\f$. 
//!
{
	double ret = 0;
	for(int i = 0; i < N; ++i) ret += a[i]*b[i];

	return ret/N;
}

class HopfieldNetwork
/**
 * @brief      Implementation of a Hopfield Network of binary neurons (spins).
 *
 * @param[in]  <unnamed>  { parameter_description }
 *
 * @return     { description_of_the_return_value }
 */
{
	private:
	double T = 0.;									// System temperature
	double alpha;	 								// System load parameter
	int N;											// Number of spins/neurons
	int M;											// Number of memories/patterns
	spin *  spins;									// Pointer to the state of the Network Spins
	std::vector<spin_pattern> patterns;				// Vector containing the patterns stored in the Network
	double ** W;									// Interaction weights between neurons
	std::uniform_int_distribution<int> spin_picker; // Distribution used for selecting spins to update
	bool initialised_weights = false;				// Keep track of weight initialisations, to avoid running simulations without uninitialised parameters in the model.

	public:

	HopfieldNetwork(int _N, int _M, double _T): T{_T}, N{_N}, M{_M}
	//!
	//! Construct a Hopfield Network with ```N``` spins, able to store ```M``` patterns of ```N``` spins, with initial temperature ```T```.
	//!
	{
		alpha = double(M/N);
		spins = new spin[N];
		patterns = std::vector<spin_pattern>(M, spin_pattern());
		W = new double * [N];
		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
		for(int i = 0; i < N; ++i) W[i] = new double[N];
		spin_picker = std::uniform_int_distribution<int>(0, N-1);
	}

	HopfieldNetwork(int _N, double _alpha, double _T): T{_T}, alpha{_alpha}, N{_N}
	//!
	//! Construct a Hopfield Network with ```N``` spins, able to store \f$ \alpha\f$ ```N``` patterns of ```N``` spins, with initial temperature ```T```.
	//! \f$ \alpha \f$ is generally known in literature as the load parameter.
	//! 
	{
		M = int(std::round(alpha * N));
		spins = new spin[N];
		patterns = std::vector<spin_pattern>(M, spin_pattern());
		W = new double * [N]; 
		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
		for(int i = 0; i < N; ++i) W[i] = new double[N];
		spin_picker = std::uniform_int_distribution<int>(0, N-1);
	}

	void init_spins_randomly()
	//!
	//! Initialise all the spins in a random configuration, using random_spin().
	//! 
	{for(int i = 0; i < N; ++i) spins[i] = random_spin();}

	void init_on_corrupted_pattern(spin_pattern pattern, double probability)
	//! 
	//! Initialise the spins with a corrupted version of the pattern pointed by pattern, where each of the spins
	//! may have been flipped with probability probability.
	//! 
	{
		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
		for(int i = 0; i < N; ++i){spins[i] = rnd(re) < probability? -1 * pattern[i] : pattern[i];}}

	void init_on_corrupted_pattern(int i, double probability)
	//! 
	//! Initialise the spins with a corrupted version of the ```i```-th pattern of the ```M``` stored in the network, where each of the spins
	//! may have been flipped with probability ```probability```. Note that pattern must point to ```N``` spins.
	//! 
	{
		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
		for(int j = 0; j < N; ++j){spins[j] = rnd(re) < probability? -1 * patterns[i][j] : patterns[i][j];}}


	// Parametric setters

	void set_temperature(double newT) 
	//!
	//! Set the network operation temperature to newT.
	//! 
	{ T = newT; }

	// Memory setters

	void build_random_patterns(int n_patterns = -1)
	//!
	//! Build and store a number n_patterns of randomly built patterns in the network. If nothing or -1 is passed, build and store M patterns.
	//! This method is useful as a benchmark since the most simple theoretical results have been proved for random patterns.
	{	
	if (n_patterns == -1) n_patterns = M; // Defaults to the number specified at construction

	#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
	
	for(int p = 0; p < M; ++p) 
	{
		patterns[p].resize(N);
		for(int n = 0; n < N; ++n) patterns[p][n] = random_spin();
	}
	initialised_weights = false;
	}

	void push_back_pattern(spin_pattern p)
	/**
	 * Store another pattern in the Network. This method falsifies the internal switch ```initialised_weights``` so that launching a simulation
	 * before calling the function build_weights() throws. This method throws if one attempts to push_back a pattern of size() different from ```N```.
	 */
	{
		if (p.size() == (unsigned)N) {this->patterns.push_back(p);}
		else {throw std::runtime_error("The length of each pattern must be equal to the number of spins in the system!");}
		this->M = this->patterns.size();
		this->alpha = double(this->M / this->N);
		initialised_weights = false;
	}

	void set_M(int newM)
	{
		int deltaM = newM - this->M;
		if(deltaM > 0)
		{
			spin_pattern tmp(N, 0);
			for(int i = 0; i < deltaM; ++i)
			{
				for(int j = 0; j < N; ++j) tmp[j] = random_spin();
				this->patterns.push_back(tmp);			
			}
			this->M = newM;
			this->alpha = double(this->M/this->N);
			initialised_weights = false;
		}
		if (deltaM == 0) return;
		if (deltaM < 0)
		{
			for(int i = 0; i < std::abs(deltaM); ++i) this->patterns.pop_back();
			this->M = newM;
			this->alpha = double(this->M/this->N);
			initialised_weights = false;
		}

	}

	void set_alpha(double newalpha)
	{
		int newM = int(std::round(newalpha));
		set_M(newM);
	}

	// Weights initialiser

	void build_weights()
	/**
	 * Construct the matrix of weights from the patterns stored in ```this->patterns``` and store it in ```this->W```. This is the matrix generally
	 * indicated in the literature with \f$ J_{ij} \f$.
	 */
	{
	#pragma omp parallel for num_threads(N_PARALLEL_THREADS)

	for(int i = 0; i < N; ++i)
	{	
		for(int j = 0; j <= i; ++j)
		{
			W[i][j] = 0;
			for(int w = 0; w < M; ++w) W[i][j] += (double)patterns[w][i] * patterns[w][j] / N;
			W[j][i] = W[i][j];
		}
	}
	initialised_weights = true;
	}

	// Evolution step function
	void glauber_evolve(unsigned int niter, unsigned int nflips = N_PARALLEL_THREADS)
	/**
	 * Evolve the network via 
	 */
	{

	for(auto it : patterns) if(it.size() != (unsigned)N) throw std::runtime_error("Simulation was launched with uninitialised patterns.");

	if (not initialised_weights) throw std::runtime_error("Simulation was launched with uninitialised weights.");

	for(unsigned int iter = 0; iter < niter; ++iter)
	{
		int flip_candidates[nflips];
		for(unsigned int i = 0; i < nflips; ++i) flip_candidates[i] = spin_picker(re);

		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)

		for(unsigned int i = 0; i < nflips; ++i)
		{	
			double lf = 0;
			for(int j = 0; j < N; ++j) lf+= W[flip_candidates[i]][j]*spins[j];

			double dE  = 2 * lf * spins[flip_candidates[i]];
			double thr = 1/(1+std::exp(dE/T));

			if(rnd(re) < thr) 
				{
					spins[flip_candidates[i]] *= -1;
				}
		}
	}
	}

	std::vector<double> overlaps()
	//!
	//! Return a std::vector of M elements, containing the memory overlaps of the current internal state of the network. The \f$ \mu \f$-th element
	//! being the overlap with pattern \f$ \mu \f$, a.k.a. the Mattis magnetisation \f$ m^{\mu} \f$.
	//! 
	{
		std::vector<double> overlaps(N, 0.);
		for(int o = 0; o < M; ++o) overlaps[o] = overlap(spins, patterns[o].data(), N);
		return overlaps;
	}

	std::pair<int, double> max_overlap()
	{
		auto overlaps = this->overlaps();
		int argmax = std::distance(overlaps.begin(), std::max_element(overlaps.begin(), overlaps.end()));
		return std::pair<int, double>(argmax, overlaps[argmax]);
	}
};
