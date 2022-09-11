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

/*! \mainpage Hopfield Networks in C++
 *
 * \section intro Introduction
 * This header-only library implements a Hopfield Network in C++, aiming for a lightweight implementation 
 * that makes use of parallelism whenever possible. The time evolution is performed according to a parallel version
 * of the conventional Glauber algorithm. Albeit being quite simple in its definition, the Hopfield Model can be regarded
 * as a minimal version of a recurrent neural network implementing an associative, content addressed memory.
 * Moreover, the study of the properties of the model in the presence of noise can be carried out in the framework of
 * Statistical Physics, thus providing a paradigmatic example of the physical properties of a disordered system.
 * 
 * The implementation of parallel procedures is achieved via the OpenMP #pragma directives, so to exploit a secure and
 * complete interface to the parallel features of C++. The default number of parallel threads can be set at compile time,
 * so to meet the specifics of the available architecture, while for some methods also runtime specification is made possible,
 * as detailed in the following documentation.
 * 
 * \section theory Theoretical Background
 *
 * The Hopfield Model or Hopfield Network is a very simple model for an associative memory. 
 * First put forward by William Little in 1974 and then developed by John Hopfield, it was devised to explain in a simplified
 * context, the associative nature of memory in the brain, e.g. the fact that the recognition of an object can be triggered also
 * by a partial or modified version of the memory itself. The model is very simplistic from the point of view of biological realism,
 * the neurons being represented by McCulloch-Pitts binary units, which can only be in one of two states: \f$+1\f$, active or \f$-1\f$, inactive.
 * Considering discretised units of time, each corresponding ideally to the average refractory period of a biological neuron, and representing
 * the internal state of neuron \f$i\f$ at time instant t as \f$\sigma_i^t\f$, the McCulloch-Pitts update rule for each neuron in a network
 * containing \f$N\f$ is defined as
 * 
 * \f[\sigma_i^{t+1} = \textrm{sgn}\left(\sum_{j=1}^{N} W_{ij} \sigma_j^t\right)\f]
 * 
 * where the \f$W_{ij}\f$ are called the synaptic weights, and \f$\textrm{sgn}\f$ is the sign function.
 * The choice of which, and how many, spins to update per time unit is part of the implementation, as in general the details tend not to affect
 * the final equilibrium state.
 * 
 * To have the Network work as an associative memory we need to store into it a number \f$p\f$ of memories, in the form of spin patterns. The idea is
 * that an associative memory, if put in a configuration close to one of the stored patterns (e.g. a corrupted version of it), should reconstruct it 
 * during its time evolution, i.e. relax to a configuration \f$\vec{\sigma}\f$ of the spins equal to the closest stored pattern.
 * 
 * To do so, with the previously defined dynamics, it can be proved that it is sufficient to define the weights \f$W_{ij}\f$ as follows. Letting the
 * patterns \f$\xi^\mu_i\f$ be indexed by a greek index such as \f$ 1 \leq \mu \leq p\f$ and the spins within each pattern with a regular latin index 
 * such as \f$1 \leq i \leq N\f$, we define the \f$W_{ij}\f$ as
 * 
 * \f[W_{ij} = \frac{1}{N}\sum_{\mu = 1}^{p} \xi^\mu_i \xi^\mu_j\f]
 * 
 * This choice of weights is called *Hebbian rule*, from the connectionist psychologist Donald Hebb, and it is generally summarised by the phrase
 * *fire together, bind together*, meaning that the synapses connecting neurons which activate together are reinforced, while those connecting neurons
 * that seldom fire together are weakened.
 * 
 * It can be proved \cite coolen2005theory that under very simple hypotheses (symmetric weights), that the function
 * 
 * \f[\mathcal{H}(\mathbf{\sigma}) = -\frac{1}{2} \sum_{i = 1}^{N} W_{ij} \sigma_i\sigma_j\f]
, * 
 * is a Lyapunov function (i.e. non-decreasing along the system's orbits) for the deterministic dynamics. From a physicist's viewpoint,
 * we are saying that the deterministic dynamics tends to minimise a Ising-like spin Hamiltonian.
 * 
 * It can be proved similarly \cite coolen2005theory, that the deterministic dynamics can lead the network into spurious minima, i.e. 
 * minima which correspond to the combination of a finite number of memories, *spurious mixtures*, or even an extensive number of memories, *glassy states* (Of course to access
 * this kind of states it is necessary to have an extensive number of available memories, i.e. we need to have the scaling \f$ p = \alpha N\f$).
 * 
 * To avoid these states, or to render them unstable (repulsive) equilibria, we can introduce a stochastic dynamics into the system. A sensible
 * choice would be to introduce a temperature and select a dynamics which is compatible with the Boltzmann-Gibbs equilibrium distribution 
 * induced by the spin Hamiltonian, this way we have a natural parametrisation for the noise level and can employ the machinery of Statistical
 * Physics to characterise the macroscopic states of the system.
 * To do so, we can select as update rule any method from the Monte Carlo techniques applied to spin systems.
 * We choose *Glauber dynamics*, which amount to selecting a spin at random and flipping it with probability
 * 
 * \f[\mathbb{P}(\sigma_i \to - \sigma_i) = \frac{1}{1+e^{\beta \Delta \mathcal{H}_i}} \f]
 * 
 * where by \f$ \Delta\mathcal{H}_i = 2 \sum_j W_{ij} \sigma_j \sigma_i \f$ we denote the energy change caused by a flip of spin \f$i\f$.
 * 
 * This dynamics is compatible with the Boltzmann-Gibbs equilibrium distribution and so we can set out to determine the phase diagram of the model
 * using the tools of Statistical Physics. In particular we are interested in the case in which the number of patterns is extensive with \f$N\f$, so
 * \f$p = \alpha N\f$, and we want to characterise the working features of the network as an associative memory in function of the noise level
 * (temperature) \f$T\f$ and the load parameter \f$\alpha\f$. The Statistical Mechanical treatment of the problem is very interesting and 
 * makes use of Replica Methods to deal with the *quenched disorder* brought on by the distribution of the memories. The resulting phase diagram is
 * presented in the figure below (from \cite coolen2005theory).
 * 
 * \image html phase_diagram.png
 * \image latex phase_diagram.png
 * 
 * The different phases correspond to the case in which the recall states are absolute minima of the system free energy (F), local minima (M)
 * or neither, due to the glassy nature of the free energy landscape (SG). The P phase corresponds to a paramagnetic fully disordered phase.
 * The regions in which the network is considered to be working as an associative memory are F and M, since provided that the initial state
 * is not too corrupted, the network will reconstruct the original pattern with its dynamics. To measure this property the pattern overlaps (also known
 * as Mattis magnetisations) are introduced, the overlap between spin configuration \f$\sigma\f$ and pattern \f$\xi^\mu\f$ being defined as
 * 
 * \f[ m^\mu(\sigma) = \frac{1}{N}\sum_{i = 1}^{N} \sigma_i \xi^\mu_i \f]
 * 
 * where if the \f$\sigma\f$ is omitted it is implied that the overlap is to be taken with the current spin state of the network.
 * 
 * Recall states, also known as pure states, correspond to a value of the appropriate \f$m^\mu = 1 \f$, and with a pure state ansatz one can 
 * obtain a self consistent equation for the amplitude \f$m^\mu\f$ as a function of \f$T,\, \alpha\f$, yielding the following plot (from \cite coolen2005theory)
 *
 * \image html magnetisation_cropped.png
 * \image latex magnetisation_cropped.png
 * 
 * Where the various lines correspond, from top to bottom, to values of \f$\alpha = 0,\, \mathellipsis,\, 0.125\f$ in increments of \f$0.025\f$, and
 * dashed lines indicate the vanishing temperature for the amplitude at the given \f$\alpha\f$.
 */


#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

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

	~HopfieldNetwork()
	/**
	 * Destroy a Hopfield Network instance. A handwritten destructor is needed for the raw pointers.
	 */
	{
		alpha = 0;
		delete[] spins;
		for(int i = 0; i < N; ++i) delete[] W[i];
		delete[] W;
	}
	void init_spins_randomly()
	//!
	//! Initialise all the spins in a random configuration, using random_spin().
	//! 
	{for(int i = 0; i < N; ++i) spins[i] = random_spin();}

	void init_on_corrupted_pattern(spin_pattern pattern, double probability)
	//! 
	//! Initialise the spins with a corrupted version of the pattern pointed by pattern, where each of the spins
	//! may have been flipped with probability ```probability```.
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
	/**
	 * Set a new value for the number of stored patterns. Depending on the value of newM, one of the following three cases can apply.
	 * 
	 * If newM = M, nothing is changed, if newM < M, (M - newM) patterns are popped back from
	 * the tail of the patterns vector. If newM > M, (newM-M) new random patterns are pushed back to the patterns vector.
	 * 
	 * Calling this method falsifies the internal switch ```initialised_weights``` so that launching a simulation without calling the function build_weights()
	 * throws an exception.
	 */
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
	/**
	 * Set a new value for the load parameter \f$ \alpha \f$. Defining newM = round(newalpha * N), one of the following three cases can apply.
	 * 
	 * If newM = M, nothing is changed, if newM < M, (M - newM) patterns are popped back from
	 * the tail of the patterns vector. If newM > M, (newM-M) new random patterns are pushed back to the patterns vector.
	 * 
	 * Calling this method falsifies the internal switch ```initialised_weights``` so that launching a simulation without calling the function build_weights()
	 * throws an exception.
	 */
	{
		int newM = int(std::round(newalpha * this->N));
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
	 * Evolve the network using a parallel version of the Glauber algorithm.
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

	double get_Energy(spin * state)
	/**
	 * Return the Energy of the Network evaluated for the spin configuration pointed by ```state```. ```state``` must point to N instances of ```spin```. 
	 */
	{
		double ret = 0;

		#pragma omp parallel for reduction (-:ret) num_threads(N_PARALLEL_THREADS)
		for(int i = 0; i < N; ++i)
		{
			for(int j = i+1; j < N; ++j) ret -= W[i][j]*state[i]*state[j];
		}

		return ret;
	}

	double get_Energy()
	/**
	 * Return the energy of the Network for the current internal state of ```this->spins```.
	 */
	{
		return get_Energy(this->spins);
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
	/**
	 * Return a std::pair<int, double> containing the index of the most condensed pattern and the corresponding Mattis magnetisation.
	 * They can be easily accessed through the .first and .second members of the std::pair.
	 */
	{
		auto overlaps = this->overlaps();
		int argmax = std::distance(overlaps.begin(), std::max_element(overlaps.begin(), overlaps.end()));
		return std::pair<int, double>(argmax, overlaps[argmax]);
	}
};
