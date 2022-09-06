#include<cmath>
#include<cstdlib>
#include<random>

#define N_PARALLEL_THREADS 12

typedef signed char spin;

// Common random engine for usage throughout the program.
static 	std::default_random_engine re;
static  std::uniform_int_distribution<short int> coin_toss(0,1);
static  std::uniform_real_distribution<double> rnd(0.0,1.0);

inline spin random_spin() {return 2 * coin_toss(re) - 1;}

double overlap (spin * s, spin * t, int N)
{
	double ret = 0;
	for(int i = 0; i < N; ++i) ret += s[i]*t[i];

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
	spin ** patterns;								// Pointer to the patterns stored in the Network
	double ** W;									// Interaction weights between neurons
	std::uniform_int_distribution<int> spin_picker; // Distribution used for selecting spins to update

	public:

	// Constructors
	HopfieldNetwork(int _N, int _M, double _T): T{_T}, N{_N}, M{_M}
	{
		alpha = double(M/N);
		spins = new spin[N];
		patterns = new spin *[M];
		for(int i = 0; i < M; ++i) patterns[i] = new spin[N];
		W = new double * [N];
		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
		for(int i = 0; i < N; ++i) W[i] = new double[N];
		spin_picker = std::uniform_int_distribution<int>(0, N-1);
	}

	HopfieldNetwork(int _N, double _alpha, double _T): T{_T}, alpha{_alpha}, N{_N}
	{
		M = int(std::round(alpha * N));
		spins = new spin[N];
		patterns = new spin *[M];
		for(int i = 0; i < M; ++i) patterns[i] = new spin[N];
		W = new double * [N]; 
		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
		for(int i = 0; i < N; ++i) W[i] = new double[N];
		spin_picker = std::uniform_int_distribution<int>(0, N-1);
	}
	
	// Spin configuration initialisers

	void init_spins() 
	{ 
		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
		for(int i = 0; i < N; ++i) spins[i] = random_spin();
	}

	void init_on_corrupted_pattern(spin * pattern, double percentage)
	{
		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
		for(int i = 0; i < N; ++i){spins[i] = rnd(re) < percentage? -1 * pattern[i] : pattern[i];}}

	void init_on_corrupted_pattern(int i, double percentage)
	{
		#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
		for(int j = 0; j < N; ++j){spins[j] = rnd(re) < percentage? -1 * patterns[i][j] : patterns[i][j];}}

	// Setters

	void set_temperature(double newT) { T = newT; }

	void build_random_patterns(int n_patterns = -1)
	{	
	if (n_patterns == -1) n_patterns = M; // Defaults to the number specified at construction

	#pragma omp parallel for num_threads(N_PARALLEL_THREADS)
	
	for(int p = 0; p < M; ++p)
	{
		for(int b = 0; b < N; ++b) patterns[p][b] = random_spin();
	}
	
	}

	void build_weights()
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
	}

	// Evolution step function
	void glauber_evolve(unsigned int niter, unsigned int nflips = N_PARALLEL_THREADS)
	{
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

	int overlap_argmax()
	{
	double * overlaps = new double[M];
	
	for(int o = 0; o < M; ++o) overlaps[o] = overlap(spins, patterns[o], N);

	int argmax = 0;
	double max = 0;

	for(int idx = 0; idx < M; ++idx) 
	{
		if (overlaps[idx] > max) 
			{
				max = overlaps[idx];
			    argmax = idx;
			}

	}
	return argmax;
	}

	double overlap_max()
	{
	double * overlaps = new double[M];
	
	for(int o = 0; o < M; ++o) overlaps[o] = overlap(spins, patterns[o], N);

	double max = 0;

	for(int idx = 0; idx < M; ++idx) 
	{
		if (overlaps[idx] > max) 
			{
				max = overlaps[idx];
			}

	}
	return max;
	}
};
