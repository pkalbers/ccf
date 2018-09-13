//
//  Cumulative Coalescent Function
//
//  Copyright © 2018 Patrick Albers. All rights reserved.
//

#include <cmath>
#include <array>
#include <fstream>
#include <limits>
#include <iostream>
#include <string>
#include <vector>

#define NSTATES 200 // hard-coded, change if needed


using states_t = std::array<double, NSTATES+1>;
using matrix_t = std::vector<states_t>;
using ct = std::vector<double>; // continuous
using dt = std::vector<size_t>; // discrete


ct infer_path(const ct obs)
{
	const size_t nobs = obs.size();
	
	
	// states
	
	states_t states;
	
	for (size_t i = 0; i <= NSTATES; ++i) 
		states[i] = double(i) / double(NSTATES);
		
	states[0] = 1E-16;
	states[NSTATES] = 1 - 1E-16;
	
	// matrices
	
	matrix_t T1(nobs);
	matrix_t T2(nobs);
	
	for (size_t i = 0; i <= NSTATES; ++i) 
	{
		T1[0][i] = std::log( (obs[0] * states[i]) + ((1.0 - obs[0]) * (1.0 - states[i])) );
		T2[0][i] = 0;
	}
	
	
	// forward
	
	for (size_t i = 1; i < nobs; ++i)
	{
		for (size_t j = 0; j <= NSTATES; ++j)
		{
			const double val = std::log( (obs[i] * states[j]) + ((1 - obs[i]) * (1 - states[j])) );
			double       max = -1 * std::numeric_limits<double>::infinity();
			size_t       arg = 0;
			
			for (size_t k = 0; k <= j; ++k)
			{
				const double tmp = T1[i-1][k];
				
				if (max < tmp)
				{
					max = tmp;
					arg = k;
				}
			}
			
			T1[i][j] = val + max;
			T2[i][j] = arg; 
		}
	}
	
	
	// traceback
	
	dt z(nobs);
	ct p(nobs);
	
	z[nobs-1] = NSTATES;
	p[nobs-1] = 1.0;
	
	for (size_t i = nobs-1; i > 0; --i)
	{
		z[i-1] = T2[ i ][ z[i] ];
		p[i-1] = states[ z[i-1] ];
	}
	
	
	return p;
}



int main(int argc, const char * argv[]) {
	
	const std::string name = argv[1];
	
	ct obs;
	
	std::ifstream file(name);
	std::string line;
	
	while (file && getline(file, line))
	{
		if (line.size() == 0)
			continue;
		
		if (line[0] == '0') 
			obs.push_back(0.0);
		else if (line[0] == '1') 
			obs.push_back(1.0);
		else return(1);
	}
	
	
	const ct path = infer_path(obs);
	
	for (const double p : path)
		std::cout << p << std::endl;
	
	return 0;
}
