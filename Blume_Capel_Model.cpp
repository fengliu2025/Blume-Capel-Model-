#include <iostream>
#include <fstream>
#include <iomanip>  
#include <vector>
#include <random>
#include <string>


class Blume_Capel_Model{
public:
	//Simulation Parameters
	int N; 
	double J;
	double Delta;
	double H;
	double Beta;

	//Simulation Steps
	int Steps;

	double Energy;
	double m;

	double E_average;
	double E2_average;
	double m_average;
	//OutPut settings
	int Output_flag;
	std::string Output_file;

	
	std::mt19937 *gen;

	std::vector< std::vector<int> > Spins;

	Blume_Capel_Model(int n,double j, double delta, double h,double beta,int steps,int output_flag, std::string output_file );
	void Print_Info();
	int  Initialize();
	double Calculate_Energy();
	double Calculate_Magnetization();
	double Calculate_Q();
	double Evolve();

};


Blume_Capel_Model::Blume_Capel_Model(int n,double j, double delta, double h,double beta,int steps,int output_flag, std::string output_file )
		: N(n), J(j), Delta(delta),H(h), Beta(beta) ,Steps(steps), Output_flag(output_flag), Output_file(output_file)
{	
	std::random_device rd;
	gen = new std::mt19937(rd());

}


void Blume_Capel_Model::Print_Info(){
	std::cout<<"*********************************************************************************************"<<std::endl;
	std::cout<<"*   Blume-Capel Model: H = -J sum_{<i,j>} S_i S_j + Delta sum_{i} S_i^2 - H sum_{i} S_i      "<<std::endl;
	std::cout<<"*   Periodic Boundary Condition                                                              "<<std::endl;
	std::cout<<"*   Size: "<<N<<"*"<<N<<"																     "<<std::endl;	
	std::cout<<"*   J = "<< std::fixed <<std::setprecision(2)<<J<<"										     "<<std::endl;
	std::cout<<"*   H = "<< std::fixed <<std::setprecision(2)<<H<<"										     "<<std::endl;
	std::cout<<"*   Delta = "<< std::fixed <<std::setprecision(2)<<Delta<<"								     "<<std::endl;
	std::cout<<"*   Beta  = "<< std::fixed <<std::setprecision(2)<<Beta <<"								     "<<std::endl;
	std::cout<<"*********************************************************************************************"<<std::endl;

																					
}

int Blume_Capel_Model::Initialize(){
	std::uniform_int_distribution<int> states_rnd_gen(-1,1);
	//initialize the all spins
	for(int i=0; i < N;i++){
		std::vector<int> temp_spin;
		for(int j=0;j < N;j++){
			temp_spin.push_back( states_rnd_gen(*gen) );
		}
		Spins.push_back(temp_spin);
	}
	std::cout<<"Successfully Initialize the Spins." <<std::endl;
	Energy = Calculate_Energy();
	m = Calculate_Magnetization();

	E_average = 0;
	E2_average = 0;
	m_average = 0;

	return 0;

}

double Blume_Capel_Model::Calculate_Energy(){
	double Energy_J     = 0.;
	double Energy_Delta = 0.;
	double Energy_H     = 0.; 
	for(int i=0; i < N ;i++){
		for(int j=0; j < N;j++){
			Energy_J     = Energy_J     + 0.5 * (-J)  *  Spins[i][j] * ( Spins[ (i-1+N)%N ][j] + Spins[ (i+1+N)%N ][j] + Spins[i][ (j+1+N)%N ] + Spins[i][ (j-1+N)%N ]  );
			Energy_Delta = Energy_Delta + Delta *  Spins[i][j] * Spins[i][j];
			Energy_H     = Energy_H     + (-H)  *  Spins[i][j];
 		}
	}

	return (Energy_J+Energy_Delta+Energy_H);

}


double Blume_Capel_Model::Calculate_Magnetization(){
	double M = 0;
	for(int i=0; i < N ;i++){
		for(int j=0; j < N ;j++){
			M = M + Spins[i][j];
		}
	}

	return M/(N*N);
}

double Blume_Capel_Model::Calculate_Q(){
	double Q = 0.;
	for(int i=0;i < N ;i++){
		for(int j=0; j < N ; j++ ){
			Q = Q + Spins[i][j]*Spins[i][j];
		}
	}

	return Q/(N*N);
}


double Blume_Capel_Model::Evolve(){
	std::ofstream *output = 0;
	if(Output_flag){
		output=new std::ofstream(Output_file);
		*output <<"# 	N 	J 	Delta 	H 	Beta \n";
		*output <<"# " <<std::setw(5)<< N << std::setw(5) << J << std::setw(5) << Delta << std::setw(5) << H << std::setw(5) << Beta <<"\n";
		*output <<"# Steps \n";
		*output <<"# "<< Steps << "\n";
		*output <<"\n";
	}
	

	std::uniform_int_distribution<int> spins_rnd_gen(0,N-1);
	std::uniform_int_distribution<int> states_rnd_gen(0,1);
	std::uniform_real_distribution<double> p_rnd_gen(0.0,1.0);

	int K = 1;
	for(int i_step = 0; i_step < Steps; i_step++){
		
		//if(i_step%1000==0) std::cout<<"i_step="<<i_step<<std::endl;
		std::vector<int> i_trial;
		std::vector<int> j_trial;
		std::vector<int> ij_state_old  ;
		std::vector<int> ij_state_trial;
		for (int k = 0; k < K ; k++){
			i_trial.push_back( spins_rnd_gen(*gen) );
			j_trial.push_back( spins_rnd_gen(*gen) );
			int change = 0;
			if ( states_rnd_gen(*gen) == 0 ) change =-1;
			else change = +1;

			ij_state_old.push_back(  Spins[i_trial[k]][j_trial[k]] );
			ij_state_trial.push_back( (ij_state_old[k] + 4 + change )%3 -1 );
		}
		
		double delta_E = 0; 

		for(int k =0; k < K; k++ ){

		delta_E = delta_E +  (-J) * (ij_state_trial[k]-ij_state_old[k]) * (  Spins[ (i_trial[k]-1+N)%N  ][j_trial[k]] 
																		  + Spins[ (i_trial[k]+1+N)%N  ][j_trial[k]] 
																		  + Spins[ i_trial[k] ][ (j_trial[k]-1+N)%N  ] 
																		  + Spins[ i_trial[k] ][ (j_trial[k]+1+N)%N  ]  )  
								+ Delta * ( ij_state_trial[k] *ij_state_trial[k] - ij_state_old[k] * ij_state_old[k] )
								+ (-H) *  ( ij_state_trial[k] - ij_state_old[k] );
		}
		double E_trial = Energy + delta_E;



		if(E_trial <= Energy ){
			//Accept this trial step 
			Energy = E_trial;
			for(int k=0; k < K; k++){
				m = m - ((double)ij_state_old[k])/(N*N) + ((double)ij_state_trial[k]*1.0)/(N*N);
				Spins[i_trial[k]][j_trial[k]] = ij_state_trial[k];
			}
			//m = Calculate_Magnetization();
		}

		else {
			double p = p_rnd_gen(*gen);
			if( p < std::exp(-Beta * (E_trial-Energy) ) ){
				//Accept this trial step 
				Energy = E_trial;
				for(int k=0;k < K ;k++){
					m = m - ((double)ij_state_old[k])/(N*N) + ((double)ij_state_trial[k])/(N*N);
					Spins[i_trial[k]][j_trial[k]] = ij_state_trial[k];
				}
				//m = Calculate_Magnetization();
			}
			else {
				//Reject this trial step 

			}

		}
		//Now the system has evolved into a new state
		//Calculate some observables and/or output the Spins......
		//OutPut 
		if(Output_flag){
			*output <<"# "<<i_step<<"\n";
			for(int i=0; i < N ; i++){
				for(int j=0;j < N ; j++){
					*output<<Spins[i][j]<<" ";
				}
				*output<<"\n";
			}
		}
		
		m_average = m_average +std::abs(m);
		E_average = E_average + Energy ;
		E2_average = E2_average + Energy*Energy;
	

	}
	if(Output_flag) output->close();
	m_average = m_average/(Steps);
	E_average = E_average/(Steps);
	E2_average = E2_average/(Steps);
	//std::cout<<"<m>="<< m_average<<std::endl;

	return 0;







}

















