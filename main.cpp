#include "Blume_Capel_Model.cpp"
#include <sstream>
#include <fstream>
int main(int argc, char* argv[]){
	if(argc !=2){
		std::cout<<"Please run by: ./BC_Model [i_d]" <<std::endl;
		return 1;
	}
	int i_d = std::stoi(argv[1]);

	int N = 20;
	double J = 1;
	double H = 0.0 ;

	//The range of Delta and T to be scanned 
	double Delta_min = 0 ;
	double Delta_max = 2.2;
	double T_min = 0.1;
	double T_max = 2  ;

	int N_Delta = 220;
	int N_T     = 200; 

	if(i_d > N_Delta){
		std::cout<<"id > N_Delta" <<std::endl;
		return 1;
	}



	double Delta_step = (Delta_max - Delta_min)/(N_Delta);
	double T_step     = (T_max - T_min)/(N_T);

	std::vector< double >  Delta; 
	std::vector< double >  Beta; 
	Delta.push_back(Delta_min);
	Beta.push_back(1./T_min);

	for(int i = 0;i<N_Delta;i++){
		Delta.push_back( Delta_min + (i+1)*Delta_step  );
	}

	for(int i = 0;i<N_T;i++){
		Beta.push_back(  1./( T_min+(i+1)*T_step )   );
	}


	int Steps = 100000000; 
	int Output_flag = 0; 
	std::stringstream s;
	s<<"phase_space_scan_i_d_"<<i_d<<".txt";
	std::ofstream phase_space_scan(s.str());
	phase_space_scan<<"#"<<std::setw(5)<<"N"<<std::setw(5)<<"J"<<std::setw(5)<<"D"<<std::setw(5)<<"H"<<std::setw(5)<<"B\n";
	phase_space_scan<<"#"<<std::setw(10)<<"<|m|>"<<std::setw(10)<<"<E>"<<std::setw(10)<<"<E^2>\n";
	phase_space_scan<<"\n";

	
	for(int i_b = 0; i_b < Beta.size();i_b++){
		std::stringstream ss;
		ss <<"N_" << N << 
		 	 "J_" << std::fixed << std::setprecision(2) << J << 
		 	 "D_" << std::fixed << std::setprecision(2) << Delta[i_d] << 
		 	 "H_" << std::fixed << std::setprecision(2) << H <<
		 	 "B_" << std::fixed << std::setprecision(2) << Beta[i_b] ; 
		std::string Output_file=ss.str();
		Blume_Capel_Model BC_Model=Blume_Capel_Model(N,J,Delta[i_d],H,Beta[i_b],Steps,Output_flag,Output_file);
		BC_Model.Print_Info();
		BC_Model.Initialize();
		BC_Model.Evolve();
			
		phase_space_scan<<"#"<<std::setw(10) <<N << std::setw(10) << J << std::setw(10) << Delta[i_d] << std::setw(10) << H << std::setw(10) << Beta[i_b] <<"\n";
		phase_space_scan<<BC_Model.m_average << std::setw(20) << BC_Model.E_average << std::setw(20) << BC_Model.E2_average << "\n";
		phase_space_scan<<"\n";
	}
	
	phase_space_scan.close();
	return 0;	
}
