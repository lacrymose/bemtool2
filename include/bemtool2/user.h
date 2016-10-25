#ifndef USER_H
#define USER_H

#include <string>
#include <sstream>
#include <vector>


////========================================================////
/////////////////===== Conversions ======///////////////////////

template <typename nbr>
std::string NbrToStr(nbr N){
	std::ostringstream strs;
	strs << N;
	std::string str = strs.str();
	return str;
}

template <typename T>
T StrToNbr ( const std::string &Text )
{
	std::istringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

inline int StrToInt(std::string str){
	std::stringstream i(str);
	int  N;
	i >> N;
	return N;
}

inline Real StrToReal(std::string str){
	std::stringstream i(str);
	Real  N;
	i >> N;
	return N;
}

inline Cplx StrToCplx(std::string str){
	std::stringstream i(str);
	Cplx  N;
	i >> N;
	return N;
}

////========================================================////
///////////////////////===== Input ======///////////////////////
inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
inline std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

// //// Charge les inputs : int N, Real lc, string 
// void GetInput(string data_name, int& N,Real& lc, bool& type,vector<Real>& kappa){
// // 	string data_name="../input/input_"+NbrToStr(N)+"_"+NbrToStr(lc)+".txt";
// 	ifstream data(data_name.c_str());
// 	
// 	// Si le fichier n'existe pas
// 	if (!data){
// 		cerr << "Input file doesn't exist" << endl;
// 		exit(1);
// 	}
// 	// Lecture du fichier
// 	else {
// 		while (data){
// 			string strInput;
// 			getline(data,strInput);
// // 			cout<<strInput<<endl;
// 			vector<string> line = split (strInput,' ');
// 			if (!line.empty()){
// 				if (line.at(0)=="Nombre"){
// 					N=StrToInt(line.back());
// 				}
// 				else if (line.at(0)=="Finesse"){
// 					lc=StrToReal(line.back());
// 				}
// 				else if(line.at(0)=="Propagatif"){
// 					type=StrToInt(line.back());
// 				}
// 				else if(line.at(0)=="Kappa"){
// 					assert(line.size()==N+3);
// 					for (int i=0;i<N+1;i++){
// 						kappa.push_back(StrToReal(line.at(2+i)));
// 					}
// 				}
// 			}
// 		}
// 	}
// }
// 
// void GetInputDDM(string data_name, int& type_mesh, Real& lc ,bool& type,vector<Real>& kappa,string& DDM, int& S, int& extension){
// // 	string data_name="../input/input_"+NbrToStr(N)+"_"+NbrToStr(lc)+".txt";
// 	ifstream data(data_name.c_str());
// 	int N=0;
// 	// Si le fichier n'existe pas
// 	if (!data){
// 		cerr << "Input file doesn't exist" << endl;
// 		exit(1);
// 	}
// 	// Lecture du fichier
// 	else {
// 		while (data){
// 			string strInput;
// 			getline(data,strInput);
// 			vector<string> line = split (strInput,' ');
// 			if (!line.empty()){
// 				if(line.at(0)=="Maillage"){
// 					type_mesh=StrToInt(line.back());
// 				}
// 				else if(line.at(0)=="Finesse"){
// 					lc=StrToReal(line.back());
// 				}
// 				else if(line.at(0)=="Propagatif"){
// 					type=StrToInt(line.back());
// 				}
// 				else if(line.at(0)=="Kappa"){
// 					assert(line.size()==N+3);
// 					for (int i=0;i<N+1;i++){
// 						kappa.push_back(StrToReal(line.at(2+i)));
// 					}
// 				}
// 				else if(line.at(0)=="DDM"){
// 					DDM=line.back();
// 					if (DDM!="none" & DDM!="asm" & DDM!="ras" & DDM!="ilu" & DDM!="diag"){
// 						cout<<DDM<<" is not an acceptable argument for DDM"<<endl;
// 						exit(1);
// 					}
// 				}
// 				else if(line.at(0)=="Sous-domaines"){
// 					S=StrToInt(line.back());
// 				}
// 				else if(line.at(0)=="Extension"){
// 					extension=StrToInt(line.back());
// 				}
// 				
// 			}
// 		}
// 	}
// }
// 
// void GetInput2(string data_name, int& type_mesh, Real& lc ,bool& type,vector<Real>& kappa){
// // 	string data_name="../input/input_"+NbrToStr(N)+"_"+NbrToStr(lc)+".txt";
// 	ifstream data(data_name.c_str());
// 	int N=0;
// 	// Si le fichier n'existe pas
// 	if (!data){
// 		cerr << "Input file doesn't exist" << endl;
// 		exit(1);
// 	}
// 	// Lecture du fichier
// 	else {
// 		while (data){
// 			string strInput;
// 			getline(data,strInput);
// 			vector<string> line = split (strInput,' ');
// 			if (!line.empty()){
// 				if(line.at(0)=="Maillage"){
// 					type_mesh=StrToInt(line.back());
// 				}
// 				else if(line.at(0)=="Finesse"){
// 					lc=StrToReal(line.back());
// 				}
// 				else if(line.at(0)=="Propagatif"){
// 					type=StrToInt(line.back());
// 				}
// 				else if(line.at(0)=="Kappa"){
// 					assert(line.size()==N+3);
// 					for (int i=0;i<N+1;i++){
// 						kappa.push_back(StrToReal(line.at(2+i)));
// 					}
// 				}
// 			}
// 		}
// 	}
// }
// 
// void GetInput(string data_name, int& N,Real& lc, bool& type,vector<Real>& kappa, Real& alpha){
// // 	string data_name="../input/input_"+NbrToStr(N)+"_"+NbrToStr(lc)+".txt";
// 	ifstream data(data_name.c_str());
// 	
// 	// Si le fichier n'existe pas
// 	if (!data){
// 		cerr << "Input file doesn't exist" << endl;
// 		exit(1);
// 	}
// 	// Lecture du fichier
// 	else {
// 		while (data){
// 			string strInput;
// 			getline(data,strInput);
// // 			cout<<strInput<<endl;
// 			vector<string> line = split (strInput,' ');
// 			if (!line.empty()){
// 				if (line.at(0)=="Nombre"){
// 					N=StrToInt(line.back());
// 				}
// 				else if (line.at(0)=="Finesse"){
// 					lc=StrToReal(line.back());
// 				}
// 				else if(line.at(0)=="Propagatif"){
// 					type=StrToInt(line.back());
// 				}
// 				else if(line.at(0)=="Kappa"){
// 					assert(line.size()==N+3);
// 					for (int i=0;i<N+1;i++){
// 						kappa.push_back(StrToReal(line.at(2+i)));
// 					}
// 				}
// 				else if(line.at(0)=="Relaxation"){
// 					alpha=StrToReal(line.back());
// 				}
// 			}
// 		}
// 	}
// }
// 
// 
// //// Charge les inputs : int N, Real lc, string 
// void GetInput_MPI(string data_name, int N, Real& lc, bool& type,vector<Real>& kappa){
// 	vector<string> test=split(data_name,'_');
// 	if (N!=StrToInt(test[1])-1){
// 		cerr << "N expected= "<<N<<endl;
// 		cerr << "N given = "<<StrToInt(test[1])-1<<endl;
// 		cerr << "Wrong input or number of procs" << endl;
// 	}
// 	
// 	ifstream data(data_name.c_str());
// 	// Si le fichier n'existe pas
// 	if (!data){
// 		cerr << "Input file doesn't exist" << endl;
// 	}
// 	// sinon : Lecture du fichier
// 	else {
// 		while (data){
// 			string strInput;
// 			getline(data,strInput);
// 			vector<string> line = split (strInput,' ');
// 			if (!line.empty()){
// 				if (line.at(0)=="Finesse"){
// 					lc=StrToReal(line.back());
// 				}
// 				else if(line.at(0)=="Propagatif"){
// 					type=StrToInt(line.back());
// 				}
// 				else if(line.at(0)=="Kappa"){
// 					if (line.size()!=N+3){
// 						cerr<<"Lenght of line : "<<line.size()<<endl;
// 						cerr<<"Lenght expected : "<<N+3<<endl;
// 						cerr<<"Error"<<endl;
// 					}
// 					for (int i=0;i<N+1;i++){
// 						kappa.push_back(StrToReal(line.at(2+i)));
// 					}
// 				}
// 				
// 			}
// 		}
// 	}
// }
// 
// ////========================================================////
// ///////////////////////===== Output ======//////////////////////
// 
// void Create_output_1(int N, Real lc, int iteration, bool type, vector<double> times){
// 	string filename="../output/output_1_"+NbrToStr(N+1)+"_"+NbrToStr(lc)+"_"+NbrToStr(type)+".txt";
// 	ofstream output(filename.c_str());
// 	if (!output){
// 		cerr<<"Output file cannot be created"<<endl;
// 		exit(1);
// 	}
// 	else{
// 		output<<N<<" "<<iteration<<endl;
// 		for (int j=0;j<times.size()-1;j++){
// 			output<<times.at(j)<<" ";
// 		}
// 		output<<times.back()<<endl;;
// 	}
// }
// 
// void Create_output_2(int N, Real lc, int iteration, bool type, vector<double> times){
// 	string filename="../output/output_2_"+NbrToStr(N+1)+"_"+NbrToStr(lc)+"_"+NbrToStr(type)+".txt";
// 	ofstream output(filename.c_str());
// 	if (!output){
// 		cerr<<"Output file cannot be created"<<endl;
// 		exit(1);
// 	}
// 	else{
// 		output<<N<<" "<<iteration<<endl;
// 		for (int j=0;j<times.size()-1;j++){
// 			output<<times.at(j)<<" ";
// 		}
// 		output<<times.back()<<endl;
// 	}
// }
// 
// void Create_output_3(int type_mesh, Real lc, bool type, vector<Real> kappa, string DDM, int S, int extension, int iteration, vector<double> times){
// 	string filename="../output/output_DDM_"+NbrToStr(type_mesh)+"_"+NbrToStr(1)+"_"+NbrToStr(lc)+"_"+NbrToStr(type)+"_"+NbrToStr(kappa.at(0))+"_"+DDM+".txt";
// 	ofstream output(filename.c_str(),ios::app);
// 	if (!output){
// 		cerr<<"Output file cannot be created"<<endl;
// 		exit(1);
// 	}
// 	else{
// 		if (DDM=="none"){
// 			output<<iteration<<endl;
// 		}
// 		else{
// 			output<<iteration<<" "<<S<<" "<<extension<<endl;
// 		}
// // 		for (int j=0;j<times.size()-1;j++){
// // 			output<<times.at(j)<<" ";
// // 		}
// // 		output<<times.back()<<endl;
// 	}
// }
// 
// void Create_output_solution(int type_mesh,Real lc, bool propagatif, vector<Real> kappa,vectC U, const mesh1D& Gamma, const vectR3& nodes){
// 	string filename="../solution/solution_"+NbrToStr(type_mesh)+"_"+NbrToStr(1)+"_"+NbrToStr(lc)+"_"+NbrToStr(propagatif)+"_"+NbrToStr(kappa.at(0))+".txt";
// 	ofstream solution(filename.c_str());
// 	if (!solution){
// 		cerr<<"Output file cannot be created"<<endl;
// 		exit(1);
// 	}
// 	else{
// 		for(int k=0; k<nodes.size(); k++){
// 			const N2& I = Gamma[k];
// 			
// 			solution << nodes[I[0]] <<"  "<<abs(U[k])<<endl;
// 
// 		}
// 		
// 	}
// 	
// }
#endif