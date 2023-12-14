#include "utilities.hpp"

using namespace std;

using namespace UTL;

void Initialization_for_N(int argc, char** argv, unsigned int &N){
    string nn = "-n";

    if(argc < 2) cout<<"Inserted by default N = "<<DEFAULT_N<<endl;
    else{
        for(int i= 0; i < argc; i++)
        {
            if(nn.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    N = std::stoi(argv[i+1]);
                }catch(exception){
                    cout<<"Please, insert a number after -n"<<endl<<"Inserted by default N = 100"<<endl;
                    exit(1);
                }
            }
            else if(nn.compare(argv[i]) == 0 && ( i+1 == argc )) {
                cout<<"Please, insert a number after -n"<<endl<<"Inserted by default N = 100"<<endl;
                exit(1);
            }
        } 
    } 
}


template<class SpMat>
void saveMatrixOnFile(SpMat A, std::string fileName){
    std::ofstream file;
    file.open(fileName, std::ofstream::trunc);
    file<<A.rows()<<" "<<A.cols()<<" "<<A.nonZeros()<<std::endl;

    for(size_t i = 0; i < A.rows(); i++){
        std::vector<size_t> row = A.nonZerosInRow(i);
        for(const auto& j : row){
            file<<i<<" "<<j<<" "<<A.coeffRef(i,j)<<std::endl;
        }
    }
    file.close();
}

template<class Vector>
void saveVectorOnFile(Vector f, std::string fileName){
    std::ofstream file;
    file.open(fileName, std::ofstream::trunc);
    file<<f.size()<<std::endl;

    for(size_t i = 0; i < f.size(); i++){
        file<<f[i]<<std::endl;
    }
    file.close();
}

std::vector<double> formatVector(std::vector<double> &in, AMG::Domain &domain){
    std::vector<double> temp;
    for(size_t i = 0; i < domain.N(); i++){
        temp.push_back(in[domain.mask(i)]);
    }
    return temp;
}