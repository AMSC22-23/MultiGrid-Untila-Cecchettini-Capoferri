#include "domain.hpp"
#include "linear_system.hpp"
#include "Main.hpp"

// to set, for the discretization of the domain, the N: -n number_of_elements

void Initialization_for_N(int argc, char** argv, unsigned int &N){
    string nn = "-n";

    if(argc < 2) cout<<"Inserted by default N = "<<DEFAULT_N<<endl;
    else{
        for(int i= 0; i < argc; i++)
        {
            if(nn.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    N = stoi(argv[i+1]);
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

using namespace std;

int main(int argc, char** argv)
{
    unsigned int N;
    Initialization_for_N(argc, argv, N);





    return 0;
}