#include <iostream>
#include <vector>

class DOMAIN1 {
    public:
        DOMAIN1(int argc, char** argv) {
            using namespace std;
        // below... 
        std::string nn = "-n";
        if(argc < 2) 
        {
            N = 100;
            cout<<"Inserted by default N = 100"<<endl;
        }
        else{
            for(int i= 0; i < argc; i++)
            {
                if(nn.compare(argv[i]) == 0 && ( i+1 < argc ) )
                {
                    try{
                        N = std::stoi(argv[i+1]);
                    }catch(std::exception){
                        cout<<"Please, insert a number after -n"<<endl;
                    }
                }
                else if(nn.compare(argv[i]) == 0 && ( i+1 == argc )) cout<<"Please, insert a number after -n"<<endl;
            }
        } 
    }

    void printN(){
        std::cout<<"And N is = "<< N<<std::endl;
    }

    private:
        std::size_t N = 0;
        std::vector<int> Boundary = {0, 1000, 0, 1000};
};