#include <iostream>
#include <vector>
#define DEFAULT_N 100

class DOMAIN {
    public:
        DOMAIN(int argc, char** argv) {
            using namespace std;
            // below... 
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

    void printN(){
        std::cout<<"And N is = "<< N<<std::endl;
    }

    private:
        std::size_t N = DEFAULT_N; // default value
        std::vector<int> Boundary = {0, 1000, 0, 1000};
};