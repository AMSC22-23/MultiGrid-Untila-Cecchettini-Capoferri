#include "allIncludes.hpp"

void Utils::Initialization_for_N(int argc, char** argv, size_t &N, double &alpha, double &width){
    //@note: the `GetPot` library automatize this task, do not reinvent the wheel if not
    //       necessary (or explicitly asked). However it is still a very nice exercise!
    std::string nn = "-n";
    std::string aa = "-a";
    std::string ww = "-w";
    alpha = DEFAULT_ALPHA;
    N = DEFAULT_N;
    width = DEFAULT_WIDTH;

    if(argc < 2)
    {
        std::cout<<"Inserted by default N = "<<DEFAULT_N<<std::endl;
        std::cout<<"Inserted by default alpha = "<<DEFAULT_ALPHA<<std::endl;
        std::cout<<"Inserted by default width = "<<DEFAULT_WIDTH<<std::endl;
    } 
    else{
        for(int i= 0; i < argc; i++)
        {
            if(nn.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    N = std::stoi(argv[i+1]);
                    std::cout<<"Inserted N = "<<N<<std::endl;
                //@note:it is generally recommended to catch exceptions by reference rather than by value. 
                // Catching exceptions by reference has several advantages:
                // - Avoid Slicing: Catching by reference prevents object slicing
                // - Avoid Object Copying: Catching by reference avoids unnecessary object copying
                // - Polymorphism: Catching by reference allows polymorphic behavior
                }catch(std::exception&){
                    std::cout<<"Please, insert a number after -n"<<std::endl;
                    std::exit(1);
                }
                
            }
            else if(aa.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    alpha = std::atof(argv[i+1]);
                    std::cout<<"Inserted alpha = "<<alpha<<std::endl;
                }catch(std::exception&){
                    std::cout<<"Please, insert a double after -a"<<std::endl;
                    std::exit(1);
                }
            }
            else if(ww.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    width = std::atof(argv[i+1]);
                    std::cout<<"Inserted width = "<<width<<std::endl;
                }catch(std::exception&){
                    std::cout<<"Please, insert a double after -w"<<std::endl;
                    std::exit(1);
                }
            }
            else if(nn.compare(argv[i]) == 0 && ( i+1 == argc )) {
                std::cout<<"Please, insert something"<<std::endl;
                exit(1);
            }
        } 
    } 
}