#include "allIncludes.hpp"

void Utils::Initialization_for_N(int argc, char** argv, size_t &N, double &alpha, double &width, int &level, int &functions_to_test){
    std::string nn = "-n";
    std::string aa = "-a";
    std::string ww = "-w";
    std::string lv = "-ml";
    std::string help = "--help";
    std::string fti = "-test";
    alpha = DEFAULT_ALPHA;
    N = DEFAULT_N;
    width = DEFAULT_WIDTH;
    level = DEFAULT_LEVEL;
    functions_to_test = DEFAULT_TEST;

    if(argc < 2)
    {
        std::cout<<"Inserted by default N = "<<DEFAULT_N<<std::endl;
        std::cout<<"Inserted by default alpha = "<<DEFAULT_ALPHA<<std::endl;
        std::cout<<"Inserted by default width = "<<DEFAULT_WIDTH<<std::endl;
        std::cout<<"Inserted by default multigrid level = "<<DEFAULT_LEVEL<<std::endl;
        std::cout<<"Inserted by default test number "<<DEFAULT_TEST<<std::endl;
    } 
    else{
        for(int i= 0; i < argc; i++)
        {
            if(nn.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    N = std::stoi(argv[i+1]);
                    std::cout<<"Inserted N = "<<N<<std::endl;
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
            else if(lv.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    level = std::stoi(argv[i+1]);
                    std::cout<<"Inserted level = "<<level<<std::endl;
                }catch(std::exception&){
                    std::cout<<"Please, insert a number after -ml"<<std::endl;
                    std::exit(1);
                }
            }
            else if(fti.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    functions_to_test = std::stoi(argv[i+1]);
                    std::cout<<"Inserted test number = "<<functions_to_test<<std::endl;
                }catch(std::exception&){
                    std::cout<<"Please, insert a double after -test"<<std::endl;
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
            else if(help.compare(argv[i]) == 0)
            {
                std::cout << "Usage: ./Multigrid [OPTIONS]\n"<<std::endl
                << "Options:"<<std::endl
                << "  -n, insert number of spaces"<<std::endl
                << "  -a, specifies differential constant"<<std::endl
                << "  -w, insert the Width of the rectangle domain"<<std::endl
                << "  -ml, insert multigrid level"<<std::endl
                << "  -test, insert type of function in input to test it"<<std::endl
                << "  --help, Display this help message"<<std::endl;
                std::exit(1);
            }
            else if(nn.compare(argv[i]) == 0 && ( i+1 == argc )) {
                std::cout<<"Please, insert something"<<std::endl;
                exit(1);
            }
        } 
    } 
}


void Utils::init_test_functions(std::function<double(const double, const double)> &f, std::function<double(const double, const double)> &g, int i)
{
    using function_to_c = std::function<double(const double, const double)>;
    std::array< function_to_c,  6> functions_to_choose_to_test ={
        [] (const double x, const double y) {  return -5.0 * exp(x) * exp(-2.0 * y); }, 
        [] (const double x, const double y) { return exp(x) * exp(-2.0 * y); }, 
        [] (const double x, const double y) { return sin(30. * sqrt(x * x + y * y)); }, 
        [] (const double x, const double y) { double k = 30.;
            double r = sqrt(x*x + y*y);
            return -k*(cos(k * r) / r - k*sin(k * r));}, 
        [] (const double x, const double y) { return 1.; },
        [] (const double x, const double y) { return 0.; }
    };  

    switch (i)
    {
    case 1:
        f = functions_to_choose_to_test[0];
        g = functions_to_choose_to_test[1];
        break;
    case 2:
        f = functions_to_choose_to_test[3];
        g = functions_to_choose_to_test[2];
        break;
    case 3:
        f = functions_to_choose_to_test[4];
        g = functions_to_choose_to_test[5];
        break;
    default:
        f = functions_to_choose_to_test[4];
        g = functions_to_choose_to_test[5];
        break;
    }

}
