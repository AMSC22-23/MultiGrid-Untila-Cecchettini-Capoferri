#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

/*the main goal of this parte is to write a function in order to discretize the square domain(esiest case).
for this case we describe the domain as 4 limits...for the most general case would be better to implement a classe DOMAIN*/


//definition of my_tuple where i save my coordinates
typedef  std::vector < std::tuple<double, double>> my_tuple;
/*paasserò un puntatore ad un vettore inizializzato vuoto di tuple*/
void Domain_qnt(my_tuple* list, int N, double B_left, double B_right, double B_top, double B_down) {

    double h;

    //trova h
    h = std::abs((B_right - B_left)/(std::round(sqrt(N))-1));

    //discretizzo il dominio per riga
    
    for(double i=B_top; i>=B_down; i -= h){
        for(double j= B_right; j<= B_left; j += h){
            list->push_back(std::make_tuple(i,j));
        }
    }
}


int main(){

my_tuple list;
Domain_qnt(&list, 576,500,450,200,150);

//stampo i tuple cioe le quantizzazioni del mio dominio
for(int k=0; k < list.size() ; k++){
        std::cout<< "(" << std::get<0>(list[k]) << " - " << std::get<1>(list[k]) << ")" <<std::endl;
}

    return 0;
}








