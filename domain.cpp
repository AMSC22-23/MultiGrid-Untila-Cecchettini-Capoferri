#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

/*the main goal of this parte is to write a function in order to discretize the square domain(esiest case).
for this case we describe the domain as 4 limits...for the most general case would be better to implement a classe DOMAIN*/


//definition of my_tuple where i save my coordinates
typedef  std::vector < std::tuple<double, double>> my_tuple;
/*paasserò un puntatore ad un vettore inizializzato vuoto di tuple*/

//funzione per la stampa del domizio quantizzato
void Print_matrix(my_tuple &tp){
    for(int i=0; i<tp.size() ;i++){

        if(std::get<1>(tp[i])!=std::get<1>(tp[i-1]) && i!=0)
            std::cout<<std::endl<<std::endl;

         std::cout<< "(" << std::get<0>(tp[i]) <<":"<< std::get<1>(tp[i]) << ")"<<"  ";
         
    }
}

void Domain_qnt(my_tuple &list, int N, double B_left, double B_right, double B_top, double B_down) {

    double h;

    //trova h
    h = std::abs((B_right - B_left)/(std::round(sqrt(N))-1));

    //discretizzo il dominio per riga
    if(list.size()!=0)
        return;
    for(double j=B_top; j>=B_down; j -= h){
        for(double i= B_left; i<= B_right; i += h){
            list.push_back(std::make_tuple(i,j));
        }
    }
}


int main(int argc, char** argv){

my_tuple list;
Domain_qnt(list,156,-10.0,10,-10,-30.0);

for(int i=0;i<argc; i++){
    std::cout<<argv[i]<<std::endl;
}

//stampo i tuple cioe le quantizzazioni del mio dominio
/*for(size_t k=0; k < list.size() ; k++){
        std::cout<< "(" << std::get<0>(list[k]) << " - " << std::get<1>(list[k]) << ")" <<std::endl;
}*/

//stampo i tuple a matrice
Print_matrix(list);
    return 0;
}








