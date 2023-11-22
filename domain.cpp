#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

/*the main goal of this parte is to write a function in order to discretize the square domain(esiest case).
for this case we describe the domain as 4 limits...for the most general case would be better to implement a classe DOMAIN*/

class Domain{
        public:
        double BLeft,BRight, BTop, BDown;
        Domain( double  B_Left,double B_Right,double B_Top,double B_Down)
        : BLeft(B_Left),BRight(B_Right),BTop(B_Top), BDown(B_Down) {}
};

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

void Domain_qnt(my_tuple &list, int N, Domain Dmn) {

    double h,k;

    //trova h
    h = std::abs((Dmn.BRight - Dmn.BLeft)/(std::round(sqrt(N))-1));
    k = std::abs((Dmn.BTop - Dmn.BDown)/(std::round(sqrt(N))-1));

    //discretizzo il dominio per riga
    if(list.size()!=0)
        return;
    for(double j=Dmn.BTop; j>=Dmn.BDown; j -= k){
        for(double i= Dmn.BLeft; i<= Dmn.BRight; i += h){
            list.push_back(std::make_tuple(i,j));
        }
    }
}


int main(int argc, char** argv){

int N=56;
my_tuple list;
Domain DMN(-100.0,10.0,-10.0,-30.0);

Domain_qnt(list,N,DMN);

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








