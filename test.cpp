#include <stdio.h>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#define epsilon 0.25


/// @brief 
/// @tparam T 
/// @param R is the row of matrix A, i.e. node i-th
/// @return a vector with position i-th equal to 1 if strongly connected, 0 otherwise

template <typename T>
std::vector<T> valueStrongConnection(std::vector<T> R, int elementI, int multi){
    // compute max
    T max = 0;
    
    for(int i=0 ; i < R.size(); i++) 
    {
        if(i != elementI && std::abs(R[i]) > max)
            max = std::abs(R[i]);
    }

    // compute strong connections
    std::vector<T> Ret;

    for(int i=0; i < R.size(); i++)
    {
        if(i != elementI && std::abs(R[i]) >= epsilon*max)
        {
            Ret.push_back(1*multi);
            
        }else if(i == elementI || multi == 1){
            Ret.push_back(0);
        }else{
            Ret.push_back(1);
        }
    }
    return Ret;
}


/// @brief this function evaluate the importance measure  of node i-th
/// @tparam T 
/// @param NodeToEval correspond to row of matrix A, is the node i-th
/// @param allNodes[i] is set to 2 if is in F and 1 if it is still undecided. 0 if in course 
/// @return 
template <typename T>
T evaluateNode(std::vector<T> NodeToEval, int NodeIndex, std::vector<T> allNodes ){
    std::vector<T> R = valueStrongConnection(NodeToEval, NodeIndex,1); // is a vector with position i-th equal to 1 if strongly connected, 0 otherwise

    if (allNodes.size() != R.size()) {
        std::cerr << "Vectors must be of the same length!" << std::endl;
        return -1;
    }

    return inner_product(R.begin(), R.end(), allNodes.begin(), 0);
}

int getRandomInit(int max){
    std::random_device rd;

    // Initialize a random number generator with the random device
    std::mt19937 gen(rd());
    // Create a uniform distribution within the specified range
    std::uniform_int_distribution<> distr(0, max);

    // Generate a random number within the range
    return distr(gen);
}

template <typename T>
void printVector(std::vector<T> result){
    std::cout << "Vector elements: ";
    for (int i : result) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

// matrix A is the correspondence of the connections of a graph, aij represent the edje between node i and j. So, if the connection exist, 
// it will be represented by a value different from 0, otherwise the two nodes are disconnected.
// So the "primi vicini" are represented by all values in a specific row differents from zero
// all spd matrix are all square graph and the AMG method correspond to geometric multigrid

template <typename T>
std::vector<T> AMG(std::vector<std::vector<T>> &A){
    std::vector<T> R(A.size(), 1); // numero righe = numero nodi
    int index = getRandomInit(A.size() - 1 );
    bool GoOn;
    do{
        GoOn = false;
        std::vector<T> Tem = valueStrongConnection(A[index],index,2);
        for(int i = 0; i<R.size(); i++)
        {
            if(R[i] == 1)
            {
                R[i] = Tem[i];
            }
        }
        //printVector(R);

        int tempMax = 0;
        
        for(int i = 0; i< R.size(); i++)
        {
            if(R[i] == 1)
            {
                int tt = evaluateNode(A[i], i, R);

                if(tt > tempMax)
                {
                    GoOn = true;
                    tempMax = tt;
                    index = i;
                }
            }
        }
    }while(GoOn);
    return R;
    // A.size() is the row number
}




int main() {

    using namespace std;
    
    vector<vector<int>> A = {{4, -1, 0, 0, -1},{-1, 4, -1,0, 0}, {0, -1, 4, -1, 0},{0, 0, -1, 4, -1},{-1, 0, 0, -1, 4}};

    int N = 1024;

    vector<vector<int>> Matrix(N, vector<int>(N, 0)); // Inizializza tutti gli elementi a 0

    // Esempio: Riempimento della matrice per renderla simmetrica
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if(i == j){
                Matrix[i][i] = 5;
            }else if (j == (i +1)){
                Matrix[i][j] = -1;
            }else if (j == (i -1)){
                Matrix[i][j] = -1;
            }
            
        }
    }

    Matrix[5][5] = -3;
    Matrix[10][7] = -4;
    Matrix[10][8] = -4;

    auto start = std::chrono::high_resolution_clock::now();
    vector<int> result = AMG(Matrix);
    //printVector(result);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> solve_time = end - start;
    std::cout<<"Solving elapsed time: "<<solve_time.count()<<" seconds"<<std::endl;
    
    return 0;
}