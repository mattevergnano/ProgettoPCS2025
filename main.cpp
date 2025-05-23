#include<iostream>
#include "PlatonicSolids.hpp"
#include "Utils.hpp"

using namespace std;
using namespace PlatonicLibrary;

int main(int argc,char  *argv[]){
    cout << "ProgettoPCS2025" << endl;
    PlatonicSolids solido;
    ImportValue(argc,argv,solido);
    //CreateSolid(solido);
    return 0;   
}