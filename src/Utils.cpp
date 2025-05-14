#include<iostream>
#include<Utils.hpp>
#include"PlatonicSolids.hpp"

namespace PlatonicLibrary{
    int ImportValue(int argc, char  *argv[],PlatonicSolids solido){

        string str = ""; //remove Project name
        str = argv[0];
        solido.p = stoi(argv[1]); //save p
        solido.q = stoi(argv[2]); //save q
        if(argc == 4){ //only one value (b)
            solido.b = stoi(argv[3]);
        } else if(argc == 5){ //two inputs
            if(stoi(argv[3]) == 0){
                solido.b = stoi(argv[4]); //only b if the third input is 0
            }
            else{ //both b and c
                solido.b = stoi(argv[3]);
                solido.c = stoi(argv[4]);
            };
            
        } else{
            cerr << "Not Valid Input b or c" << endl; //every other case
            return 1;
        }
        //check for p and q
        if(solido.p < 3 || solido.q <3){
            cerr << "Not Valid Input p or q" << endl;
            return 1;
        }
        cout << "p: " << solido.p << "\nq: " << solido.q << "\nb: " << solido.b << "\nc: " << solido.c << endl;
        return 0;
    }
}