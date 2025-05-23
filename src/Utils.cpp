#include<iostream>
#include<math.h>
#include<tgmath.h>
#include<Utils.hpp>
#include"PlatonicSolids.hpp"
#include "Eigen/Eigen"
#include<string>
#include<vector>

using namespace std;
using namespace Eigen;

namespace PlatonicLibrary{
    int ImportValue(int argc, char  *argv[],PlatonicSolids& solido){

        string str = ""; //remove Project name
        str = argv[0];
        solido.p = stoi(argv[1]); //save p
        solido.q = stoi(argv[2]); //save q
        if(argc == 4){ //only one value (b)
            solido.b = stoi(argv[3]);
        } else if(argc == 5){ //four inputs
            if(stoi(argv[3])!=0 && stoi(argv[4])==0){ //input like p q b 0
                solido.b = stoi(argv[3]);
            }else if(stoi(argv[3]) == 0 && stoi(argv[4])!=0){ //input like p q 0 b
                solido.b = stoi(argv[4]); //only b if the third input is 0
            }else if(stoi(argv[3]) == 0 && stoi(argv[4]) == 0){ //input like p q b c
                cerr << "Not valid input, b and c are 0" << endl;
                return 1;
            }else{ //both b and c
                solido.b = stoi(argv[3]);
                solido.c = stoi(argv[4]);
            };
        } else{
            cerr << "Not Valid input" << endl; //every other case
            return 1;
        }
        //check for p and q
        if(solido.p < 3 || solido.q <3){
            cerr << "Not Valid input" << endl;
            return 1;
        }
        cout << "p: " << solido.p << "\nq: " << solido.q << "\nb: " << solido.b << "\nc: " << solido.c << endl;
        return 0;
    };
    int CreateSolid(PlatonicSolids solido){
        //classify type of polyedra
        if(solido.p == 3 && solido.q == 3){
            cout << "Tetraedro" << endl;
            solido.NumCells0Ds = 4;
            solido.Cells0DsId = vector<unsigned int>(solido.NumCells0Ds,0);
            solido.NumCells1Ds=6;
			for(unsigned int i=0;i<solido.NumCells0Ds;i++)
                solido.Cells0DsId[i] = i;
            solido.Cells0DsCoordinates = MatrixXd::Zero(solido.NumCells0Ds,3);
            solido.Cells0DsCoordinates << 0.0,0.0,1.0,
                                        -sqrt(3)/4,sqrt(3)/4,1.0/2,
                                        -sqrt(3)/4,-sqrt(3)/4,-(1.0/2),
                                        0.0,0.0,-1.0;
        }else if(solido.p == 4 && solido.q == 3){
            cout << "Cubo" << endl;
			solido.NumCells0Ds=8;
			solido.NumCells1Ds=12;
        }else if(solido.p == 3 && solido.q == 4){
            cout << "Ottaedro" << endl;
            solido.NumCells0Ds = 6;
			solido.NumCells1Ds=12;
            solido.Cells0DsId = vector<unsigned int>(solido.NumCells0Ds,0);
            for(unsigned int i=0;i<solido.NumCells0Ds;i++)
                solido.Cells0DsId[i] = i;
            solido.Cells0DsCoordinates = MatrixXd::Zero(solido.NumCells0Ds,3);
            solido.Cells0DsCoordinates << 1.0,0.0,0.0,
                                         -1.0,0.0,0.0,
                                          0.0,1.0,0.0,
                                          0.0,-1.0,0.0,
                                          0.0,0.0,1.0,
                                          0.0,0.0,-1.0;
        }else if(solido.p == 5 && solido.q == 3){
            cout << "Dodecaedro" << endl;
			solido.NumCells0Ds=20;
			solido.NumCells1Ds=30;
        }else if(solido.p == 3 && solido.q == 5){
            cout << "Icosaedro" << endl;
            double phi = (1+sqrt(5))/2;
            solido.NumCells0Ds = 12;
			solido.NumCells1Ds=30;
            solido.Cells0DsId = vector<unsigned int>(solido.NumCells0Ds,0);
            for(unsigned int i=0;i<solido.NumCells0Ds;i++)
                solido.Cells0DsId[i] = i;
            solido.Cells0DsCoordinates = MatrixXd::Zero(solido.NumCells0Ds,3);
            solido.Cells0DsCoordinates << 0.0,1.0,phi,
                                          0.0,-1.0,phi,
                                          0.0,1.0,-phi,
                                          0.0,-1.0,-phi,
                                          1.0,phi,0.0,
                                          -1.0,phi,0.0,
                                          1.0,-phi,0.0,
                                          -1.0,-phi,0.0,
                                          phi,0.0,1.0,
                                          -phi,0.0,1.0,
                                          phi,0.0,-1.0,
                                          -phi,0.0,-1.0;
        }
        return 0;
    };
}