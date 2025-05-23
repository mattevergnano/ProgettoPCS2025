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
			solido.Cells1DsVertices={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        }else if(solido.p == 4 && solido.q == 3){
            cout << "Cubo" << endl;
			solido.NumCells0Ds=8;
			solido.NumCells1Ds=12;
			
			solido.Cells1DsVertices={{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
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
			solido.Cells1DsVertices ={{4,0},{4,1},{4,2},{4,3},{5,0},{5,1},{5,2},{5,3},{0,1},{1,2},{2,3},{3,0}}; 
        }else if(solido.p == 5 && solido.q == 3){
            cout << "Dodecaedro" << endl;
			solido.NumCells0Ds=20;
			solido.NumCells1Ds=30;
			solido.Cells1DsVertices={{0, 8}, {0, 10}, {0, 16},{1, 9}, {1, 11}, {1, 18},{2, 10}, {2, 12}, {2, 17},{3, 11}, {3, 13}, {3, 19},
			                        {4, 8}, {4, 14}, {4, 16},{5, 9},{5, 14}, {5, 20},
									{6, 12}, {6, 15}, {6, 17},{7, 13}, {7, 15}, {7, 20},{8, 14}, {9, 14},
                                    {10, 12}, {11, 13},{16, 17}, {18, 19},{18, 20}, {19, 20}};
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
		    solido.Cells1DsVertices= {{0,1}, {0,4}, {0,5}, {0,8}, {0,9},{1,4}, {1,6}, {1,8}, {1,9},
			                        {2,3}, {2,4}, {2,5}, {2,10}, {2,11},{3,6}, {3,7}, {3,10}, {3,11},
                                    {4,5}, {4,8}, {4,10},{5,9}, {5,11},{6,7}, {6,8}, {6,10},
                                    {7,9}, {7,11},{8,10}, {9,11}};
        }
        return 0;
    };
}