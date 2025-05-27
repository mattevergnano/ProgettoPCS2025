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
   
	int CreateSolid(PlatonicSolids& solido){
        //classify type of polyedra
        if(solido.p == 3 && solido.q == 3){
            cout << "Tetraedro" << endl;
            solido.NumCells0Ds = 4;
			solido.NumCells1Ds=6;
			solido.NumCells2Ds = 4;
            solido.Cells0DsId = vector<unsigned int>(solido.NumCells0Ds,0);
			solido.Cells1DsId = vector<unsigned int>(solido.NumCells1Ds,0);
			solido.Cells2DsId = vector<unsigned int>(solido.NumCells2Ds,0);
			for(unsigned int i=0;i<solido.NumCells0Ds;i++)
                solido.Cells0DsId[i] = i;
            solido.Cells0DsCoordinates = MatrixXd::Zero(3,solido.NumCells0Ds);
            solido.Cells0DsCoordinates << 1/sqrt(3),-1/sqrt(3),-1/sqrt(3),1/sqrt(3),
                                        1/sqrt(3),-1/sqrt(3),1/sqrt(3),-1/sqrt(3),
                                        1/sqrt(3),1/sqrt(3),-1/sqrt(3),-1/sqrt(3);
            solido.Cells1DsExtrema = MatrixXi::Zero(2,solido.NumCells1Ds);
			solido.Cells1DsExtrema << 0,0,0,1,1,2,
                                      1,2,3,2,3,3;
            cout << solido.Cells1DsExtrema << endl;
			for(unsigned int j=0; j<solido.NumCells2Ds; j++)
				solido.Cells2DsId[j] = j;
			solido.Cells2DsVertices = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
            solido.Cells2DsEdges = {{0, 1, 3}, {0, 2, 4}, {1, 2, 5}, {3, 4, 5} };
            solido.Cells2DsNumEdges = {3, 3, 3, 3};
        }else if(solido.p == 4 && solido.q == 3){
            cout << "Cubo" << endl;
			solido.NumCells0Ds=8;
			solido.NumCells1Ds=12;	
            solido.NumCells2Ds = 6;	
            solido.Cells2DsId = vector<unsigned int>(solido.NumCells2Ds,0);			
			solido.Cells1DsExtrema={{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
			for(unsigned int j=0; j<solido.NumCells2Ds; j++)
				solido.Cells2DsId[j] = j;
			solido.Cells2DsVertices = {{0,1,2,3}, {4,5,6,7}, {0,4,5,1}, {1,5,6,2}, {2,6,7,3}, {3,7,4,0}};
            solido.Cells2DsEdges = {{0,1,2,3}, {4,5,6,7}, {8,4,9,0}, {9,5,10,1}, {10,6,11,2}, {11,7,8,3}};
			solido.Cells2DsNumEdges = {4, 4, 4, 4, 4, 4};
        }else if(solido.p == 3 && solido.q == 4){
            cout << "Ottaedro" << endl;
            solido.NumCells0Ds = 6;
			solido.NumCells1Ds=12;
			solido.NumCells2Ds = 8;
            solido.Cells0DsId = vector<unsigned int>(solido.NumCells0Ds,0);
            solido.Cells1DsId = vector<unsigned int>(solido.NumCells1Ds,0);
			for(unsigned int i=0;i<solido.NumCells0Ds;i++)
                solido.Cells0DsId[i] = i;
            solido.Cells0DsCoordinates = MatrixXd::Zero(solido.NumCells0Ds,3);
            solido.Cells0DsCoordinates << 1.0,0.0,0.0,
                                         -1.0,0.0,0.0,
                                          0.0,1.0,0.0,
                                          0.0,-1.0,0.0,
                                          0.0,0.0,1.0,
                                          0.0,0.0,-1.0;
			solido.Cells1DsExtrema ={{4,0},{4,1},{4,2},{4,3},{5,0},{5,1},{5,2},{5,3},{0,1},{1,2},{2,3},{3,0}}; 
            solido.Cells2DsId = vector<unsigned int>(solido.NumCells2Ds, 0);
            for (unsigned int j = 0; j < solido.NumCells2Ds; j++)
                  solido.Cells2DsId[j] = j;
            solido.Cells2DsVertices = {{4, 0, 2}, {4, 2, 1}, {4, 1, 3}, {4, 3, 0}, {5, 0, 2}, {5, 2, 1}, {5, 1, 3}, {5, 3, 0} };
            solido.Cells2DsEdges = {{0, 2, 10}, {2, 1, 9}, {1, 3, 8}, {3, 0, 11},{4, 6, 10}, {6, 5, 9}, {5, 7, 8}, {7, 4, 11}};
            solido.Cells2DsNumEdges = {3, 3, 3, 3, 3, 3, 3, 3};
        }else if(solido.p == 5 && solido.q == 3){
            cout << "Dodecaedro" << endl;
			solido.NumCells0Ds=20;
			solido.NumCells1Ds=30;
			solido.NumCells2Ds = 12;
			solido.Cells1DsExtrema={{0, 8}, {0, 10}, {0, 16},{1, 9}, {1, 11}, {1, 18},{2, 10}, {2, 12}, {2, 17},{3, 11}, {3, 13}, {3, 19},
			                        {4, 8}, {4, 14}, {4, 16},{5, 9},{5, 14}, {5, 20},
									{6, 12}, {6, 15}, {6, 17},{7, 13}, {7, 15}, {7, 20},{8, 14}, {9, 14},
                                    {10, 12}, {11, 13},{16, 17}, {18, 19},{18, 20}, {19, 20}};
			solido.Cells2DsId = vector<unsigned int>(solido.NumCells2Ds, 0);
            for (unsigned int j = 0; j < solido.NumCells2Ds; j++)
                  solido.Cells2DsId[j] = j;						
			solido.Cells2DsVertices = {{0,8,4,14,5}, {0,10,2,12,8}, {1,9,5,14,4}, {1,11,3,13,9}, {2,10,0,16,17}, {2,17,6,12,10}, 
			                           {3,11,1,18,19}, {3,19,7,13,11}, {4,8,12,6,15}, {5,9,13,7,20}, {6,17,16,0,8}, {7,20,18,1,9}};
            solido.Cells2DsEdges = {{0,12,13,24,16}, {1,6,7,26,0}, {3,15,16,25,12}, {4,9,10,27,3}, {1,6,28,29,2}, {7,18,19,20,28}, {4,5,30,31,29}, 
			                        {10,11,23,22,27}, {12,24,13,19,21}, {15,25,10,22,32}, {18,28,2,0,12}, {23,31,5,3,15}};
            solido.Cells2DsNumEdges = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};						
        }else if(solido.p == 3 && solido.q == 5){
            cout << "Icosaedro" << endl;
            double phi = (1+sqrt(5))/2;
            solido.NumCells0Ds = 12;
			solido.NumCells1Ds=30;
			solido.NumCells2Ds = 20;
            solido.Cells0DsId = vector<unsigned int>(solido.NumCells0Ds,0);
			solido.Cells1DsId = vector<unsigned int>(solido.NumCells1Ds,0);
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
		    solido.Cells1DsExtrema= {{0,1}, {0,4}, {0,5}, {0,8}, {0,9},{1,4}, {1,6}, {1,8}, {1,9},
			                        {2,3}, {2,4}, {2,5}, {2,10}, {2,11},{3,6}, {3,7}, {3,10}, {3,11},
                                    {4,5}, {4,8}, {4,10},{5,9}, {5,11},{6,7}, {6,8}, {6,10},
                                    {7,9}, {7,11},{8,10}, {9,11}};
			solido.Cells2DsId = vector<unsigned int>(solido.NumCells2Ds, 0);
            for(unsigned int j = 0; j < solido.NumCells2Ds; j++)
                 solido.Cells2DsId[j] = j;
		    solido.Cells2DsVertices = {{0,1,8}, {0,8,4}, {0,4,5}, {0,5,9}, {0,9,1},{1,6,8}, {1,9,7}, {1,7,6}, {2,3,10}, {2,10,4},{2,4,5}, {2,5,11}, {2,11,3}, {3,6,10}, {3,7,6},
                                      {3,11,7}, {4,8,10}, {5,9,11}, {6,10,8}, {7,11,9}};
		    solido.Cells2DsEdges = {{0, 3, 7}, {3, 1, 19}, {1, 2, 18}, {2, 4, 21}, {4, 0, 8},{5, 6, 7}, {8, 22, 23}, {6, 23, 24}, {9, 10, 26}, {10, 11, 20},
                                    {11, 18, 19}, {12, 13, 22}, {13, 9, 28}, {14, 15, 26}, {15, 24, 25}, {13, 25, 27}, {19, 20, 28}, {21, 22, 29}, {6, 20, 24},{25, 27, 29}};
			solido.Cells2DsNumEdges = {3, 3, 3, 3};	
            solido.Cells2DsNumEdges = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};		
		
        
    return 0;
    };
	
	
}
}