#include<iostream>
#include<math.h>
#include<tgmath.h>
#include "Utils.hpp"
#include "PlatonicSolids.hpp"
#include "Eigen/Eigen"
#include<string>
#include<iomanip>
#include<fstream>
#include<vector>
#include<queue>
#include<unordered_map>
#include<set>

using namespace std;
using namespace Eigen;

bool comp(unsigned int a,unsigned int b) {
    return a < b;
}
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
	    } else if(argc == 7){
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
            solido.id_vertice1 = stoi(argv[5]);
		 	solido.id_vertice2 = stoi(argv[6]);
            if((stoi(argv[3]) == 0 && stoi(argv[4]) == 0) || (stoi(argv[5]) == 0 && stoi(argv[6])== 0)){
		        cerr << "Not Valid input" << endl; //every other case
		   		return 1;
            }
        }
        //check for p and q
        if(solido.p < 3 || solido.q <3){
            cerr << "Not Valid input" << endl;
            return 1;
        }
        return 0;
    }
  
	int CreateSolid(PlatonicSolids& solido){
        //classify type of polyhedra
        if(solido.p == 3 && solido.q == 3){
            //TETRAEDRO
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
			for(unsigned int j=0; j<solido.NumCells2Ds; j++)
				solido.Cells2DsId[j] = j;
			solido.Cells2DsVertices = MatrixXi::Zero(solido.p, solido.NumCells2Ds);
			solido.Cells2DsVertices <<  0, 0, 0, 1,  
                                        1, 1, 2, 2,  
                                        2, 3, 3, 3;  
			solido.Cells2DsEdges = MatrixXi::Zero(solido.p, solido.NumCells2Ds);		
            solido.Cells2DsEdges << 0, 0, 1, 3,  
                                    1, 2, 2, 4,  
                                    3, 4, 5, 5;  					   
            solido.Cells2DsNumEdges = VectorXi::Zero(solido.NumCells2Ds);
            solido.Cells2DsNumEdges << 3, 3, 3, 3;

            //creare un vettore di vettori, ogni vettore piccolo contiene gli id delle facce che sono adiacenti alla faccia, ovvero quelli che hanno un lato in comune
            solido.Cells2DsNeighborhood.reserve(solido.NumCells2Ds);
            for(unsigned int i : solido.Cells2DsId){
                vector<unsigned int> vettore;
                for(unsigned int latoi : solido.Cells2DsEdges.col(i)){
                    for(unsigned int j : solido.Cells2DsId){
                        for(unsigned int latoj : solido.Cells2DsEdges.col(j)){
                            if(i != j && latoi == latoj){
                                vettore.push_back(j);
                            }
                        }
                }
                }
                sort(vettore.begin(), vettore.end(), comp);
                solido.Cells2DsNeighborhood.insert(solido.Cells2DsNeighborhood.begin() + i,vettore);
            }

            solido.VerticeFaces.resize(solido.NumCells0Ds);// ho messo resize , al posto di reserve
            //riservo memoria
            for (auto& sottovettore : solido.VerticeFaces) {
                sottovettore.reserve(solido.q); // riserva spazio per 5 elementi per ogni sottovettore
            }
            for(unsigned int i = 0;i<solido.NumCells2Ds;i++){
                for(unsigned int j = 0;j<solido.q;j++){
                    solido.VerticeFaces[solido.Cells2DsVertices(j,i)].push_back(i);
                }
            }
            for(auto& el : solido.VerticeFaces){
                sort(el.begin(),el.end(),comp);
            }
			solido.Cells3DsVertices.resize(solido.NumCells0Ds);
			for(unsigned int i = 0; i<solido.NumCells0Ds;i++)
				solido.Cells3DsVertices[i]=i;
			solido.Cells3DsEdges.resize(solido.NumCells1Ds);
			for(unsigned int i = 0; i<solido.NumCells1Ds;i++)
				solido.Cells3DsEdges[i]=i;
			solido.Cells3DsFaces.resize(solido.NumCells2Ds);
			for(unsigned int i = 0; i<solido.NumCells2Ds;i++)
				solido.Cells3DsFaces[i]=i;
        }else if(solido.p == 4 && solido.q == 3){
            //CUBO, duale dell'ottaedro
            //creo una nuova istanza di solido, solido1, in cui genero l'ottaedro, e poi passo solido e solido1 a Duale per creare il cubo
			PlatonicSolids solido1;
            solido1.p=solido.q;
            solido1.q=solido.p;
            CreateSolid(solido1);
			DualPolyhedron(solido,solido1);
            return 0;
        }else if(solido.p == 3 && solido.q == 4){
            //OTTAEDRO
            solido.NumCells0Ds = 6;
			solido.NumCells1Ds = 12;
			solido.NumCells2Ds = 8;
            solido.Cells0DsId = vector<unsigned int>(solido.NumCells0Ds,0);
            solido.Cells1DsId = vector<unsigned int>(solido.NumCells1Ds,0);
			for(unsigned int i=0;i<solido.NumCells0Ds;i++)
                solido.Cells0DsId[i] = i;
            solido.Cells0DsCoordinates = MatrixXd::Zero(3,solido.NumCells0Ds);
            solido.Cells0DsCoordinates << 1.0,-1.0,0.0,0.0,0.0,0.0,
                                          0.0,0.0,1.0,-1.0,0.0,0.0,
                                          0.0,0.0,0.0,0.0,1.0,-1.0;
			solido.Cells1DsExtrema = MatrixXi::Zero(2,solido.NumCells1Ds);
			solido.Cells1DsExtrema << 0,0,0,0,1,1,1,1,2,2,3,3,
                                      2,3,4,5,2,3,4,5,4,5,4,5; 
            solido.Cells2DsId = vector<unsigned int>(solido.NumCells2Ds, 0);
            for (unsigned int j = 0; j < solido.NumCells2Ds; j++)
                  solido.Cells2DsId[j] = j;
			solido.Cells2DsVertices = MatrixXi::Zero(solido.p, solido.NumCells2Ds);
            solido.Cells2DsVertices << 0,0,0,0,1,1,1,1,
                                       2,3,2,3,2,3,3,2,
                                       4,4,5,5,4,4,5,5;
			solido.Cells2DsEdges = MatrixXi::Zero(solido.p, solido.NumCells2Ds);
            solido.Cells2DsEdges << 8,  1, 9,  11,  4, 10,  7, 9,
                                    2,  2,  3,  3,  8,  5,  5,  7,
                                    0,  10,  0,  1,  6,  6,  11,  4;
			solido.Cells2DsNumEdges = VectorXi::Zero(solido.NumCells2Ds);
            solido.Cells2DsNumEdges << 3, 3, 3, 3, 3, 3, 3, 3;

            solido.Cells2DsNeighborhood.reserve(solido.NumCells2Ds);
            for(unsigned int i : solido.Cells2DsId){
                vector<unsigned int> vettore;
                for(unsigned int latoi : solido.Cells2DsEdges.col(i)){
                    for(unsigned int j : solido.Cells2DsId){
                        for(unsigned int latoj : solido.Cells2DsEdges.col(j)){
                            if(i != j && latoi == latoj){
                                // vettore.insert(vettore.begin(),j);
                                vettore.push_back(j);
                            }
                        }
                }
                }
                sort(vettore.begin(), vettore.end(), comp);
                solido.Cells2DsNeighborhood.insert(solido.Cells2DsNeighborhood.begin() + i,vettore);
            }

            solido.VerticeFaces.resize(solido.NumCells0Ds); //un vettore per ogni vertice
            //riservo memoria
            for (auto& sottovettore : solido.VerticeFaces) {
                sottovettore.reserve(solido.q); // riserva spazio per 5 elementi per ogni sottovettore
            }
            for(unsigned int i = 0;i<solido.NumCells2Ds;i++){
                for(unsigned int j = 0;j<solido.p;j++){
                    solido.VerticeFaces[solido.Cells2DsVertices(j,i)].push_back(i); //al vettore corrispondente al vertice aggiungo la faccia i
                }
            }
            for(auto& el : solido.VerticeFaces){
                sort(el.begin(),el.end(),comp);
            }	
			solido.Cells3DsVertices.resize(solido.NumCells0Ds);
			for(unsigned int i = 0; i<solido.NumCells0Ds;i++)
				solido.Cells3DsVertices[i]=i;
			solido.Cells3DsEdges.resize(solido.NumCells1Ds);
			for(unsigned int i = 0; i<solido.NumCells1Ds;i++)
				solido.Cells3DsEdges[i]=i;
			solido.Cells3DsFaces.resize(solido.NumCells2Ds);
			for(unsigned int i = 0; i<solido.NumCells2Ds;i++)
				solido.Cells3DsFaces[i]=i; 
        }else if(solido.p == 5 && solido.q == 3){
            //DODECAEDRO, duale dell'icosaedro
            PlatonicSolids solido1;
            solido1.p=solido.q;
            solido1.q=solido.p;
            CreateSolid(solido1);
            DualPolyhedron(solido,solido1);
            return 0;			
        }else if(solido.p == 3 && solido.q == 5){
            //ICOSAEDRO
            double phi = (1+sqrt(5))/2;
            solido.NumCells0Ds = 12;
			solido.NumCells1Ds=30;
			solido.NumCells2Ds = 20;
            solido.Cells0DsId = vector<unsigned int>(solido.NumCells0Ds,0);
			solido.Cells1DsId = vector<unsigned int>(solido.NumCells1Ds,0);
            for(unsigned int i=0;i<solido.NumCells0Ds;i++)
                solido.Cells0DsId[i] = i;
            solido.Cells0DsCoordinates = MatrixXd::Zero(3,solido.NumCells0Ds);
            solido.Cells0DsCoordinates << 0.0,0.0,0.0,0.0,1.0,-1.0,1.0,-1.0,phi,-phi,phi,-phi,
                                          1.0,-1.0,1.0,-1.0,phi,phi,-phi,-phi,0.0,0.0,0.0,0.0,
                                          phi,phi,-phi,-phi,0.0,0.0,0.0,0.0,1.0,1.0,-1.0,-1.0;
            for(unsigned int i=0;i<solido.NumCells0Ds;i++){
                double x = solido.Cells0DsCoordinates(0,i);
                double y = solido.Cells0DsCoordinates(1,i);
                double z = solido.Cells0DsCoordinates(2,i);
                double norm = sqrt(x*x+y*y+z*z);
                solido.Cells0DsCoordinates(0,i) /= norm;
                solido.Cells0DsCoordinates(1,i) /= norm;
                solido.Cells0DsCoordinates(2,i) /= norm;
            };

            solido.Cells1DsExtrema = MatrixXi::Zero(2,solido.NumCells1Ds);
            solido.Cells1DsExtrema << 0,0,0,0,0,1,1,1,1,2,2,2, 2, 2,3,3, 3, 3,4, 4,4, 5,5,6, 6,6, 7,7, 8, 9,
                                      1,4,5,8,9,6,7,8,9,3,4,5,10,11,6,7,10,11,5,10,8,11,9,7,10,8,11,9,10,11;
		
			solido.Cells2DsId = vector<unsigned int>(solido.NumCells2Ds, 0);
            for(unsigned int j = 0; j < solido.NumCells2Ds; j++)
                 solido.Cells2DsId[j] = j;

			solido.Cells2DsVertices = MatrixXi::Zero(solido.p, solido.NumCells2Ds);
            solido.Cells2DsVertices << 0,0,0,0,0,1,1,1, 8, 8, 9, 9,7, 7, 6,4, 4, 2, 2, 2,
                                       1,1,4,4,5,9,6,6, 6, 4, 5, 7,6,11,10,5, 2, 5, 3, 3,
                                       8,9,8,5,9,7,7,8,10,10,11,11,3, 3, 3,2,10,11,11,10;


			solido.Cells2DsEdges = MatrixXi::Zero(solido.p, solido.NumCells2Ds);						
			solido.Cells2DsEdges <<  0, 0, 1, 1, 2, 8, 5, 5,25,20,22,27,23,26,24,18,10,11, 9, 9,
                                      7, 8, 3,18,22,27,23,25,24,19,21,26,14,17,16,11,12,21,17,16,
                                      3, 4,20, 2, 4, 6, 6, 7,28,28,29,29,15,15,14,10,19,13,13,12;				
		    solido.Cells2DsNumEdges = VectorXi::Zero(solido.NumCells2Ds);
            solido.Cells2DsNumEdges << 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3;

            solido.Cells2DsNeighborhood.reserve(solido.NumCells2Ds);
            for(unsigned int i : solido.Cells2DsId){
                vector<unsigned int> vettore;
                for(unsigned int latoi : solido.Cells2DsEdges.col(i)){
                    for(unsigned int j : solido.Cells2DsId){
                        for(unsigned int latoj : solido.Cells2DsEdges.col(j)){
                            if(i != j && latoi == latoj){
                                // vettore.insert(vettore.begin(),j);
                                vettore.push_back(j);
                            }
                        }
                }
                }
                sort(vettore.begin(), vettore.end(), comp);
                solido.Cells2DsNeighborhood.insert(solido.Cells2DsNeighborhood.begin() + i,vettore);
            }

            solido.VerticeFaces.resize(solido.NumCells0Ds); //un vettore per ogni vertice
            //riservo memoria
            for (auto& sottovettore : solido.VerticeFaces) {
                sottovettore.reserve(solido.q); // riserva spazio per 5 elementi per ogni sottovettore
            }
            for(unsigned int i = 0;i<solido.NumCells2Ds;i++){
                for(unsigned int j = 0;j<solido.p;j++){
                    solido.VerticeFaces[solido.Cells2DsVertices(j,i)].push_back(i); //al vettore corrispondente al vertice aggiungo la faccia i
                }
            }
            for(auto& el : solido.VerticeFaces){
                sort(el.begin(),el.end(),comp);
            }
            solido.Cells3DsVertices.resize(solido.NumCells0Ds);
			for(unsigned int i = 0; i<solido.NumCells0Ds;i++)
				solido.Cells3DsVertices[i]=i;
				solido.Cells3DsEdges.resize(solido.NumCells1Ds);
			for(unsigned int i = 0; i<solido.NumCells1Ds;i++)
				solido.Cells3DsEdges[i]=i;
			solido.Cells3DsFaces.resize(solido.NumCells2Ds);
			for(unsigned int i = 0; i<solido.NumCells2Ds;i++)
				solido.Cells3DsFaces[i]=i;   			
			}	
    return 0;
    };
	
	void FileCell0Ds(PlatonicSolids& solido){
		
		ofstream outputFile;
		outputFile.open("Cells0Ds.txt");
		if (!outputFile) {
			cerr << "Error opening file: Cells0Ds.txt" << endl;
		return;
		}

		outputFile << "ID    x              y              z" << endl;
		for (unsigned int  i = 0; i < solido.Cells0DsId.size(); ++i) {
			outputFile << i << "    "
				   << fixed << setprecision(8)
				   << solido.Cells0DsCoordinates(0, i) << "    "
				   << solido.Cells0DsCoordinates(1, i) << "    "
				   << solido.Cells0DsCoordinates(2, i)
				   << endl;
		}

		outputFile.close();
    }

	void FileCell1Ds(PlatonicSolids& solido){
		
		ofstream outputFile;
		outputFile.open("Cells1Ds.txt");
		if (!outputFile) {
			cerr << "Error opening file: Cells1Ds.txt" << endl;
        return;
		}

		outputFile << "ID\tInitial Vertices\tFinal Vertices" << endl;

    
		for (unsigned int i = 0; i < solido.Cells1DsId.size(); ++i) {
			outputFile << i << "\t"
                   << solido.Cells1DsExtrema(0, i) << "\t\t"
                   << solido.Cells1DsExtrema(1, i) << endl;
		}

		outputFile.close();
	}
	
	void FileCell2Ds(PlatonicSolids& solido){
		
		ofstream outputFile;
		outputFile.open("Cells2Ds.txt");
		if (!outputFile) {
			cerr << "Error opening file: Cells2Ds.txt" << endl;
        return;
		}
       
		outputFile << "ID\tNumVertices\tNumEdges\tVerticesIDs\t\tEdgesIDs" << endl;

		for (unsigned int i = 0; i < solido.Cells2DsId.size(); ++i) {
			outputFile << i << "\t"                              
                   << solido.Cells2DsVertices.rows() << "\t\t" 
                   << solido.Cells2DsNumEdges(i) << "\t\t";    

        
			for (int j = 0; j < solido.Cells2DsVertices.rows(); ++j) {
				outputFile << solido.Cells2DsVertices(j, i) << " ";
			}

			outputFile << "\t\t\t";  

        
			for (int j = 0; j < solido.Cells2DsNumEdges(i); ++j) {
				outputFile << solido.Cells2DsEdges(j, i) << " ";
			}

			outputFile << endl;
		}

        outputFile.close();
    }
	
    void FileCell3Ds(PlatonicSolids& solido){
		ofstream outputFile("Cells3Ds.txt");
		if (!outputFile) {
			cerr << "Error opening file: Cells3Ds.txt" << endl;
			return;
		}

		outputFile << "FinalVertices:\n";
		for (unsigned int i = 0; i < solido.Cells3DsVertices.size(); ++i) {
			outputFile << solido.Cells3DsVertices[i] << "\n";
		}

		outputFile << "Edges:\n";
		for (unsigned int i = 0; i < solido.Cells3DsEdges.size(); ++i) {
			outputFile << solido.Cells3DsEdges[i] << "\n";
		}

		outputFile << "Faces:\n";
		for (unsigned int i = 0; i < solido.Cells3DsFaces.size(); ++i) {
			outputFile << solido.Cells3DsFaces[i] << "\n";
		}

		outputFile.close();
    }

    int DualPolyhedron(PlatonicSolids& solido,PlatonicSolids& solido1){
			solido.NumCells0Ds = solido1.NumCells2Ds;
			solido.NumCells1Ds = solido1.NumCells1Ds;
			solido.NumCells2Ds = solido1.NumCells0Ds;
			solido.Cells0DsId = vector<unsigned int>(solido.NumCells0Ds,0);
			solido.Cells1DsId = vector<unsigned int>(solido.NumCells1Ds,0);
			solido.Cells2DsId = vector<unsigned int>(solido.NumCells2Ds,0);
			for(unsigned int i=0;i<solido.NumCells0Ds;i++)
				solido.Cells0DsId[i] = i;
			for(unsigned int i=0;i<solido.NumCells1Ds;i++)
				solido.Cells1DsId[i] = i;
			for(unsigned int i=0;i<solido.NumCells2Ds;i++)
				solido.Cells2DsId[i] = i;
			MatrixXd punti = MatrixXd::Zero(3,solido.NumCells0Ds);
			MatrixXi lati = MatrixXi::Zero(2,solido.NumCells1Ds);
			solido.Cells2DsVertices = MatrixXi::Zero(3, solido.NumCells2Ds);
			solido.Cells2DsEdges = MatrixXi::Zero(3, solido.NumCells2Ds);
			solido.Cells2DsNumEdges = VectorXi::Zero(solido.NumCells2Ds);
			// for (int i = 0; i < 10; ++i) {
			//     cout << "Faccia " << i << ": ";
			//     for (int j = 0; j < 3; ++j)
			//         cout << solido1.Cells2DsEdges(j, i) << " ";
			//     cout << endl;
			// }
			for(unsigned int j=0; j<solido1.NumCells2Ds;j++){
				int Pid = solido1.Cells2DsVertices(0,j);
				int Qid = solido1.Cells2DsVertices(1,j);
				int Sid = solido1.Cells2DsVertices(2,j);
				cout << Pid << " " << Qid << " " << Sid << endl;
				double Px = solido1.Cells0DsCoordinates(0,Pid);
				double Py = solido1.Cells0DsCoordinates(1,Pid);
				double Pz = solido1.Cells0DsCoordinates(2,Pid);
				double Qx = solido1.Cells0DsCoordinates(0,Qid);
				double Qy = solido1.Cells0DsCoordinates(1,Qid);
				double Qz = solido1.Cells0DsCoordinates(2,Qid);
				double Sx = solido1.Cells0DsCoordinates(0,Sid);
				double Sy = solido1.Cells0DsCoordinates(1,Sid);
				double Sz = solido1.Cells0DsCoordinates(2,Sid);
				punti(0,j) = (Px+Qx+Sx)/3;
				punti(1,j)= (Py+Qy+Sy)/3;
				punti(2,j)= (Pz+Qz+Sz)/3;
				double x = punti(0,j);
				double y = punti(1,j);
				double z = punti(2,j);
				double norm = sqrt(x*x+y*y+z*z);
				punti(0,j) /= norm;
				punti(1,j) /= norm;
				punti(2,j) /= norm;
				cout << punti(0,j) << " " << punti(1,j) << " " << punti(2,j) << endl;
			} 
			//accedo a Cells2dNeighborhood e collego ciascun vertice di solido ai vertici corrispondenti alle facce solido.adjacency a quello (inizio ad avere il duplicato dei lati)
			// solido.Cells1DsExtrema = MatrixXi::Zero(2,solido.NumCells1Ds);
			// solido.Cells1DsExtrema.setConstant(1000);
			unsigned int nlati = 0;
			for(unsigned int idfaccia = 0;idfaccia<solido.NumCells0Ds;idfaccia++){
				vector<unsigned int> vettore = solido1.Cells2DsNeighborhood[idfaccia];
				// cout << idfaccia << ": ";
				// for(auto& el:vettore){
				//     cout << el << " ";
				// }
				// cout << endl;
				if(solido1.Cells2DsNeighborhood[idfaccia].size()!=0){
					for(unsigned int adj = 0;adj<vettore.size();adj++){
						if(idfaccia<vettore[adj]){
							lati(0,nlati) = idfaccia;
							lati(1,nlati) = vettore[adj];
							nlati ++;
							// cout << nlati -1 << endl;
						}
					}
				}
			}
			solido.Cells0DsCoordinates = MatrixXd::Zero(3,solido.NumCells0Ds);
			solido.Cells1DsExtrema = MatrixXi::Zero(2,solido.NumCells1Ds);
			solido.Cells0DsCoordinates = punti;
			solido.Cells1DsExtrema = lati;

			solido.Cells2DsVertices = MatrixXi::Zero(solido.p, solido.NumCells2Ds);
			solido.Cells2DsEdges = MatrixXi::Zero(solido.p, solido.NumCells2Ds);
			for(unsigned int i = 0;i<solido1.NumCells0Ds;i++){
				for(unsigned int j = 0;j<solido1.VerticeFaces[i].size();j++){
					solido.Cells2DsVertices(j,i)=solido1.VerticeFaces[i][j];
				}
			}
			for(unsigned int i=0;i<solido.NumCells2Ds;i++){
				unsigned int nlati=0;
				for(unsigned int j=0;j<solido.p;j++){ //prendo gli id di un punto della faccia i
					for(unsigned int k=0;k<solido.p;k++){ //adesso scorro su tutti i segmenti
						for(unsigned int s=0;s<solido.NumCells1Ds;s++){
							if((solido.Cells2DsVertices(j,i)==solido.Cells1DsExtrema(0,s) && solido.Cells2DsVertices(k,i)==solido.Cells1DsExtrema(1,s))){
								solido.Cells2DsEdges(nlati,i)=s;
								nlati++;
							}
						}
					}
				}
			}
			// cout << solido.Cells2DsEdges << endl;
			solido.Cells3DsVertices.resize(solido.NumCells0Ds);
			for(unsigned int i = 0; i<solido.NumCells0Ds;i++)
				solido.Cells3DsVertices[i]=i;
			solido.Cells3DsEdges.resize(solido.NumCells1Ds);
			for(unsigned int i = 0; i<solido.NumCells1Ds;i++)
				solido.Cells3DsEdges[i]=i;
			solido.Cells3DsFaces.resize(solido.NumCells2Ds);
			for(unsigned int i = 0; i<solido.NumCells2Ds;i++)
				solido.Cells3DsFaces[i]=i;
        return 0;
    }
	 
    int CreateMesh(PlatonicSolids& solido){
        // unsigned int npunti = solido.NumCells0Ds + solido.NumCells1Ds * (solido.b-1) + solido.NumCells2Ds * (solido.b-1)*(solido.b-2)/2;
        // unsigned int nlati = solido.NumCells1Ds * solido.b + 3 * solido.NumCells2Ds * (solido.b-1)*(solido.b-2)/2;
        // unsigned int npunti = solido.NumCells0Ds + solido.NumCells1Ds * (solido.b-1) + solido.NumCells2Ds * (solido.b-1)*(solido.b-2)/2;
        // unsigned int nlati1 = solido.NumCells1Ds * solido.b + 3 * solido.NumCells2Ds * (solido.b-1)*(solido.b-2)/2;
        unsigned int nlati = 0;
        unsigned int npunti = 0;
        if(solido.q == 3){
            nlati = 6*solido.b*solido.b;
            npunti = 2*solido.b*solido.b+2;
        }
        if(solido.q == 4){
            nlati = 12*solido.b*solido.b;
            npunti = 4*solido.b*solido.b+2;
        }
        if(solido.q == 5){
            nlati = 30*solido.b*solido.b;
            npunti = 10*solido.b*solido.b+2;
        }
        unsigned int ntriangoli = 1000;
        MatrixXd punti(3,npunti);
        punti.setConstant(2.0);
        MatrixXi lati(2,nlati);
        lati.setConstant(1000);
        MatrixXi latiOriginali(2,solido.NumCells1Ds);
        latiOriginali.setConstant(1000);
        VectorXi verticiOriginali(solido.NumCells0Ds);
        verticiOriginali.setConstant(1000);
        unsigned int counter = 0;
        unsigned int idlato = 0;
        // unsigned int idt = 0; //id triangolo
        // MatrixXi latiTriangoli(3,ntriangoli);
        // MatrixXi verticiTriangoli(3,ntriangoli);
        // cout << "vertici originali 0: " << verticiOriginali << endl;
        for(unsigned int lato = 0;lato<solido.NumCells1Ds;lato++){ //scorro ogni lato
            double x1 = solido.Cells0DsCoordinates(0,solido.Cells1DsExtrema(0,lato));
            double y1 = solido.Cells0DsCoordinates(1,solido.Cells1DsExtrema(0,lato));
            double z1 = solido.Cells0DsCoordinates(2,solido.Cells1DsExtrema(0,lato));
            double x2 = solido.Cells0DsCoordinates(0,solido.Cells1DsExtrema(1,lato));
            double y2 = solido.Cells0DsCoordinates(1,solido.Cells1DsExtrema(1,lato));
            double z2 = solido.Cells0DsCoordinates(2,solido.Cells1DsExtrema(1,lato));
            if(verticiOriginali(solido.Cells1DsExtrema(0,lato))==1000){ // non ho ancora salvato il punto
                verticiOriginali(solido.Cells1DsExtrema(0,lato))=counter;
                punti(0, counter) = x1;
                punti(1, counter) = y1;
                punti(2, counter) = z1;
                counter ++;
                // cout << solido.Cells1DsExtrema(0,lato) << endl;
            }
            // cout << solido.Cells1DsExtrema(0,lato) << " " << verticiOriginali(solido.Cells1DsExtrema(0,lato))<< endl;
            latiOriginali(0,lato)=counter; //primo punto che metterò sul lato
            lati(0,idlato)=verticiOriginali(solido.Cells1DsExtrema(0,lato));
            lati(1,idlato)=counter;
            idlato++;
            for (unsigned int i = 1; i < solido.b; i++) {
                double t = static_cast<double>(i) / (solido.b);
                punti(0, counter) = (1 - t) * x1 + t * x2;
                punti(1, counter) = (1 - t) * y1 + t * y2;
                punti(2, counter) = (1 - t) * z1 + t * z2;
                if(i<solido.b-1){
                    lati(0,idlato) = counter;
                    lati(1,idlato) = counter+1; 
                    idlato++;
                }
                counter++;
            }
            latiOriginali(1,lato)=counter-1;
            if(verticiOriginali(solido.Cells1DsExtrema(1,lato))==1000){ // non ho ancora salvato il punto
                verticiOriginali(solido.Cells1DsExtrema(1,lato))=counter;
                punti(0, counter) = x2;
                punti(1, counter) = y2;
                punti(2, counter) = z2;
                counter++;
                lati(0,idlato)=counter-2;
                lati(1,idlato)=verticiOriginali(solido.Cells1DsExtrema(1,lato));
                idlato++;
                // cout << solido.Cells1DsExtrema(1,lato) << endl;
            }else{
                lati(0,idlato)=counter-1;
                lati(1,idlato)=verticiOriginali(solido.Cells1DsExtrema(1,lato));
                idlato++;
            }
            // cout << solido.Cells1DsExtrema(1,lato) << " " << verticiOriginali(solido.Cells1DsExtrema(1,lato))<< endl;
            
        }
        // cout << "vertici originali 1: " << verticiOriginali << endl;


            //prendere 2 lati. fare divisione per i primi b-2 punti.
        for(unsigned int nfaccia=0;nfaccia<solido.NumCells2Ds;nfaccia++){ //guardo una faccia per volta
            unsigned int lato0 = solido.Cells2DsEdges(0,nfaccia);//salvo i 3 lati della faccia
            unsigned int lato1 = solido.Cells2DsEdges(1,nfaccia);
            unsigned int lato2 = solido.Cells2DsEdges(2,nfaccia);
            unsigned int counterlato0_iniziale = latiOriginali(0,lato0); //salvo counter iniziale e finale per ogni lato (esclusi i vertici) sono dal vertice 0 al vertice 1
            unsigned int counterlato0_finale = latiOriginali(1,lato0);
            unsigned int counterlato1_iniziale = latiOriginali(0,lato1);
            unsigned int counterlato1_finale = latiOriginali(1,lato1);
            unsigned int counterlato2_iniziale = latiOriginali(0,lato2);
            unsigned int counterlato2_finale = latiOriginali(1,lato2);
            if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(0,lato1)){ //estr1 lato0 = estr0 lato1
                for(unsigned int i=0;i<solido.b-2;i++){ //se b=3 scelgo solo il primo punto, numero di segmenti interni da creare e partizionare
                    double x1 = punti(0,counterlato0_iniziale+i);
                    double y1 = punti(1,counterlato0_iniziale+i);
                    double z1 = punti(2,counterlato0_iniziale+i);
                    double x2 = punti(0,counterlato1_finale-i);
                    double y2 = punti(1,counterlato1_finale-i);
                    double z2 = punti(2,counterlato1_finale-i);
                    // cout  << "1 " << x1 << " " << y1 << " " << z1 << endl;
                    // cout  << "2 "  << x2 << " " << y2 << " " << z2 << endl;
                    // cout << counter << endl;
                    for(unsigned int j=1;j<solido.b-(i+1);j++){ //numero di punti per ogni segmento interno
                        double t = static_cast<double>(j) / (solido.b-i-1);
                        punti(0, counter) = (1-t) * x1 + t * x2;
                        punti(1, counter) = (1-t) * y1 + t * y2;
                        punti(2, counter) = (1-t) * z1 + t * z2;
                        if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(0,lato2)){
                            if(j==1){ //giallo
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_iniziale+i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_iniziale+i+1;
                                idlato++;
                            }
                            if(j==1 && j!=solido.b-2-i){//verde
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                            }
                            if(i==0){ //viola
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_iniziale+j-1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_iniziale+j;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i<solido.b-3){ //blu
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j-1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i==solido.b-3){ //rosa
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j-1;
                                idlato++;
                            }
                            if(j<solido.b-2-i&&j>1){//grigio
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                        }
                        if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(1,lato2)){
                            if(j==1){ //giallo
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_iniziale+i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_iniziale+i+1;
                                idlato++;
                            }
                            if(j==1 && j!=solido.b-2-i){//verde
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                            }
                            if(i==0){ //viola
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_finale-j+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_finale-j;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i<solido.b-3){ //blu
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j-1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i==solido.b-3){ //rosa
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j-1;
                                idlato++;
                            }
                            if(j<solido.b-2-i&&j>1){//grigio
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                        }
                        counter++;
                    }
                }
            }
            if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(1,lato1)){ //estr1 lato0 = estr0 lato1
                for(unsigned int i=0;i<solido.b-2;i++){ //se b=3 scelgo solo il primo punto, numero di segmenti interni da creare e partizionare
                    double x1 = punti(0,counterlato0_iniziale+i);
                    double y1 = punti(1,counterlato0_iniziale+i);
                    double z1 = punti(2,counterlato0_iniziale+i);
                    double x2 = punti(0,counterlato1_iniziale+i);
                    double y2 = punti(1,counterlato1_iniziale+i);
                    double z2 = punti(2,counterlato1_iniziale+i);
                    // cout  << "1 " << x1 << " " << y1 << " " << z1 << endl;
                    // cout  << "2 "  << x2 << " " << y2 << " " << z2 << endl;
                    // cout << counter << endl;
                    for(unsigned int j=1;j<solido.b-(i+1);j++){ //numero di punti per ogni segmento interno
                        double t = static_cast<double>(j) / (solido.b-i-1);
                        punti(0, counter) = (1-t) * x1 + t * x2;
                        punti(1, counter) = (1-t) * y1 + t * y2;
                        punti(2, counter) = (1-t) * z1 + t * z2;
                        if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(0,lato2)){
                            if(j==1){ //giallo
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_iniziale+i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_iniziale+i+1;
                                idlato++;
                            }
                            if(j==1 && j!=solido.b-2-i){//verde
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                            }
                            if(i==0){ //viola
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_iniziale+j-1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_iniziale+j;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i<solido.b-3){ //blu
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i==solido.b-3){ //rosa
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j+1;
                                idlato++;
                            }
                            if(j<solido.b-2-i&&j>1){//grigio
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                        }
                        if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(1,lato2)){
                            if(j==1){ //giallo
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_iniziale+i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_iniziale+i+1;
                                idlato++;
                            }
                            if(j==1 && j!=solido.b-2-i){//verde
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                            }
                            if(i==0){ //viola
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_finale-j+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_finale-j;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i<solido.b-3){ //blu
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i==solido.b-3){ //rosa
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j+1;
                                idlato++;
                            }
                            if(j<solido.b-2-i&&j>1){//grigio
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                        }
                        counter++;
                    }
                }
            }
            ////////
            if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(0,lato1)){ //estr0 lato0 = estr0 lato1
                for(unsigned int i=0;i<solido.b-2;i++){ //se b=3 scelgo solo il primo punto, numero di segmenti interni da creare e partizionare
                    double x1 = punti(0,counterlato0_finale-i);
                    double y1 = punti(1,counterlato0_finale-i);
                    double z1 = punti(2,counterlato0_finale-i);
                    double x2 = punti(0,counterlato1_finale-i);
                    double y2 = punti(1,counterlato1_finale-i);
                    double z2 = punti(2,counterlato1_finale-i);
                    // cout  << "1 " << x1 << " " << y1 << " " << z1 << endl;
                    // cout  << "2 "  << x2 << " " << y2 << " " << z2 << endl;
                    // cout << counter << endl;
                    for(unsigned int j=1;j<solido.b-(i+1);j++){ //numero di punti per ogni segmento interno
                        double t = static_cast<double>(j) / (solido.b-i-1);
                        punti(0, counter) = (1-t) * x1 + t * x2;
                        punti(1, counter) = (1-t) * y1 + t * y2;
                        punti(2, counter) = (1-t) * z1 + t * z2;
                        if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(0,lato2)){
                            if(j==1){ //giallo
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_finale-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_finale-i-1;
                                idlato++;
                            }
                            if(j==1 && j!=solido.b-2-i){//verde
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                            }
                            if(i==0){ //viola
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_iniziale+j-1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_iniziale+j;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i<solido.b-3){ //blu
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j-1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i==solido.b-3){ //rosa
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j-1;
                                idlato++;
                            }
                            if(j<solido.b-2-i&&j>1){//grigio
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                        }
                        if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(1,lato2)){
                            if(j==1){ //giallo
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_finale-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_finale-i-1;
                                idlato++;
                            }
                            if(j==1 && j!=solido.b-2-i){//verde
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                            }
                            if(i==0){ //viola
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_finale-j+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_finale-j;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i<solido.b-3){ //blu
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j-1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i==solido.b-3){ //rosa
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_iniziale+j-1;
                                idlato++;
                            }
                            if(j<solido.b-2-i&&j>1){//grigio
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                        }
                        counter++;
                    }
                }
            }
            if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(1,lato1)){ //estr1 lato0 = estr0 lato1
                for(unsigned int i=0;i<solido.b-2;i++){ //se b=3 scelgo solo il primo punto, numero di segmenti interni da creare e partizionare
                    double x1 = punti(0,counterlato0_finale-i);
                    double y1 = punti(1,counterlato0_finale-i);
                    double z1 = punti(2,counterlato0_finale-i);
                    double x2 = punti(0,counterlato1_iniziale+i);
                    double y2 = punti(1,counterlato1_iniziale+i);
                    double z2 = punti(2,counterlato1_iniziale+i);
                    // cout  << "1 " << x1 << " " << y1 << " " << z1 << endl;
                    // cout  << "2 "  << x2 << " " << y2 << " " << z2 << endl;
                    // cout << counter << endl;
                    for(unsigned int j=1;j<solido.b-(i+1);j++){ //numero di punti per ogni segmento interno
                        double t = static_cast<double>(j) / (solido.b-i-1);
                        punti(0, counter) = (1-t) * x1 + t * x2;
                        punti(1, counter) = (1-t) * y1 + t * y2;
                        punti(2, counter) = (1-t) * z1 + t * z2;
                        if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(0,lato2)){
                            if(j==1){ //giallo
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_finale-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_finale-i-1;
                                idlato++;
                            }
                            if(j==1 && j!=solido.b-2-i){//verde
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                            }
                            if(i==0){ //viola
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_iniziale+j-1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_iniziale+j;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i<solido.b-3){ //blu
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i==solido.b-3){ //rosa
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j+1;
                                idlato++;
                            }
                            if(j<solido.b-2-i&&j>1){//grigio
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                        }
                        if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(1,lato2)){
                            if(j==1){ //giallo
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_finale-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato0_finale-i-1;
                                idlato++;
                            }
                            if(j==1 && j!=solido.b-2-i){//verde
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                            }
                            if(i==0){ //viola
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_finale-j+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato2_finale-j;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i<solido.b-3){ //blu
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                            if(j==solido.b-2-i&&i==solido.b-3){ //rosa
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counterlato1_finale-j+1;
                                idlato++;
                            }
                            if(j<solido.b-2-i&&j>1){//grigio
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+1;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-2-i;
                                idlato++;
                                lati(0,idlato)=counter;
                                lati(1,idlato)=counter+solido.b-3-i;
                                idlato++;
                            }
                        }
                        counter++;
                    }
                }
            }
            if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(1,lato1)){
                lati(0,idlato)=counterlato0_iniziale;
                lati(1,idlato)=counterlato1_finale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(1,lato1)){
                lati(0,idlato)=counterlato0_finale;
                lati(1,idlato)=counterlato1_finale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(0,lato1)){
                lati(0,idlato)=counterlato0_iniziale;
                lati(1,idlato)=counterlato1_iniziale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(0,lato1)){
                lati(0,idlato)=counterlato0_finale;
                lati(1,idlato)=counterlato1_iniziale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(1,lato2)){
                lati(0,idlato)=counterlato0_iniziale;
                lati(1,idlato)=counterlato2_finale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(1,lato2)){
                lati(0,idlato)=counterlato0_finale;
                lati(1,idlato)=counterlato2_finale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(0,lato0)==solido.Cells1DsExtrema(0,lato2)){
                lati(0,idlato)=counterlato0_iniziale;
                lati(1,idlato)=counterlato2_iniziale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(1,lato0)==solido.Cells1DsExtrema(0,lato2)){
                lati(0,idlato)=counterlato0_finale;
                lati(1,idlato)=counterlato2_iniziale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(0,lato2)==solido.Cells1DsExtrema(1,lato1)){
                lati(0,idlato)=counterlato2_iniziale;
                lati(1,idlato)=counterlato1_finale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(1,lato2)==solido.Cells1DsExtrema(1,lato1)){
                lati(0,idlato)=counterlato2_finale;
                lati(1,idlato)=counterlato1_finale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(0,lato2)==solido.Cells1DsExtrema(0,lato1)){
                lati(0,idlato)=counterlato2_iniziale;
                lati(1,idlato)=counterlato1_iniziale;
                idlato++;
            }
            if(solido.Cells1DsExtrema(1,lato2)==solido.Cells1DsExtrema(0,lato1)){
                lati(0,idlato)=counterlato2_finale;
                lati(1,idlato)=counterlato1_iniziale;
                idlato++;
            }
        }
        solido.Cells0DsCoordinates.resize(3,counter);
        solido.Cells0DsCoordinates=punti.leftCols(counter);
        solido.Cells1DsExtrema.resize(2,idlato);
        solido.Cells1DsExtrema=lati.leftCols(idlato);

        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double norm = 0.0;
        for (unsigned int i = 0; i < solido.Cells0DsCoordinates.cols(); ++i) {
            x = solido.Cells0DsCoordinates(0, i);
            y = solido.Cells0DsCoordinates(1, i);
            z = solido.Cells0DsCoordinates(2, i);
            norm = sqrt(x*x + y*y + z*z);
            if (norm != 0.0) {
                solido.Cells0DsCoordinates(0, i) = x / norm;
                solido.Cells0DsCoordinates(1, i) = y / norm;
                solido.Cells0DsCoordinates(2, i) = z / norm;
            } else {
                cerr << "Warning: punto nullo alla colonna " << i << endl;
            }
        }


       /*
		MatrixXi triangles_vertices(3, idt); 


		MatrixXi triangles_adjacency = MatrixXi::Zero(idt, idt);

		
		for (unsigned int i = 0; i < idt; i++) {
			
			vector<unsigned int> t1 = {
				static_cast<unsigned int>(triangles_vertices(0, i)),
				static_cast<unsigned int>(triangles_vertices(1, i)),
				static_cast<unsigned int>(triangles_vertices(2, i))
			};

			for (unsigned int j = i + 1; j < idt; j++) {
				
				vector<unsigned int> t2 = {
					static_cast<unsigned int>(triangles_vertices(0, j)),
					static_cast<unsigned int>(triangles_vertices(1, j)),
					static_cast<unsigned int>(triangles_vertices(2, j))
				};

				
				unsigned int common = 0;
				for (unsigned int a : t1)
					for (unsigned int b : t2)
						if (a == b) common++;

				
				if (common == 2) {
					triangles_adjacency(i, j) = 1;
					triangles_adjacency(j, i) = 1; 
				}
			}
		}

		
		cout << "\nAdiacenze triangoli:\n";
		for (unsigned int i = 0; i < idt; i++) {
			cout << "Triangoli adjacency al triangolo " << i << ": ";
			for (unsigned int j = 0; j < idt; j++) {
				if (triangles_adjacency(i, j) == 1)
					cout << j << " ";
			}
			cout << endl;
		}
        
		*/
		// unsigned int n_punti = counter;
		// unsigned int n_lati = idlato;
   
		
		// solido.adjacency.resize(n_punti); //contiene i vertici adiacenti (che comunicano con un lato)
		// for(auto& k : solido.adjacency){
		//     k = {};	
		// }	
		// for(unsigned int i = 0; i < n_lati; i++) {
		// 	unsigned int v0 = lati(0, i);
		// 	unsigned int v1 = lati(1, i);
		// 	solido.adjacency[v0].push_back(v1);
		// 	solido.adjacency[v1].push_back(v0); 
		// }
        // for(auto& vec : solido.adjacency) {
        //     sort(vec.begin(), vec.end(), comp);
        //     vec.erase(unique(vec.begin(), vec.end()), vec.end());
        // }
        VerticeAdjacency(solido);
        // set<pair<unsigned int, unsigned int>> added_edges;
        // for (unsigned int i = 0; i < n_lati; i++) {
        //     unsigned int v0 = lati(0, i);
        //     unsigned int v1 = lati(1, i);
        //     if (v0 > v1) swap(v0, v1);
        //     if (added_edges.insert({v0, v1}).second) {
        //         solido.adjacency[v0].push_back(v1);
        //         solido.adjacency[v1].push_back(v0); 
        //     }
        // }


		// cout << "adjacency" << endl;
        // for(auto& sottovettore: solido.adjacency){
        //         for(auto& el : sottovettore){
        //             cout << el << " ";
        //         }
        //         cout << endl;
        // }   
        // controllo per lati duplicati
        // set<pair<unsigned int, unsigned int>> lato_set;
        // for(unsigned int i = 0; i < n_lati; i++) {
        //     unsigned int u = min(lati(0, i), lati(1, i));
        //     unsigned int v = max(lati(0, i), lati(1, i));
        //     if (!lato_set.insert({u, v}).second) {
        //         cout << "Lato duplicato trovato: " << u << " - " << v << endl;
        //     }
        // }
        
        set<tuple<unsigned int, unsigned int, unsigned int>> triangoli;

        for(unsigned int u = 0; u < solido.adjacency.size(); ++u) {
            for(unsigned int v : solido.adjacency[u]) {
                if(v <= u) 
                    continue; // evita duplicati simmetrici

                for(unsigned int w : solido.adjacency[v]) {
                    if(w <= v || w == u)
                        continue; // evita duplicati e cicli
                    // controlla che w sia adiacente anche a u
                    if(find(solido.adjacency[u].begin(), solido.adjacency[u].end(), w) != solido.adjacency[u].end()) {
                        bool has_uv = false, has_vw = false, has_wu = false;
                        for (int i = 0; i < solido.Cells1DsExtrema.cols(); ++i) {
                            unsigned int a = solido.Cells1DsExtrema(0, i);
                            unsigned int b = solido.Cells1DsExtrema(1, i);
                            if ((a == u && b == v) || (a == v && b == u)) has_uv = true;
                            if ((a == v && b == w) || (a == w && b == v)) has_vw = true;
                            if ((a == w && b == u) || (a == u && b == w)) has_wu = true;
                        }
                        if (has_uv && has_vw && has_wu && u!=v && v!=w && u!=w) {
                            vector<unsigned int> vert = {u, v, w};
                            sort(vert.begin(), vert.end(), comp);
                            triangoli.insert({vert[0], vert[1], vert[2]});
                        }
                    }
                }
            }
        }
        // Stampa per verificare triangoli duplicati (oltre a quelli già ordinati)
        cout << "Numero triangoli trovati: " << triangoli.size() << endl;
        // // Rimuove triangoli duplicati con stessi vertici in ordine diverso
        // set<tuple<unsigned int, unsigned int, unsigned int>> triangoli_unici;
        // for (const auto& t : triangoli) {
        //     vector<unsigned int> vert = {get<0>(t), get<1>(t), get<2>(t)};
        //     sort(vert.begin(), vert.end());
        //     triangoli_unici.insert({vert[0], vert[1], vert[2]});
        // }
        // triangoli = triangoli_unici;
        // stampo triangoli per controllo
        // for(const auto& t : triangoli) {
        //     unsigned int a = get<0>(t);
        //     unsigned int b = get<1>(t);
        //     unsigned int c = get<2>(t);
        //     cout << "Triangolo: " << a << " " << b << " " << c << endl;
        // }
        MatrixXi verticiFaccia = MatrixXi::Zero(3,0);
        MatrixXi latiFaccia = MatrixXi::Zero(3,0);
        // solido.Cells2DsVertices.resize(3,10); //inizializzare matrice a 0
        // solido.Cells2DsEdges.resize(3,10);
        for(const auto& [v1, v2, v3] : triangoli) {
            Vector3i vertici(v1, v2, v3);
            Vector3i lati = Vector3i::Zero();
            bool triangolo_valido = true;
            for(unsigned int i = 0; i < 3; ++i) {
                unsigned int a = vertici(i);
                unsigned int b = vertici((i + 1) % 3);
                bool edge_found = false;
                // Controlla che il lato (a, b) esista in solido.Cells1DsExtrema
                for(unsigned int idlato = 0; idlato < solido.Cells1DsExtrema.cols(); ++idlato) {
                    unsigned int u = solido.Cells1DsExtrema(0, idlato);
                    unsigned int v = solido.Cells1DsExtrema(1, idlato);
                    if((a == u && b == v) || (a == v && b == u)) {
                        lati(i) = idlato;
                        edge_found = true;
                        break;
                    }
                }
                if (!edge_found) {
                    lati(i) = std::numeric_limits<unsigned int>::max();
                    triangolo_valido = false;
                    std::cerr << "Errore: lato non trovato per triangolo (" << a << ", " << b << ")" << std::endl;
                }
            }
            // Inserisci solo se tutti i lati sono validi
            if (triangolo_valido &&
                lati(0) != std::numeric_limits<unsigned int>::max() &&
                lati(1) != std::numeric_limits<unsigned int>::max() &&
                lati(2) != std::numeric_limits<unsigned int>::max()) {
                // solido.Cells2DsVertices.conservativeResize(3, solido.Cells2DsVertices.cols() + 1);
                // solido.Cells2DsEdges.conservativeResize(3, solido.Cells2DsEdges.cols() + 1);
                // solido.Cells2DsVertices.col(solido.Cells2DsVertices.cols() - 1) = vertici;
                // solido.Cells2DsEdges.col(solido.Cells2DsEdges.cols() - 1) = lati;
                verticiFaccia.conservativeResize(3, verticiFaccia.cols() + 1);
                latiFaccia.conservativeResize(3, latiFaccia.cols() + 1);
                verticiFaccia.col(verticiFaccia.cols()-1) = vertici;
                latiFaccia.col(latiFaccia.cols()-1) = lati;
            }
            // controllo area triangoli
            // double x1 = solido.Cells0DsCoordinates(0,vertici[0]);
            // double y1 = solido.Cells0DsCoordinates(1,vertici[0]);
            // double z1 = solido.Cells0DsCoordinates(2,vertici[0]);
            // double x2 = solido.Cells0DsCoordinates(0,vertici[1]);
            // double y2 = solido.Cells0DsCoordinates(1,vertici[1]);
            // double z2 = solido.Cells0DsCoordinates(2,vertici[1]);
            // double x3 = solido.Cells0DsCoordinates(0,vertici[2]);
            // double y3 = solido.Cells0DsCoordinates(1,vertici[2]);
            // double z3 = solido.Cells0DsCoordinates(2,vertici[2]);
            // Vector3d A(x1, y1, z1);
            // Vector3d B(x2, y2, z2);
            // Vector3d C(x3, y3, z3);
            // double area = 0.5 * (B - A).cross(C - A).norm();
            // cout << area << endl;
        }
        solido.Cells2DsEdges = latiFaccia;
        solido.Cells2DsVertices = verticiFaccia;
        solido.NumCells2Ds = solido.Cells2DsEdges.cols();
        // cout << solido.Cells2DsVertices.transpose() << endl;
        // cout << endl;
        solido.Cells2DsNeighborhood.resize(solido.NumCells2Ds); //facce adiacenti tra loro (con un lato in comune)
        for(unsigned int i=0;i<solido.NumCells2Ds;i++){
            set<unsigned int> vicini;
            for(unsigned int latoi : solido.Cells2DsEdges.col(i)){
                if (latoi == std::numeric_limits<unsigned int>::max()) continue;
                for(unsigned int j=0;j<solido.NumCells2Ds;j++){
                    if(i==j) continue;
                    for(unsigned int latoj : solido.Cells2DsEdges.col(j)){
                        if (latoj == std::numeric_limits<unsigned int>::max()) continue;
                        if(i < j && latoi == latoj && vicini.size()<4){
                            vicini.insert(j);
                        }
                    }
                }
            }
            solido.Cells2DsNeighborhood[i] = vector<unsigned int>(vicini.begin(),vicini.end());
        }
            solido.NumCells0Ds = solido.Cells0DsCoordinates.cols();
            solido.NumCells1Ds = solido.Cells1DsExtrema.cols();
            solido.VerticeFaces.resize(solido.NumCells0Ds);
            // //riservo memoria
            for (auto& sottovettore : solido.VerticeFaces) {
                sottovettore.reserve(10);
            }
            
            // for(unsigned int i = 0;i<solido.Cells2DsVertices.cols();i++){
            //     for(unsigned int j = 0;j<solido.Cells2DsVertices.col(i).size();j++){
            //         // cout << "j " << j << " i " << i << endl;
            //         solido.VerticeFaces[solido.Cells2DsVertices(j,i)].push_back(i);
            //         // cout << solido.Cells2DsVertices(j,i) << endl;
            //     }
            // }
            // for(auto& el : solido.VerticeFaces){
            //     sort(el.begin(),el.end(),comp);
            // }
            for (unsigned int i = 0; i < solido.Cells2DsVertices.cols(); ++i) {
                set<unsigned int> inseriti;
                for (unsigned int j = 0; j < solido.Cells2DsVertices.rows(); ++j) {
                    unsigned int v = solido.Cells2DsVertices(j, i);
                    if (inseriti.insert(v).second) {
                        solido.VerticeFaces[v].push_back(i);
                    }
                }
            }

            // Controllo: verifica correttezza di VerticeFaces
            // cout << "\nVerticeFaces:" << endl;
            // for(unsigned int i = 0; i < solido.VerticeFaces.size(); ++i){
            //     cout << "Vertice " << i << ": ";
            //     for(unsigned int f : solido.VerticeFaces[i]){
            //         cout << f << " ";
            //     }
            //     cout << endl;
            // }
            // for(auto& sottovettore: solido.VerticeFaces){
            //     for(auto& el : sottovettore){
            //         cout << el << " ";
            //     }
            //     cout << endl;
            // }

		/*
		solido.Cells3DsVertices.resize(solido.NumCells0Ds);
        for(unsigned int i = 0; i<solido.NumCells0Ds;i++)
            solido.Cells3DsVertices[i]=i;
        solido.Cells3DsEdges.resize(solido.NumCells1Ds);
        for(unsigned int i = 0; i<solido.NumCells1Ds;i++)
            solido.Cells3DsEdges[i]=i;
        solido.Cells3DsFaces.resize(solido.NumCells2Ds);
        for(unsigned int i = 0; i<solido.NumCells2Ds;i++)
            solido.Cells3DsFaces[i]=i;
		*/
        return 0;
    }
	
    int VerticeAdjacency(PlatonicSolids& solido){
        solido.adjacency.resize(solido.Cells0DsCoordinates.cols()); //contiene i vertici adiacenti (che comunicano con un lato)
		for(auto& k : solido.adjacency){
		    k = {};	
		}	

        set<pair<unsigned int, unsigned int>> added_edges;
        for (unsigned int i = 0; i < solido.Cells0DsCoordinates.cols(); i++) {
            unsigned int v0 = solido.Cells1DsExtrema(0, i);
            unsigned int v1 = solido.Cells1DsExtrema(1, i);
            if (v0 > v1) swap(v0, v1);
            if (added_edges.insert({v0, v1}).second) {
                solido.adjacency[v0].push_back(v1);
                solido.adjacency[v1].push_back(v0); 
            }
        }
        return 0;
    }
	
    int ShortestPath(PlatonicSolids& solido){
	    int n = solido.adjacency.size(); //numero dei nodi
        // for(unsigned int i=0; i < n; i++){
		//     for (unsigned int j = 0; j < n; j++){
        //         if (solido.id_vertice1 == i && solido.id_vertice2 == j) {
        //             cout << "trovato" <<  endl; 
        //         }
        //     }
        // }
		
		if(solido.id_vertice1 < n  && solido.id_vertice2 < n){
			queue<int> Q;
			vector<int> distanza(n,-1);
			vector<bool> visited(n, false);
            vector<int> predecessore(n, -1);
			
			Q.push(solido.id_vertice1);
			visited[solido.id_vertice1]=true;
			distanza[solido.id_vertice1]=0;

			while(!Q.empty()){
				unsigned int u=Q.front(); //o Q.front()
				Q.pop();

				for(unsigned int w : solido.adjacency[u]){
				    if(!visited[w]){
					    visited[w] = true;
					    distanza[w] = distanza[u] + 1;
                        predecessore[w] = u;
					    Q.push(w);
					    if(w == solido.id_vertice2){
                            cout << "Cammino trovato: " << distanza[w] << " passi." << endl;
                            vector<int> path;
                            vector<int> edgepath;
                            for(int at = w; at != -1; at = predecessore[at]) {
                                path.push_back(at);
                            }
                            reverse(path.begin(), path.end());
                            cout << "Cammino vertici: ";
                            for(unsigned int v : path){
                                cout << v << " ";
                                const auto it = solido.ShortPathVertices.find(1);
                                if(it == solido.ShortPathVertices.end())
                                {
                                    solido.ShortPathVertices.insert({1, {v}});
                                }
                                else
                                {
                                    it->second.push_back(v);
                                }
                            }
                            cout << endl;
                            for(unsigned int i = 0;i<solido.NumCells1Ds;i++){
                                for(unsigned int k = 0;k<path.size();k++){
                                    if((solido.Cells1DsExtrema(0,i)==path[k] && (solido.Cells1DsExtrema(1,i)==path[(k+1)%(path.size())])) || (solido.Cells1DsExtrema(1,i)==path[k] && (solido.Cells1DsExtrema(0,i)==path[(k+1)%(path.size())]))){
                                        edgepath.push_back(i);
                                    }
                                }
                            }
                            cout << "Cammino lati: ";
                            for(unsigned int e : edgepath){
                                cout << e << " ";
                                const auto it = solido.ShortPathEdges.find(1);
                                if(it == solido.ShortPathEdges.end())
                                {
                                    solido.ShortPathEdges.insert({1, {e}});
                                }
                                else
                                {
                                    it->second.push_back(e);
                                }
                            }
                            cout << endl;
                            return 0;
                        }
				    }
			    }   
		    }
        }else{
            cerr << "Not valid vertice indexes" << endl;
            return 1;
        }
      return 0;
	}
}
	

