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
                 solido.b = stoi(argv[3]);
                 solido.c = stoi(argv[4]);
                 solido.id_vertice1 = stoi(argv[5]);
		 		 solido.id_vertice2 = stoi(argv[6]);
                } if (stoi(argv[3]) == 0 || stoi(argv[4]) == 0 || stoi(argv[5]) == 0 || stoi(argv[6])== 0){
		        	cerr << "Not Valid input" << endl; //every other case
		   		return 1;
                }
        //check for p and q
        if(solido.p < 3 || solido.q <3){
            cerr << "Not Valid input" << endl;
            return 1;
        }
        return 0;
    };
   
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

            //creare un vettore di vettori, ogni vettore piccolo contiene gli id delle facce che sono solido.adjacency alla faccia, ovvero quelli che hanno un lato in comune
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

            solido.VerticeFaces.reserve(solido.NumCells0Ds);
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

            solido.VerticeFaces.reserve(solido.NumCells0Ds); //un vettore per ogni vertice
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

            solido.VerticeFaces.reserve(solido.NumCells0Ds); //un vettore per ogni vertice
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
		}	

    return 0;
    };
	
    void FileCell0Ds(const PlatonicSolids& solido) {
		
		ofstream outputFile;
		outputFile.open("Cells0Ds.txt");
		if (!outputFile) {
			cerr << "Error opening file: Cells0Ds.txt" << endl;
		return;
		}

		outputFile << "ID    x              y              z" << endl;
		for (unsigned int i = 0; i < solido.NumCells0Ds; ++i) {
			outputFile << i << "    "
				   << fixed << setprecision(8)
				   << solido.Cells0DsCoordinates(0, i) << "    "
				   << solido.Cells0DsCoordinates(1, i) << "    "
				   << solido.Cells0DsCoordinates(2, i)
				   << endl;
		}

		outputFile.close();
		cout << "Cells0Ds.txt file created successfully with " << solido.NumCells0Ds << " vertices." << endl;
    }
	
	void FileCell1Ds(const PlatonicSolids& solido) {
		
		ofstream outputFile;
		outputFile.open("Cells1Ds.txt");
		if (!outputFile) {
			cerr << "Error opening file: Cells1Ds.txt" << endl;
        return;
		}

    
		outputFile << "ID\tInitial Vertices\tFinal Vertices" << endl;

    
		for (unsigned int i = 0; i < solido.NumCells1Ds; ++i) {
			outputFile << i << "\t"
                   << solido.Cells1DsExtrema(0, i) << "\t\t"
                   << solido.Cells1DsExtrema(1, i) << endl;
		}

		outputFile.close();
		cout << "Cells1Ds.txt file created successfully with " << solido.NumCells1Ds << " edges." << endl;
	}
	
	void FileCell2Ds(const PlatonicSolids& solido) {
		
		ofstream outputFile;
		outputFile.open("Cells2Ds.txt");
		if (!outputFile) {
			cerr << "Error opening file: Cells2Ds.txt" << endl;
        return;
		}

    
		outputFile << "ID\tNumVertices\tNumEdges\tVerticesIDs\t\tEdgesIDs" << endl;

		for (unsigned int i = 0; i < solido.NumCells2Ds; ++i) {
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
        cout << "Cells2Ds.txt file created successfully with " << solido.NumCells2Ds << " faces." << endl;
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
        solido.Cells0DsCoordinates = MatrixXd::Zero(3,solido.NumCells0Ds);
        solido.Cells1DsExtrema = MatrixXi::Zero(2,solido.NumCells1Ds);
        solido.Cells2DsVertices = MatrixXi::Zero(3, solido.NumCells2Ds);
        solido.Cells2DsEdges = MatrixXi::Zero(3, solido.NumCells2Ds);
        solido.Cells2DsNumEdges = VectorXi::Zero(solido.NumCells2Ds);
        for(unsigned int j=0; j<solido1.NumCells2Ds;j++){
            int Pid = solido1.Cells2DsVertices(0,j);
            int Qid = solido1.Cells2DsVertices(1,j);
            int Sid = solido1.Cells2DsVertices(2,j);
            double Px = solido1.Cells0DsCoordinates(0,Pid);
            double Py = solido1.Cells0DsCoordinates(1,Pid);
            double Pz = solido1.Cells0DsCoordinates(2,Pid);
            double Qx = solido1.Cells0DsCoordinates(0,Qid);
            double Qy = solido1.Cells0DsCoordinates(1,Qid);
            double Qz = solido1.Cells0DsCoordinates(2,Qid);
            double Sx = solido1.Cells0DsCoordinates(0,Sid);
            double Sy = solido1.Cells0DsCoordinates(1,Sid);
            double Sz = solido1.Cells0DsCoordinates(2,Sid);
            solido.Cells0DsCoordinates(0,j) = (Px+Qx+Sx)/3;
            solido.Cells0DsCoordinates(1,j)= (Py+Qy+Sy)/3;
            solido.Cells0DsCoordinates(2,j)= (Pz+Qz+Sz)/3;
            double x = solido.Cells0DsCoordinates(0,j);
            double y = solido.Cells0DsCoordinates(1,j);
            double z = solido.Cells0DsCoordinates(2,j);
            double norm = sqrt(x*x+y*y+z*z);
            solido.Cells0DsCoordinates(0,j) /= norm;
            solido.Cells0DsCoordinates(1,j) /= norm;
            solido.Cells0DsCoordinates(2,j) /= norm;
        } 
        
        //accedo a Cells2dNeighborhood e collego ciascun vertice di solido ai vertici corrispondenti alle facce solido.adjacency a quello (inizio ad avere il duplicato dei lati)
        solido.Cells1DsExtrema = MatrixXi::Zero(2,solido.NumCells1Ds);
        solido.Cells1DsExtrema.setConstant(1000);
        unsigned int nlati = 0;
        for(unsigned int idfaccia = 0;idfaccia<solido1.NumCells2Ds;idfaccia++){
            vector<unsigned int> vettore = solido1.Cells2DsNeighborhood[idfaccia];
            for(unsigned int adj = 0;adj<solido.q;adj++){
                if(idfaccia<vettore[adj]){
                    solido.Cells1DsExtrema(0,nlati) = idfaccia;
                    solido.Cells1DsExtrema(1,nlati) = vettore[adj];
                    nlati ++;
                }
            }
        }
        solido.Cells2DsVertices = MatrixXi::Zero(solido.p, solido.NumCells2Ds);
        solido.Cells2DsEdges = MatrixXi::Zero(solido.p, solido.NumCells2Ds);
        for(unsigned int i = 0;i<solido1.NumCells0Ds;i++){
            for(unsigned int j = 0;j<solido1.q;j++){
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

        return 0;
    }
    // int SphereProjection(double& x, double& y, double& z){
    //     double norm = sqrt(x*x+y*y+z*z);
    //     x /= norm;
    //     y /= norm;
    //     z /= norm;
    //     return 0;
    // }
    int CreateMesh(PlatonicSolids& solido){
        cout << "create mesh" << endl;
        unsigned int npunti = 10000;//3*solido.b*solido.b*5;
        unsigned int nlati = 10000;//solido.b*solido.b*solido.NumCells2Ds*3/2+solido.NumCells1Ds*5;
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
        unsigned int idt = 0; //id triangolo
        MatrixXi latiTriangoli(3,ntriangoli);
        MatrixXi verticiTriangoli(3,ntriangoli);
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
                lati(0,idlato)=counterlato0_finale;
                lati(1,idlato)=counterlato1_iniziale;
                idlato++;
            }
        }



        solido.Cells0DsCoordinates.resize(3,counter);
        solido.Cells0DsCoordinates=punti.leftCols(counter);
        solido.Cells1DsExtrema.resize(2,idlato);
        solido.Cells1DsExtrema=lati.leftCols(idlato);

        // for(unsigned int j=0;j<counter;j++){
        //     double x = solido.Cells0DsCoordinates(0,j);
        //     double y = solido.Cells0DsCoordinates(1,j);
        //     double z = solido.Cells0DsCoordinates(2,j);
        //     double norm = sqrt(x*x+y*y+z*z);
        //     solido.Cells0DsCoordinates(0,j) /= norm;
        //     solido.Cells0DsCoordinates(1,j) /= norm;
        //     solido.Cells0DsCoordinates(2,j) /= norm;
        // }
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
		unsigned int n_punti = counter-1;
		unsigned int n_lati = idlato-1;
   
		
		solido.adjacency.resize(n_punti);
		for(auto& k : solido.adjacency){
		    k = {};	
		}	
		for (unsigned int i = 0; i < n_lati; i++) {
			unsigned int v0 = lati(0, i);
			unsigned int v1 = lati(1, i);
			solido.adjacency[v0].push_back(v1);
			solido.adjacency[v1].push_back(v0); 
		}
		
	
		
		cout << "\nNumero punti generati: " << n_punti << endl;
		for (unsigned int i = 0; i < n_punti; i++) {
			cout << "Punto " << i << ": (" << punti(0,i) << ", " << punti(1,i) << ", " << punti(2,i) << ")" << endl;
		}

		cout << "\nNumero lati generati: " << n_lati << endl;
		for (unsigned int i = 0; i < n_lati; i++) {
			cout << "Lato " << i << ": " << lati(0,i) << " <-> " << lati(1,i) << endl;
		}

     /*
		cout << "\nAdiacenza vertici:\n";
		for (unsigned int i = 0; i < n_punti; i++) {
			cout << "Vertice " << i << ": ";
			for (auto v : solido.adjacency[i])
				cout << v << " ";
			cout << endl;
		}
		*/
		
        return 0;
    }


 
	int ShortestPath(PlatonicSolids& solido){
		
		cout << solido.Cells0DsCoordinates.size() << endl;
		
	    int n= solido.Cells1DsExtrema.size(); //numero dei nodi
		/*
        for (unsigned int i=0; i < n; i++) {
		   for (unsigned int j =0; j < n; j++){
             if (solido.id_vertice1 == i && solido.id_vertice2 == j) {
                 cout << "trovato" <<  endl; 
             }
			 //else {
              //  cerr << "Error" << endl;
             //}
           }
        }
		*/
		MatrixXi AdjacencyVertices = MatrixXi::Zero(n,n);
		
		if(solido.id_vertice1 < n  && solido.id_vertice2 < n){
		   
			queue<int> Q;
			vector<int> distanza(n,-1);
			vector<bool> visited(n, false);
			
			Q.push(solido.id_vertice1);
			visited[solido.id_vertice1]=true;
			distanza[solido.id_vertice1]=0;

			while(!Q.empty()){
				int u=Q.back(); //o Q.front()
				Q.pop();
				for(int w : solido.adjacency[u]){
				  if(!visited[w]){
					visited[w] = true;
					distanza[w] = distanza[u] + 1;
					Q.push(w);
					 if(w == solido.id_vertice2) {
						cout << "Cammino minimo trovato: " << distanza[w] << " passi." << endl;
						return distanza[w];
					}

				  }
			    }   
			}
			
			
			
		 cout << "Nessun cammino tra " << solido.id_vertice1 << " e " << solido.id_vertice2 << endl;
	    } 
	 
      return 0;
	}	
 
}	
	
