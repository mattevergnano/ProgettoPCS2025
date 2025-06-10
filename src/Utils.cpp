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
        } else{
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
        
        //accedo a Cells2dNeighborhood e collego ciascun vertice di solido ai vertici corrispondenti alle facce adiacenti a quello (inizio ad avere il duplicato dei lati)
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
    int CreateMesh(PlatonicSolids& solido){
        cout << "create mesh" << endl;
        unsigned int npunti = solido.NumCells2Ds * (solido.b+1)*(solido.b+2)/2;
        unsigned int nlati = solido.b*solido.b*solido.NumCells2Ds*3/2+solido.NumCells1Ds;
        MatrixXd punti = MatrixXd::Zero(3,npunti);
        MatrixXi lati = MatrixXi::Zero(2,nlati);
        unsigned int counter = 0;
        for(unsigned int nfaccia=0;nfaccia<solido.NumCells2Ds;nfaccia++){
            // cout << solido.Cells0DsCoordinates(0,v1) << solido.Cells0DsCoordinates(1,v1) << solido.Cells0DsCoordinates(2,v1) << endl;
            //divido ogni lato in b segmenti
            for(unsigned int lato = 0;lato<3;lato++){
                double x1 = solido.Cells0DsCoordinates(0,solido.Cells1DsExtrema(0,solido.Cells2DsEdges(lato,nfaccia)));
                double y1 = solido.Cells0DsCoordinates(1,solido.Cells1DsExtrema(0,solido.Cells2DsEdges(lato,nfaccia)));
                double z1 = solido.Cells0DsCoordinates(2,solido.Cells1DsExtrema(0,solido.Cells2DsEdges(lato,nfaccia)));
                double x2 = solido.Cells0DsCoordinates(0,solido.Cells1DsExtrema(1,solido.Cells2DsEdges(lato,nfaccia)));
                double y2 = solido.Cells0DsCoordinates(1,solido.Cells1DsExtrema(1,solido.Cells2DsEdges(lato,nfaccia)));
                double z2 = solido.Cells0DsCoordinates(2,solido.Cells1DsExtrema(1,solido.Cells2DsEdges(lato,nfaccia)));
                for (unsigned int i = 0; i < solido.b; i++) {
                    double t = static_cast<double>(i) / (solido.b);
                    punti(0, counter) = (1 - t) * x1 + t * x2;
                    punti(1, counter) = (1 - t) * y1 + t * y2;
                    punti(2, counter) = (1 - t) * z1 + t * z2;
                    counter++;
                }
            }
        }
        cout << counter << endl;
        cout << npunti << endl;
        solido.Cells0DsCoordinates.resize(3,npunti);
        solido.Cells0DsCoordinates=punti;
        return 0;
    }
}