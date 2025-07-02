#pragma once

#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PlatonicLibrary{
    struct PlatonicSolids
    {
        //define numbers of cells
        unsigned int NumCells0Ds = 0;
        unsigned int NumCells1Ds = 0;
        unsigned int NumCells2Ds = 0;
        unsigned int NumCells3Ds = 0;

        //define ids
        vector<unsigned int> Cells0DsId = {};
        vector<unsigned int> Cells1DsId = {};
        vector<unsigned int> Cells2DsId = {};
        vector<unsigned int> Cells3DsId = {};

        //points
        MatrixXd Cells0DsCoordinates;
        //segments extrema
        MatrixXi Cells1DsExtrema = {};
		MatrixXi Cell1DVertices={};
        //polygons
        //vertices

        MatrixXi Cells2DsVertices = {};
		vector<vector<unsigned int>> Cells1DsVertices = {};

        vector<vector<unsigned int>> VerticeFaces = {};
        //vector<vector<unsigned int>> Cells2DsVertices = {{}};

        //edges
        MatrixXi Cells2DsEdges = {};
        //edges number
        VectorXi Cells2DsNumEdges = {};
        
        vector<vector<unsigned int>> Cells2DsNeighborhood;

        VectorXi Cells3DsVertices = {};
        VectorXi Cells3DsEdges = {};
        VectorXi Cells3DsFaces = {};
        //parameters
        unsigned int p = 0;
        unsigned int q = 0;
        unsigned int b = 0;
        unsigned int c = 0;
		unsigned int id_vertice1 = 0;
		unsigned int id_vertice2 = 0;

        vector<vector<unsigned int>> adjacency = {}; 

        map<unsigned int, list<unsigned int>> ShortPathVertices = {};
        map<unsigned int, list<unsigned int>> ShortPathEdges = {};

    };
}