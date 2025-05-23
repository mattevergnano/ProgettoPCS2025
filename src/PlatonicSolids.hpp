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
        //polygons
        //vertices
        vector<vector<unsigned int>> Cells2DsVertices = {};
		vector<vector<unsigned int, unsigned int>> Cells1DsVertices={};
        //edges
        vector<vector<unsigned int>> Cells2DsEdges = {};
        //edges number
        vector<unsigned int> Cells2DsNumEdges = {};

        //parameters
        unsigned int p = 0;
        unsigned int q = 0;
        unsigned int b = 0;
        unsigned int c = 0;

    };
}