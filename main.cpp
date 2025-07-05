#include<iostream>
#include "PlatonicSolids.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"

using namespace std;
using namespace PlatonicLibrary;

int main(int argc,char  *argv[]){
    PlatonicSolids solido;
    ImportValue(argc,argv,solido);
    if(argc==3){
        CreateSolid(solido);
    }else{
	if(solido.p==3){
        CreateSolid(solido);
        CreateMesh(solido);
    }
    if(solido.q==3 && solido.p != 3){
        PlatonicSolids solido1;
        solido1.p=solido.q;
        solido1.q=solido.p;
        solido1.b=solido.b;
        solido1.id_vertice1 = solido.id_vertice1;
        solido1.id_vertice2 = solido.id_vertice2;
        CreateSolid(solido1);
        CreateMesh(solido1);
        DualPolyhedron(solido,solido1);
        VerticeAdjacency(solido);
    }}
    FileCell0Ds(solido);
	FileCell1Ds(solido);
	FileCell2Ds(solido);
	FileCell3Ds(solido);
    ShortestPath(solido);
    
    Gedim::UCDUtilities utilities;
    {
        vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

        cell0Ds_properties[0].Label = "Marker";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(solido.NumCells0Ds, 0.0);

        for(const auto &m : solido.ShortPathVertices)
            for(const unsigned int id: m.second)
                cell0Ds_marker.at(id) = m.first;

        cell0Ds_properties[0].Data = cell0Ds_marker.data();
        utilities.ExportPoints("./Cell0Ds.inp",
                               solido.Cells0DsCoordinates,
                               cell0Ds_properties
                            );
    }

    {

        vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].Label = "Marker";
        cell1Ds_properties[0].UnitLabel = "-";
        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(solido.NumCells1Ds, 0.0);

        for(const auto &m : solido.ShortPathEdges)
            for(const unsigned int id: m.second)
                cell1Ds_marker.at(id) = m.first;

        cell1Ds_properties[0].Data = cell1Ds_marker.data();

        utilities.ExportSegments("./Cell1Ds.inp",
                                 solido.Cells0DsCoordinates,
                                 solido.Cells1DsExtrema,
                                 {},
                                cell1Ds_properties);
    }
    return 0;   
}