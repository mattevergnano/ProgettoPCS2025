#include<iostream>
#include "PlatonicSolids.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"

using namespace std;
using namespace PlatonicLibrary;

int main(int argc,char  *argv[]){
    cout << "ProgettoPCS2025" << endl;
    PlatonicSolids solido;
    ImportValue(argc,argv,solido);
    CreateSolid(solido);
    Gedim::UCDUtilities utilities;
    {
        vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

        cell0Ds_properties[0].Label = "Marker";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(solido.NumCells0Ds, 0.0);
        cell0Ds_properties[0].Data = cell0Ds_marker.data();
        utilities.ExportPoints("./Cell0Ds.inp",
                               solido.Cells0DsCoordinates
                            );
    }

    {

        vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].Label = "Marker";
        cell1Ds_properties[0].UnitLabel = "-";
        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(solido.NumCells1Ds, 0.0);

        cell1Ds_properties[0].Data = cell1Ds_marker.data();

        utilities.ExportSegments("./Cell1Ds.inp",
                                 solido.Cells0DsCoordinates,
                                 solido.Cells1DsExtrema,
                                 {});
    }
    

    return 0;   
}