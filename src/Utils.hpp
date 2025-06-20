#pragma once
#include"PlatonicSolids.hpp"
#include<iostream>
#include<string>

using namespace std;

namespace PlatonicLibrary{
    int ImportValue(int argc, char  *argv[],PlatonicSolids& solido);
    int CreateSolid(PlatonicSolids& solido);
	void FileCell0Ds(const PlatonicSolids& solido);
	void FileCell1Ds(const PlatonicSolids& solido);
	void FileCell2Ds(const PlatonicSolids& solido);
    int DualPolyhedron(PlatonicSolids& solido,PlatonicSolids& solido1);
    int CreateMesh(PlatonicSolids& solido);
    int SphereProjection(double& x, double& y, double& z);
}

