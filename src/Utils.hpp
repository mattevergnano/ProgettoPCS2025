#pragma once
#include"PlatonicSolids.hpp"
#include<iostream>
#include<string>

using namespace std;

namespace PlatonicLibrary{
    int ImportValue(int argc, char  *argv[],PlatonicSolids& solido);
    int CreateSolid(PlatonicSolids& solido);
	void FileCell0Ds(PlatonicSolids& solido);
	void FileCell1Ds(PlatonicSolids& solido);
	void FileCell2Ds(PlatonicSolids& solido);
	void FileCell3Ds(PlatonicSolids& solido);
    int DualPolyhedron(PlatonicSolids& solido,PlatonicSolids& solido1);
    int CreateMesh(PlatonicSolids& solido);
	int ShortestPath(PlatonicSolids& solido);
}

