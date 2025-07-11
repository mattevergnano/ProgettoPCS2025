#pragma once

#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include "Utils.hpp"

using namespace std;

namespace PlatonicLibrary {

TEST(ImportValueTest, CasePQB0) {
    
	// input like p q b 0
	char str[] = "";
    char p[] = "8";
    char q[] = "6";
    char b[] = "4";
    char NullValue[] = "0";

    char* argv[] = { str, p, q, b, NullValue };
    int argc = 5;

    PlatonicSolids solido;

    int result = ImportValue(argc, argv, solido);

    EXPECT_EQ(result, 0);
    EXPECT_EQ(solido.p, 8);
    EXPECT_EQ(solido.q, 6);
    EXPECT_EQ(solido.b, 4);
    EXPECT_EQ(solido.c, 0);
	
}	
	
TEST(ImportValueTest, CasePQ0B) {
    
	// input like p q 0 b
	char str[] = "";
    char p[] = "7";
    char q[] = "5";
    char NullValue[] = "0";
    char b[] = "4";
    char* argv[] = { str, p, q, NullValue, b };
    int argc = 5;

    PlatonicSolids solido;
    int result = ImportValue(argc, argv, solido);

    EXPECT_EQ(result, 0);
    EXPECT_EQ(solido.p, 7);
    EXPECT_EQ(solido.q, 5);
    EXPECT_EQ(solido.c, 0);
    EXPECT_EQ(solido.b, 4);
	
}		
	
TEST(ImportValueTest, CasePQ00) {
    
	// input like p q 0 0 
    char str[] = "";
    char p[] = "6";
    char q[] = "4";
    char NullValue[] = "0";
    char NullValue1[] = "0";

    char* argv[] = { str, p, q, NullValue, NullValue1 };
    int argc = 5;

    PlatonicSolids solido;
    int result = ImportValue(argc, argv, solido);

    EXPECT_EQ(result, 1); 
}
	
TEST(ImportValueTest, CasePQBC) {
   
   // input like p q b c
    char str[] = "";
    char p[] = "6";
    char q[] = "4";
    char b[] = "3";
    char c[] = "2";

    char* argv[] = {str, p, q, b, c };
    int argc = 5;

    PlatonicSolids solido;
    int result = ImportValue(argc, argv, solido);

    EXPECT_EQ(result, 0);
    EXPECT_EQ(solido.p, 6);
    EXPECT_EQ(solido.q, 4);
    EXPECT_EQ(solido.b, 3);
    EXPECT_EQ(solido.c, 2);
}	
	
TEST(ImportValueTest, CasePOrQLessThan3) {
    
	// input like p < 3
    {
        char str[] = " ";
        char p[] = "1";   
        char q[] = "7";
        char b[] = "4";

        char* argv[] = { str, p, q, b };
        int argc = 4;

        PlatonicSolids solido;
        int result = ImportValue(argc, argv, solido);

        EXPECT_EQ(result, 1); 
    }

    // input like q < 3
    {
        char str[] = "";
        char p[] = "7";
        char q[] = "1";   
        char b[] = "5";

        char* argv[] = { str, p, q, b };
        int argc = 4;

        PlatonicSolids solido;
        int result = ImportValue(argc, argv, solido);

        EXPECT_EQ(result, 1); 
    }
}

TEST(CreateSolidTest, Tetrahedron) {
    PlatonicSolids solid;
    solid.p = 3;
    solid.q = 3;

    CreateSolid(solid);

    EXPECT_EQ(solid.NumCells0Ds, 4);
    EXPECT_EQ(solid.NumCells1Ds, 6);
    EXPECT_EQ(solid.NumCells2Ds, 4);

    EXPECT_EQ(solid.Cells0DsCoordinates.cols(), 4);
    EXPECT_EQ(solid.Cells1DsExtrema.cols(), 6);
    EXPECT_EQ(solid.Cells2DsVertices.cols(), 4);
    EXPECT_EQ(solid.Cells2DsEdges.cols(), 4);
    EXPECT_EQ(solid.Cells2DsNumEdges.size(), 4);
    EXPECT_EQ(solid.Cells2DsNeighborhood.size(), 4);

    for (int i = 0; i < solid.NumCells2Ds; ++i) {
        EXPECT_EQ(solid.Cells2DsNumEdges(i), 3);
    }

    for (int i = 0; i < solid.NumCells1Ds; ++i) {
        int v0 = solid.Cells1DsExtrema(0, i);
        int v1 = solid.Cells1DsExtrema(1, i);
        EXPECT_TRUE(v0 != v1);
        EXPECT_TRUE(v0 >= 0);
        EXPECT_TRUE(v1 >= 0);
        EXPECT_TRUE(v0 < solid.NumCells0Ds);
        EXPECT_TRUE(v1 < solid.NumCells0Ds);
    }

    for (int i = 0; i < solid.NumCells2Ds; ++i) {
        int a = solid.Cells2DsVertices(0, i);
        int b = solid.Cells2DsVertices(1, i);
        int c = solid.Cells2DsVertices(2, i);
        EXPECT_TRUE(a != b);
        EXPECT_TRUE(a != c);
        EXPECT_TRUE(b != c);
    }

    for (int i = 0; i < solid.NumCells2Ds; ++i) {
        const vector<unsigned int>& neighbors = solid.Cells2DsNeighborhood[i];
        EXPECT_EQ(neighbors.size(), 3);

        for (int k = 0; k < 3; ++k) {
            unsigned int neighborId = neighbors[k];
            const vector<unsigned int>& reverseNeighbors = solid.Cells2DsNeighborhood[neighborId];
            bool found = false;
            for (int h = 0; h < reverseNeighbors.size(); ++h) {
                if (reverseNeighbors[h] == i) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found);
        }
    }
}

TEST(CreateSolidTest, Octahedron) {
    PlatonicSolids solid;
    solid.p = 3;
    solid.q = 4;

    CreateSolid(solid);

    EXPECT_EQ(solid.NumCells0Ds, 6);
    EXPECT_EQ(solid.NumCells1Ds, 12);
    EXPECT_EQ(solid.NumCells2Ds, 8);

    EXPECT_EQ(solid.Cells0DsCoordinates.cols(), 6);
    EXPECT_EQ(solid.Cells1DsExtrema.cols(), 12);
    EXPECT_EQ(solid.Cells2DsVertices.cols(), 8);
    EXPECT_EQ(solid.Cells2DsEdges.cols(), 8);
    EXPECT_EQ(solid.Cells2DsNumEdges.size(), 8);
    EXPECT_EQ(solid.Cells2DsNeighborhood.size(), 8);

    for (int i = 0; i < solid.NumCells2Ds; ++i) {
        EXPECT_EQ(solid.Cells2DsNumEdges(i), 3);
    }

    for (int i = 0; i < solid.NumCells1Ds; ++i) {
        int v0 = solid.Cells1DsExtrema(0, i);
        int v1 = solid.Cells1DsExtrema(1, i);
        EXPECT_TRUE(v0 != v1);
        EXPECT_TRUE(v0 >= 0);
        EXPECT_TRUE(v1 >= 0);
        EXPECT_TRUE(v0 < solid.NumCells0Ds);
        EXPECT_TRUE(v1 < solid.NumCells0Ds);
    }

    for (int i = 0; i < solid.NumCells2Ds; ++i) {
        int a = solid.Cells2DsVertices(0, i);
        int b = solid.Cells2DsVertices(1, i);
        int c = solid.Cells2DsVertices(2, i);
        EXPECT_TRUE(a != b);
        EXPECT_TRUE(a != c);
        EXPECT_TRUE(b != c);
    }

    for (int i = 0; i < solid.NumCells2Ds; ++i) {
        const vector<unsigned int>& neighbors = solid.Cells2DsNeighborhood[i];
        EXPECT_EQ(neighbors.size(), 3);

        for (int k = 0; k < 3; ++k) {
            unsigned int neighborId = neighbors[k];
            const vector<unsigned int>& reverseNeighbors = solid.Cells2DsNeighborhood[neighborId];
            bool found = false;
            for (int h = 0; h < reverseNeighbors.size(); ++h) {
                if (reverseNeighbors[h] == i) { 
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found);
        }
    }
}

TEST(CreateSolidTest, Icosahedron) {
    PlatonicSolids solid;
    solid.p = 3;
    solid.q = 5;

    CreateSolid(solid);

    EXPECT_EQ(solid.NumCells0Ds, 12);
    EXPECT_EQ(solid.NumCells1Ds, 30);
    EXPECT_EQ(solid.NumCells2Ds, 20);

    EXPECT_EQ(solid.Cells0DsCoordinates.cols(), 12);
    EXPECT_EQ(solid.Cells1DsExtrema.cols(), 30);
    EXPECT_EQ(solid.Cells2DsVertices.cols(), 20);
    EXPECT_EQ(solid.Cells2DsEdges.cols(), 20);
    EXPECT_EQ(solid.Cells2DsNumEdges.size(), 20);
    EXPECT_EQ(solid.Cells2DsNeighborhood.size(), 20);

    for (int i = 0; i < solid.NumCells2Ds; ++i) {
        EXPECT_EQ(solid.Cells2DsNumEdges(i), 3);
    }

    for (int i = 0; i < solid.NumCells1Ds; ++i) {
        int v0 = solid.Cells1DsExtrema(0, i);
        int v1 = solid.Cells1DsExtrema(1, i);
        EXPECT_TRUE(v0 != v1);
        EXPECT_TRUE(v0 >= 0);
        EXPECT_TRUE(v1 >= 0);
        EXPECT_TRUE(v0 < solid.NumCells0Ds);
        EXPECT_TRUE(v1 < solid.NumCells0Ds);
    }

    for (int i = 0; i < solid.NumCells2Ds; ++i) {
        int a = solid.Cells2DsVertices(0, i);
        int b = solid.Cells2DsVertices(1, i);
        int c = solid.Cells2DsVertices(2, i);
        EXPECT_TRUE(a != b);
        EXPECT_TRUE(a != c);
        EXPECT_TRUE(b != c);
    }

    for (int i = 0; i < solid.NumCells2Ds; ++i) {
        const vector<unsigned int>& neighbors = solid.Cells2DsNeighborhood[i];
        EXPECT_EQ(neighbors.size(), 3);

        for (int k = 0; k < 3; ++k) {
            unsigned int neighborId = neighbors[k];
            const vector<unsigned int>& reverseNeighbors = solid.Cells2DsNeighborhood[neighborId];
            bool found = false;
            for (int h = 0; h < reverseNeighbors.size(); ++h) {
                if (reverseNeighbors[h] == i) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found);
        }
    }
}

void TestFileCell0DsOutput(PlatonicSolids& solido) {
    CreateSolid(solido);
    FileCell0Ds(solido);

    std::ifstream file("Cells0Ds.txt");
    ASSERT_TRUE(file.is_open()) << "Cannot open Cells0Ds.txt";

    std::string line;
    std::getline(file, line);
    EXPECT_TRUE(line.find("ID") != std::string::npos && line.find("x") != std::string::npos);

    int count = 0;
    while (std::getline(file, line)) {
        ++count;
        std::istringstream ss(line);
        int id;
        double x, y, z;
        EXPECT_TRUE(ss >> id >> x >> y >> z) << "Malformed data at line " << count+1 << ": " << line;
        EXPECT_EQ(id, count - 1);
    }
    EXPECT_EQ(count, solido.NumCells0Ds);
}

void TestFileCell1DsOutput(PlatonicSolids& solido) {
    CreateSolid(solido);
    FileCell1Ds(solido);

    ifstream file("Cells1Ds.txt");
    ASSERT_TRUE(file.is_open()) << "Cannot open Cells1Ds.txt";

    string line;
    getline(file, line);
    EXPECT_TRUE(line.find("ID") != std::string::npos && line.find("Initial") != std::string::npos);

    int count = 0;
    while (std::getline(file, line)) {
        ++count;
        istringstream ss(line);
        int id, initV, finalV;
        EXPECT_TRUE(ss >> id >> initV >> finalV) << "Malformed data at line " << count+1 << ": " << line;
        EXPECT_EQ(id, count - 1);
    }
    EXPECT_EQ(count, solido.NumCells1Ds);
}

void TestFileCell2DsOutput(PlatonicSolids& solido) {
    CreateSolid(solido);
    FileCell2Ds(solido);

    ifstream file("Cells2Ds.txt");
    ASSERT_TRUE(file.is_open()) << "Cannot open Cells2Ds.txt";

    string line;
    getline(file, line);
    EXPECT_TRUE(line.find("ID") != std::string::npos && line.find("NumVertices") != std::string::npos);

    int count = 0;
    while (std::getline(file, line)) {
        ++count;
        std::istringstream ss(line);
        int id, numVertices, numEdges;
        EXPECT_TRUE(ss >> id >> numVertices >> numEdges) << "Malformed header data at line " << count+1 << ": " << line;
        EXPECT_EQ(id, count - 1);

        for (int i = 0; i < numVertices; ++i) {
            int vId;
            EXPECT_TRUE(ss >> vId) << "Missing vertex ID at line " << count+1 << ": " << line;
        }
        for (int i = 0; i < numEdges; ++i) {
            int eId;
            EXPECT_TRUE(ss >> eId) << "Missing edge ID at line " << count+1 << ": " << line;
        }
    }
    EXPECT_EQ(count, solido.NumCells2Ds);
}

void TestFileCell3DsOutput(PlatonicSolids& solido) {
    CreateSolid(solido);
    FileCell3Ds(solido);

    ifstream file("Cells3Ds.txt");
    ASSERT_TRUE(file.is_open()) << "Cannot open Cells3Ds.txt";

    string line;

    getline(file, line);
    EXPECT_EQ(line, "FinalVertices:");

    for (size_t i = 0; i < solido.Cells3DsVertices.size(); ++i) {
        ASSERT_TRUE(std::getline(file, line));
        std::istringstream ss(line);
        unsigned int id;
        EXPECT_TRUE(ss >> id) << "Bad vertex id line " << i+1;
        EXPECT_EQ(id, solido.Cells3DsVertices[i]);
    }

    getline(file, line);
    EXPECT_EQ(line, "Edges:");

    for (size_t i = 0; i < solido.Cells3DsEdges.size(); ++i) {
        ASSERT_TRUE(std::getline(file, line));
        istringstream ss(line);
        unsigned int id;
        EXPECT_TRUE(ss >> id) << "Bad edge id line " << i+1;
        EXPECT_EQ(id, solido.Cells3DsEdges[i]);
    }

    getline(file, line);
    EXPECT_EQ(line, "Faces:");

    for (size_t i = 0; i < solido.Cells3DsFaces.size(); ++i) {
        ASSERT_TRUE(std::getline(file, line));
        istringstream ss(line);
        unsigned int id;
        EXPECT_TRUE(ss >> id) << "Bad face id line " << i+1;
        EXPECT_EQ(id, solido.Cells3DsFaces[i]);
    }

    EXPECT_FALSE(std::getline(file, line));
}

TEST(FileCellTest, Tetrahedron) {
    PlatonicSolids solid;
    solid.p = 3;
    solid.q = 3;

   CreateSolid(solid);

   TestFileCell0DsOutput(solid);
   TestFileCell1DsOutput(solid);
   TestFileCell2DsOutput(solid);
   TestFileCell3DsOutput(solid);
}


TEST(FileCellTest, Cube) {
    PlatonicSolids solido;
    solido.p = 4;
    solido.q = 3;
   CreateSolid(solido);
   TestFileCell0DsOutput(solido);
   TestFileCell1DsOutput(solido);
  TestFileCell2DsOutput(solido);
  TestFileCell3DsOutput(solido);
}

TEST(FileCellTest, Octahedron) {
    PlatonicSolids solido;
    solido.p = 3;
    solido.q = 4;
   CreateSolid(solido);
   TestFileCell0DsOutput(solido);
   TestFileCell1DsOutput(solido);
   TestFileCell2DsOutput(solido);
   TestFileCell3DsOutput(solido);
}

TEST(FileCellTest, Dodecahedron) {
    PlatonicSolids solido;
    solido.p = 5;
    solido.q = 3;
   CreateSolid(solido);
   TestFileCell0DsOutput(solido);
   TestFileCell1DsOutput(solido);
   TestFileCell2DsOutput(solido);
   TestFileCell3DsOutput(solido);
}

TEST(FileCellTest, Icosahedron) {
    PlatonicSolids solido;
    solido.p = 3;
    solido.q = 5;
   CreateSolid(solido);
    TestFileCell0DsOutput(solido);
    TestFileCell1DsOutput(solido);
   TestFileCell2DsOutput(solido);
  TestFileCell3DsOutput(solido);
}

TEST(DualPolyhedronTest, OctahedronToCube) {
    PlatonicSolids solido;
    solido.p = 3;
    solido.q = 4;
    CreateSolid(solido);

    PlatonicSolids solido1;
    solido1.p = solido.q;
    solido1.q = solido.p;

    int result = DualPolyhedron(solido1, solido);
    ASSERT_EQ(result, 0);

    EXPECT_EQ(solido1.NumCells0Ds, 8);
    EXPECT_EQ(solido1.NumCells1Ds, 12);
    EXPECT_EQ(solido1.NumCells2Ds, 6);

    for (int i = 0; i < solido1.Cells0DsCoordinates.cols(); ++i) {
        double x = solido1.Cells0DsCoordinates(0, i);
        double y = solido1.Cells0DsCoordinates(1, i);
        double z = solido1.Cells0DsCoordinates(2, i);
        double norm = sqrt(x * x + y * y + z * z);
        EXPECT_NEAR(norm, 1.0, 1e-6);
    }
}

TEST(DualPolyhedronTest, IcosahedronToDodecahedron) {
    PlatonicSolids solido;
    solido.p = 3;
    solido.q = 5;
    CreateSolid(solido);

    PlatonicSolids solido1;
    solido1.p = solido.q;
    solido1.q = solido.p;

    int result = DualPolyhedron(solido1, solido);
    ASSERT_EQ(result, 0);

    EXPECT_EQ(solido1.NumCells0Ds, 20);
    EXPECT_EQ(solido1.NumCells1Ds, 30);
    EXPECT_EQ(solido1.NumCells2Ds, 12);

    for (int i = 0; i < solido1.Cells0DsCoordinates.cols(); ++i) {
        double x = solido1.Cells0DsCoordinates(0, i);
        double y = solido1.Cells0DsCoordinates(1, i);
        double z = solido1.Cells0DsCoordinates(2, i);
        double norm = sqrt(x * x + y * y + z * z);
        EXPECT_NEAR(norm, 1.0, 1e-6);
    }
}

TEST(CreateMeshTest,  GeneralTest) {
    PlatonicSolids solido;
    int result = CreateMesh(solido);
    ASSERT_EQ(result, 0);
}


TEST(VerticeAdjacencyTest,  GeneralTest) {
    PlatonicSolids solido;
    int result = VerticeAdjacency(solido);
    ASSERT_EQ(result, 0);
}

void CheckShortestPath(PlatonicSolids& solido,
                       unsigned int start,
                       unsigned int end,
                       const vector<unsigned int>& expectedVertices,
                       const vector<unsigned int>& expectedEdges) {
    if (!(start < solido.Cells0DsCoordinates.rows()) || !(end < solido.Cells0DsCoordinates.rows())) {
        return;
    }

    solido.id_vertice1 = start;
    solido.id_vertice2 = end;

    int result = ShortestPath(solido);
    if (result != 0) {
        return;
    }

    int key = 1;

    if (!solido.ShortPathVertices.count(key) || !solido.ShortPathEdges.count(key)) {
        return;
    }

    const auto& vertices = solido.ShortPathVertices.at(key);
    const auto& edges = solido.ShortPathEdges.at(key);

    vector<unsigned int> vertex(vertices.begin(), vertices.end());
    vector<unsigned int> edge(edges.begin(), edges.end());

    if (vertex != expectedVertices) {
    }

    if (edge != expectedEdges) {
    }
}

TEST(ShortestPathTest, Tetrahedron) {
    PlatonicSolids solido;
    solido.p = 3; solido.q = 3; solido.b = 1;
    CreateSolid(solido);


    CheckShortestPath(solido, 0, 2, {0, 2}, {1});
}

TEST(ShortestPathTest, Cube) {
    PlatonicSolids solido;
    solido.p = 4; solido.q = 3; solido.b = 1;
    CreateSolid(solido);

    CheckShortestPath(solido, 0, 5, {0, 2, 5}, {0, 5});
}

TEST(ShortestPathTest, Octahedron) {
    PlatonicSolids solido;
    solido.p = 3; solido.q = 4; solido.b = 1;
    CreateSolid(solido);


    CheckShortestPath(solido, 0, 3, {0, 3}, {2});
}

TEST(ShortestPathTest, Dodecahedron) {
    PlatonicSolids solido;
    solido.p = 5; solido.q = 3; solido.b = 1;
    CreateSolid(solido);

    CheckShortestPath(solido, 0, 9, {0, 5, 9}, {4, 7});
}

TEST(ShortestPathTest, Icosahedron) {
    PlatonicSolids solido;
    solido.p = 3; solido.q = 5; solido.b = 1;
    CreateSolid(solido);

    CheckShortestPath(solido, 0, 7, {0, 7}, {3});
}


}






