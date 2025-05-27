#pragma once

#include <iostream>

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
    
	PlatonicSolids solido;
    solido.p = 3;
    solido.q = 3;

    CreateSolid(solido);

    EXPECT_EQ(solido.NumCells0Ds, 4);
    EXPECT_EQ(solido.NumCells1Ds, 6);
    EXPECT_EQ(solido.NumCells2Ds, 4);
    EXPECT_EQ(solido.Cells0DsId.size(), 4);
    EXPECT_EQ(solido.Cells1DsExtrema.size(), 6);
    EXPECT_EQ(solido.Cells2DsVertices.size(), 4);
    EXPECT_EQ(solido.Cells2DsEdges.size(), 4);
    EXPECT_EQ(solido.Cells2DsNumEdges.size(), 4);
}

TEST(CreateSolidTest, Octahedron) {
    PlatonicSolids solid;
    solid.p = 3;
    solid.q = 4;

    CreateSolid(solid);

    EXPECT_EQ(solid.NumCells0Ds, 6);
    EXPECT_EQ(solid.NumCells1Ds, 12);
    EXPECT_EQ(solid.NumCells2Ds, 8);
    EXPECT_EQ(solid.Cells1DsExtrema.size(), 12);
    EXPECT_EQ(solid.Cells2DsVertices.size(), 8);
    EXPECT_EQ(solid.Cells2DsEdges.size(), 8);
    EXPECT_EQ(solid.Cells2DsNumEdges.size(), 8);
}

TEST(CreateSolidTest, GeneratesIcosahedron) {
    PlatonicSolids solid;
    solid.p = 3;
    solid.q = 5;

    CreateSolid(solid);

    EXPECT_EQ(solid.NumCells0Ds, 12);
    EXPECT_EQ(solid.NumCells1Ds, 30);
    EXPECT_EQ(solid.NumCells2Ds, 20);
    EXPECT_EQ(solid.Cells1DsExtrema.size(), 30);
    EXPECT_EQ(solid.Cells2DsVertices.size(), 20);
    EXPECT_EQ(solid.Cells2DsEdges.size(), 20);
    EXPECT_EQ(solid.Cells2DsNumEdges.size(), 20);
}





	
	
}






