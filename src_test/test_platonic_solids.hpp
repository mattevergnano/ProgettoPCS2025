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
	

	
	
}






