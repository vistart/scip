/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sudoku_utils.cpp
 * @brief  A set of utilities that are used to read the puzzle and display the puzzle
 * @author Naga V C Gudapati
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

namespace sudoku
{
   /** reads in the sudoku puzzle from filepath
    *
     * Reads the string of sudoku puzzle into a 9x9 grid represented by a vector
     * of a vector of ints. The actual number is stored as itself and the blanks are stored as -1.
     *
     */
   std::vector<std::vector<int>> getSudokuPuzzle( std::string &filepath )
   {
      /* setting up a 9x9 grid forstoring the sudoku puzzle. */
      std::vector<std::vector<int>> puzzle(9, std::vector<int>(9));

      /* Reading the puzzle into a stringstream */
      std::ifstream infile(filepath);

      std::string puzzledata = "";

      if( infile.is_open() )
      {
         std::getline( infile, puzzledata );
         if( puzzledata.length() != 81 ) /* The puzzle should have 81 characters */
         {
            std::cerr << "Please check the puzzle file forinconsistencies"
                      << "\n";
            exit(1);
         }
      }

      int idx = 0; /* This variable will be used to access the numbers in the puzzle string */

      for( int i = 0; i < 9; ++i )
      {
         for( int j = 0; j < 9; ++j )
         {
            /* We will only convert the numeric string to an integer if it is not '.' or '0'. */
            if( (puzzledata.substr(idx, 1) != ".") && (puzzledata.substr(idx, 1) != "0") )
            {
               puzzle[i][j] = std::stoi( puzzledata.substr(idx, 1) );
            }
            else
            {
               /* If we are currently reading a '.' or '0' make it -1. */
               puzzle[i][j] = -1;
            }
            idx++;
         }
      }

      return puzzle;
   }

   /** prints the sudoku puzzle to console */
   void printSudoku( const std::vector<std::vector<int>> &sudokupuzzle )
   {
      std::cout << "+----------+-----------+-----------+" << "\n";
      for( int i = 0; i < 9; ++i )
      {
         std::cout << "|";
         for( int j = 0; j < 9; ++j )
         {
            if( sudokupuzzle[i][j] > 0 )
            {

               if( j == 2 || j == 5 || j == 8 )
               {
                  std::cout << sudokupuzzle[i][j] << " | ";
               }
               else
               {
                  std::cout << sudokupuzzle[i][j] << "   ";
               }
            }
            else
            {
               if( j == 2 || j == 5 || j == 8 )
               {
                  std::cout << "*" << " | ";
               }
               else
               {
                  std::cout << "*" << "   ";
               }
            }
         }
         std::cout << "\n";

         if( i == 2 || i == 5 || i == 8 )
         {
            std::cout << "+----------+-----------+-----------+" << "\n";
         }
      }
   }
} /* namespace sudoku */
