using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using MscNumericalLinearAlgebra.Matrices;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex2
    {
        public static double[,] SolveLinearSystemGaussElimination(double[,] matrixA, double[,] matrixB)
        {

            int lenColMatrixA = matrixA.GetLength(1);
            int lenRowMatrixA = matrixA.GetLength(0);
            int lenColMatrixB = matrixB.GetLength(1);
            int lenRowMatrixB = matrixB.GetLength(0);

            double[,] vectorX = new double[lenRowMatrixA, lenColMatrixB];

            int maxRow = 0;

            Console.WriteLine("System of Equations Before Gaussian Elimination with Partial Pivoting:");
            PrintSystem(matrixA, matrixB);

            for (int k = 0; k < lenRowMatrixA - 1; k++) // k is the column we are making pivoting
            {
                // Find the row with the largest pivot in the current column (k)
                double maxPivot = Math.Abs(matrixA[k, k]);  // Set maxPivot as the first Point 
                maxRow = k;  // Set the k index row as the maxPivot row

                for (int i = k + 1; i < lenRowMatrixA; i++) // Run the rest of the columns to compare maxPivot with the rest of the rows
                {
                    if (Math.Abs(matrixA[i, k]) > maxPivot)  // If the current maxPivot is less than the next row
                    {
                        maxPivot = Math.Abs(matrixA[i, k]);  // Assign the new maxPivot
                        maxRow = i;  // Assign the row index of the maxPivot
                    }
                }

                // Swap the current row with the row containing the largest pivot
                if (maxRow != k)  // If the maxRow is not the initial one
                {
                    for (int j = 0; j < lenColMatrixA; j++)  // run all the elements of the maxRow
                    {
                        // A[k, j] <-- Change elements --> A[maxRow, j]
                        double tempA = matrixA[k, j];
                        matrixA[k, j] = matrixA[maxRow, j];
                        matrixA[maxRow, j] = tempA;
                    }

                    for (int j = 0; j < lenColMatrixB; j++)  // run all the elements of the maxRow
                    {
                        // B[k, j] < --Change elements-- > B[maxRow, j]
                        double tempB = matrixB[k, j];
                        matrixB[k, j] = matrixB[maxRow, j];
                        matrixB[maxRow, j] = tempB;
                    }
                }

                // Perform Gaussian elimination
                for (int i = k + 1; i < lenRowMatrixA; i++)  // Run all rows under the 1st one
                {
                    double factor = matrixA[i, k] / matrixA[k, k];  // find the coefficient factor (Aik / Akk)
                    for (int j = k; j < lenColMatrixA; j++)  // Run all the columns from k-th element of the ith-row
                    {
                        matrixA[i, j] =  - matrixA[i, j] +  factor * matrixA[k, j]; // multiply with this so we get 0 in the current one
                    }

                    for (int j = 0; j < lenColMatrixB; j++)  // Run all the columns from k-th elemnt of the ith-row
                    {
                        matrixB[i, j] = - matrixB[i, j] + factor * matrixB[k, j]; // make the same multiplication for Matrix B
                    }
                }
            }

            // Back substitution
            for (int i = lenRowMatrixA - 1; i >= 0; i--) // from (m-1) to 0 for matrix A
            {
                for (int j = 0; j < lenColMatrixB; j++)  // j running all the colums of matrix B
                {
                    double sum = 0.0;  // Initialize sum
                    for (int k = i + 1; k < lenRowMatrixA; k++)  // 
                    {
                        sum += matrixA[i, k] * vectorX[k, j];
                    }
                    vectorX[i, j] = (matrixB[i, j] - sum) / matrixA[i, i];
                }
            }

            Console.WriteLine("Solution Vector (vectorX):");

            MatrixMethods.PrintMatrix(vectorX);

            // Function to print the system of equations
            void PrintSystem(double[,] A, double[,] B)
            {
                int rows = A.GetLength(0);
                int cols = A.GetLength(1);

                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        Console.Write(A[i, j] + "x" + j + " ");
                        if (j < cols - 1)
                            Console.Write("+ ");
                    }
                    Console.Write("= ");
                    for (int j = 0; j < B.GetLength(1); j++)
                    {
                        Console.Write(B[i, j] + " ");
                    }
                    Console.WriteLine();
                }
            }

            return vectorX;
        }
    }
}
