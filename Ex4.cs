using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex4
    {
        /// <summary>
        /// Performs LU decomposition in-place on a given square matrix, storing the results in the original matrix.
        /// </summary>
        /// <param name="matrixA">The input square matrix to be decomposed.</param>
        /// <returns>The input matrix, transformed into the combined form of lower and upper triangular matrices (LU decomposition).</returns>
        public static double[,] FactorizationLUStoredInMatrixA(double[,] matrixA)
        {
            int n = matrixA.GetLength(1);  // Length of Columns (n)
            int m = matrixA.GetLength(0);  // Length of Rows (m)

            if (n != m)  // Check if the matrix is square
            {
                throw new Exception("The matrix should be square!");
            }

            for (int i = 0; i < n; i++)             // i = 0, 1, 2,..., n
            {
                for (int j = i + 1; j < n; j++)     // j = 1, 2, 3,..., n 
                {
                    if (matrixA[i, i] == 0)  // If the element is 0, do not divide
                    {
                        throw new Exception("Cannot divide with 0, Singular Matrix!");
                    }

                    double factor = matrixA[j, i] / matrixA[i, i];  // i.e. A10 / A00

                    for (int k = i; k < n; k++)  // Run all the row elements to subtract the quantity (factor * above element)   
                    {
                        matrixA[j, k] = matrixA[j, k] - factor * matrixA[i, k];  // i.e. A10 = A10 - factor * A00
                    }

                    matrixA[j, i] = factor;  // Assign the Aji as the factor of the L Matrix
                }
            }

            Console.WriteLine("-----------------------");
            Console.WriteLine("A into LU = ");
            MatrixMethods.PrintMatrix(matrixA);
            Console.WriteLine("-----------------------");

            return matrixA;


        }
    }
}
