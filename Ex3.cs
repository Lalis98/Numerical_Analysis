using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MscNumericalLinearAlgebra.Matrices;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex3
    {
        /// <summary>
        /// Perform LU decomposition on a given square matrix.
        /// </summary>
        /// <param name="matrixA">The input square matrix.</param>
        /// <returns>
        /// Tuple containing the lower triangular matrix (L) and the upper triangular matrix (U).
        /// </returns>
        public static (double[,] matrixL, double[,] matrixU) FactorizationLU(double[,] matrixA)
        {
            int n = matrixA.GetLength(1);  // Length of Columns (n)
            int m = matrixA.GetLength(0);  // Length of Rows (m)

            double[,] matrixU = matrixA;  // Create matrix U
            double[,] matrixL = new double[m, n];  // Create matrix L
            double factor = 0;  // Create factor variable

            if (n != m)  // Check if the matrix is square
            {
                throw new Exception("The matrix should be square!");  
            }

            matrixL = MatrixMethods.CreateIdentityMatrix(n);  // Create Identity Matrix of size n

            for (int i = 0; i < n; i++)  
            {
                for (int j = i + 1; j < n ; j++)
                {
                    if (matrixU[i, i] == 0)  // If the element is 0, do not divide
                    {
                        throw new Exception("Cannot divide with 0, Singular Matrix!");
                    }
                    factor = matrixU[j, i] / matrixU[i, i];  // i.e. A10 / A00
                    matrixL[j, i] = factor;  // Assign the Aji as the factor

                    for (int k = i ;  k < n ; k++)  // Run all the row elements to subtract the quantity (factor * above element)   
                    {
                        matrixU[j, k] = matrixU[j, k] - factor * matrixU[i, k];  // i.e. A10 = A10 - factor * A00
                    }
                }
            }

/*            Console.WriteLine("L = ");
            MatrixMethods.PrintMatrix(matrixL);

            Console.WriteLine("U = ");
            MatrixMethods.PrintMatrix(matrixU);*/

            return (matrixL, matrixU);


        }

        /// <summary>
        /// Solves a linear system using LU factorization.
        /// </summary>
        /// <param name="matrixL">Lower triangular matrix from LU decomposition.</param>
        /// <param name="matrixU">Upper triangular matrix from LU decomposition.</param>
        /// <param name="vectorB">Input vector in the equation Ax = b.</param>
        /// <returns>The solution vector x for the linear system Ax = b.</returns>
        public static double[] SolveFactorizationLULinearSystem(double[,] matrixL, double[,] matrixU, double[] vectorB)
        {

            //  A * x = b
            //  1. L * y = b
            //  2. U * x = y
            //
            // x = U^-1 (L^-1 * b)
            int n = matrixL.GetLength(0);  // Length matrix L & U
            double[] y = new double[n];  // Create y vector
            double[] x = new double[n];  // Create x vector

            if (matrixL.GetLength(0) != matrixL.GetLength(1))  // Check if the matrix is square
            {
                throw new Exception("The matrix should be square!");
            }

            // Forward Substitution (L * y = b), with b known vector
            for (int i = 0; i < n; i++)  // i = 0, 1, 2, ..., n
            {
                y[i] = vectorB[i];

                for (int j = 0; j < i; j++)
                {
                    y[i] = y[i] - matrixL[i, j] * y[j];
                }
            }

            // Backward Substitution (U * x = y), with y known vector
            for (int i = n - 1; i >= 0; i--)  // i = n-1, ..., 1, 0
            {
                x[i] = y[i];  // for a 3x3, x[3] = y[3]
                for (int j = i + 1; j < n; j++)  // for j = 1, 2, ..., n
                {
                    x[i] = x[i] - matrixU[i, j] * x[j];  //  U[1,1]*x1 = y[i](=x[i]) - U[1,2]*x2 - U[1,3]*x3
                }
                x[i] = x[i] / matrixU[i, i];  // x1 = (...) / U[1,1]
            }

/*            VectorMethods.PrintVector(x);*/

            return x;
        }


    }
}
