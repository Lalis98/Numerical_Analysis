using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex5
    {
        /// <summary>
        /// Performs Cholesky Factorization on a given positive definite matrix to produce a lower triangular matrix L.
        /// </summary>
        /// <param name="matrixA">The input positive definite matrix to be decomposed.</param>
        /// <returns>The lower triangular matrix resulting from the Cholesky decomposition.</returns>
        public static double[,] FactorizationCholesky(double[,] matrixA)
        {
            int n = matrixA.GetLength(0);
            double[,] matrixL = new double[n, n]; // Initialize the lower triangular matrix L

            // Validate if the input matrix is square
            if (matrixA.GetLength(0) != matrixA.GetLength(1))
            {
                throw new Exception("The Matrix A must be a square matrix!");
            }

            for (int i = 0; i < n; i++)       // Iterate over rows (i-th Row)
            {
                for (int j = 0; j <= i; j++)  // Iterate over columns (j-th Column)
                {
                    double sum = 0.0;

                    // Calculate the sum for the elements of the lower triangular matrix
                    for (int k = 0; k < j; k++)  // k < j
                    {
                        sum = sum + matrixL[i, k] * matrixL[j, k];
                    }

                    if (i == j)
                    {
                        // Diagonal elements calculation: L[i,i] = sqrt(A[i,i] - Σ L[i,k]*L[i,k])
                        matrixL[i, j] = Math.Sqrt(matrixA[i, i] - sum);
                    }
                    else
                    {
                        // Off-diagonal elements calculation: L[i,j] = (1.0 / L[j,j]) * (A[i,j] - Σ L[i,k]*L[j,k])
                        matrixL[i, j] = (1.0 / matrixL[j, j]) * (matrixA[i, j] - sum);
                    }
                }
            }

            MatrixMethods.PrintMatrix(matrixL); // Print Matrix L transpose

            return matrixL; // Return Matrix L
        }

        /// <summary>
        /// Solves a linear system using Cholesky factorization.
        /// </summary>
        /// <param name="matrixL">The lower triangular matrix L from the Cholesky factorization.</param>
        /// <param name="vectorB">The known vector in the equation Ax = b.</param>
        /// <returns>The solution vector x for the equation Ax = b.</returns>
        public static double[] SolveFactorizationCholeskyLinearSystem(double[,] matrixL, double[] vectorB)
        {
            int n = matrixL.GetLength(0); // Number of Elements in a row
            double[] x = new double[n];   // Initilize result vector x
            double[] y = new double[n];   // Initilize vector y

            double[,] matrixLT = Matrices.MatrixMethods.TransposeMatrix(matrixL);
            
            if (n != vectorB.Length)
            {
                throw new Exception("Size of matrix L and vecotr b should be equal!");
            }

            // Forward Substitution (L * y = b), with b known vector
            for (int i = 0; i < n; i++)    // i = 0, 1, 2,..., n
            {

                y[i] = vectorB[i];

                for (int j = 0; j < i; j++)   // j = 0, 1,..., i - 1
                {

                    y[i] = y[i] - matrixL[i, j] * y[j];  // Subtract the known mebmers of the equation

                }

                y[i] = y[i] / matrixL[i, i];  // Calculate the final value of y
            }

            // Backward Substitution (U * x = y), with y known vector
            for (int i = n - 1; i >= 0; i--)  // i = n-1, ..., 0
            {
                x[i] = y[i];

                for (int j = i + 1; j < n; j++)   // j = i+1, ..., n
                {
                    x[i] = x[i] - matrixLT[i, j] * x[j];  // Subtract the known mebmers of the equation
                }

                x[i] = x[i] / matrixLT[i, i];  // Calculate the final value of x
            }

/*            Console.WriteLine("Matrix L:");
            MatrixMethods.PrintMatrix(matrixL); // Print Matrix L transpose
            Console.WriteLine("Matrix L^T:");
            MatrixMethods.PrintMatrix(matrixLT); // Print Matrix L transpose
            Console.WriteLine("Vector y:");
            Matrices.VectorMethods.PrintVector(y);
            Console.WriteLine("Vector x:");
            Matrices.VectorMethods.PrintVector(x);*/


            return x;
        }

    }
}
