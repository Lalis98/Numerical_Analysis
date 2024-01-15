using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex11
    {
        /// <summary>
        /// Solves a linear system of equations Ax = b using the Jacobi iterative method.
        /// </summary>
        /// <param name="matrixA">The coefficient matrix A in the system Ax = b.</param>
        /// <param name="vectorB">The constant vector b in the system Ax = b.</param>
        /// <param name="epsilon">The tolerance for convergence. Iterations stop when the relative residual is less than or equal to this value.</param>
        /// <param name="maxIterations">The maximum number of iterations to perform.</param>
        /// <returns>The solution vector x.</returns>
        /// <exception cref="Exception">Thrown when the dimensions of the input matrices do not match the requirements of a linear system.</exception>
        public static double[] Jacobi(double[,] matrixA, double[] vectorB, double epsilon, int maxIterations)
        {

            int current_t = maxIterations;  // In case the for loop makes all the iterations
            int n = matrixA.GetLength(0);  // dimension of the matrix A (also for vectorB)

            // Condition to check the dimensions of the matrices
            if (n != matrixA.GetLength(1) || n != vectorB.Length)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b");  // Print Message to fix the dimensions
            }

            // Initilization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_norm_new = 0;

            for (int t = 0; t < maxIterations; t++)  // t: Current Number of Iteration
            {

                for (int i = 0; i < n; i++)  // run all the rows of Aij
                {
                    // x(t+1) = [D^-1 * L * x(t)] + [D^-1 * U * x(t)] + [D^-1 * b]
                    // x(t+1) =        [s1]       +        [s2]       + [D^-1 * b]
                    double s1 = 0;
                    double s2 = 0;

                    for (int j = 0; j < i; j++)  // run all the elements of a row of Aij
                    {
                        // x(i) = D^-1 * L * x(j)
                        s1 = s1 + 1 / matrixA[i, i] * (-matrixA[i, j] * xt[j]);  // Sum of xi for L
                    }

                    for (int j = i + 1; j < n; j++)  // run all the elements of a row of Aij
                    {
                        // x(i) = D^-1 * U * x(j)
                        s2 = s2 + 1 / matrixA[i, i] * (-matrixA[i, j] * xt[j]);  // Sum of xi for D
                    }

                    // x = D^-1 * (L + U) x + D^-1 * b
                    xt_new[i] = s1 + s2 + 1 / matrixA[i, i] * vectorB[i];  // xi = (s1 + s2) + 1 / Aii * bi

                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                // Calculate the r(t) = b - A * x(t)
                for (int row = 0; row < n; row++)
                {
                    double sum_row = 0;

                    for (int col = 0; col < n; col++)
                    {
                        sum_row = sum_row + matrixA[row, col] * xt_new[col];
                    }
                    rt_new[row] = vectorB[row] - sum_row;
                }

                rt_norm_new = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_norm_new / rt_norm <= epsilon)  // If ||r(t)|| / ||r0|| < ε ---> Break Loop
                {
                    current_t = t;  // Save current Iteration
                    break;
                }

                xt = xt_new; // Assign x(t) = x(t+1), for the next loop
            }

            Console.WriteLine("Residual ||rt|| = " + rt_norm_new);
            Console.WriteLine("Number of Iterations N = " + (current_t + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);
           
            return xt_new; // Return Solution x
        }

    }
}
