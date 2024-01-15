using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class test
    {
        public static double[] SOR(double[,] matrixA, double[] vectorB, double epsilon, int maxIterations, double weight)
        {
            int currentIteration = maxIterations;
            int n = matrixA.GetLength(0);

            // Condition to check the dimensions of the matrices
            if (n != matrixA.GetLength(1) || n != vectorB.Length)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b");  // Print Message to fix the dimensions
            }

            // Initilization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_new_norm = 0;

            for (int t = 0; t < maxIterations; t++)
            {
                for (int i = 0; i < n; i++)
                {
                    double sum1 = 0;
                    double sum2 = 0;

                    for (int j = 0; j < i; j++)
                    {
                        sum1 = sum1 + matrixA[i, j] * xt_new[j];
                    }

                    for (int j = i + 1; j < n; j++)
                    {
                        sum2 = sum2 + matrixA[i, j] * xt[j];
                    }

                    xt_new[i] = (1 - weight) * xt[i] + (weight / matrixA[i, i]) * (vectorB[i] - sum1 - sum2);
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                // Calculate the r(t) = b - A * x(t)
                for (int row = 0; row < n; row++)
                {
                    double sum_row = 0;

                    for (int col = 0; col < n; col++)
                    {
                        sum_row += matrixA[row, col] * xt_new[col];
                    }

                    rt_new[row] = vectorB[row] - sum_row;
                }

                rt_new_norm = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_new_norm / rt_norm <= epsilon)  // If ||r(t)|| / ||r0|| < ε ---> Break Loop
                {
                    currentIteration = t;  // Save current Break Iteration
                    break;
                }

                xt = xt_new;
            }

            Console.WriteLine("Residual ||rt|| = " + rt_new_norm);
            Console.WriteLine("Number of Iterations N = " + (currentIteration + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);

            return xt_new;
        }


    }

}
