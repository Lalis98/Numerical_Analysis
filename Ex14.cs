using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex14
    {
        public static double[] GradientDescentWithSkyline(double[] valuesMatrixA, int[] diagonalOffsetsMatrixA, double[] vectorB, int maxIterations, double epsilon)
        {

            int currentIteration = maxIterations;   // Iteration Number
            double a;                               // Step Size of Method
            int n = vectorB.Length;                 // Dimension of Matrix
            double[] x = new double[n];             // Vector x(t)
            double[] x_new = new double[n];         // Vector x(t+1)
            double[] r = new double[n];             // Residual r(t) = b - A * x(t)

            double[] r0 = vectorB;
            double r0_norm = Matrices.VectorMethods.VectorEuclideanNorm(r0);
            double rt_new_norm = 0;

            for (int t = 0; t < maxIterations; t++) // t --> Number of Iterations
            {

                // r(t) = b - A * x(t)
                double[] Ax = MatrixMethods.MatrixVectorMultiplicationWithSkyline(valuesMatrixA, diagonalOffsetsMatrixA, x);
                for (int i = 0; i < n; i++)
                {
                    r[i] = vectorB[i] - Ax[i];  // r = b - A * x
                }

                // a = (r^T * r) / (r^T * A * r)
                // Numerator: r^T * r
                // Denominator: r^T * A * r

                double[] Ar = MatrixMethods.MatrixVectorMultiplicationWithSkyline(valuesMatrixA, diagonalOffsetsMatrixA, r); // A * r

                double num = VectorMethods.VectorDotProduct(r, r); // r^T * r
                double denom = VectorMethods.VectorDotProduct(r, Ar); // r^T * Ar

                a = num / denom; // a = Numerator / Denominator

                for (int i = 0; i < n; i++)
                {
                    x_new[i] = x[i] + a * r[i];  // x(t+1) = x(t) + a * r(t)
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||
                double[] Ax_new = MatrixMethods.MatrixVectorMultiplicationWithSkyline(valuesMatrixA, diagonalOffsetsMatrixA, x_new); //  A* x(t)

                // Calculate the Residual r(t) = b - A * x(t)
                for (int i = 0; i < n; i++)
                {
                    rt_new[i] = vectorB[i] - Ax_new[i];
                }

                rt_new_norm = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_new_norm / r0_norm <= epsilon)  // If ||r(t)|| / ||r0|| < ε ---> Break Loop
                {
                    currentIteration = t;  // Save current Break Iteration
                    break;
                }

                x = x_new;

            }

            Console.WriteLine("Residual ||rt|| = " + rt_new_norm);
            Console.WriteLine("Number of Iterations N = " + (currentIteration + 1));
            Console.Write("Solution x = ");
            Matrices.VectorMethods.PrintVector(x_new);

            return x_new;
        }


        public static double[] ConjugateGradientWithSkyline(double[] valuesMatrixA, int[] diagonalOffsetsMatrixA, double[] vectorB, int maxIterations, double epsilon)
        {
            int currentIteration = maxIterations;   // Current Number of Iterations of the for loop
            double a;                               // Step Size of Method
            double b;
            int n = vectorB.Length;                 // Dimension of Matrix
            double[] x = new double[n];             // Vector x(t)
            double[] x_new = new double[n];         // Vector x(t+1)
            double[] r = new double[n];             // Residual r(t) = b - A * x(t)
            double[] r_new = new double[n];         // Residual r(t+1)
            double[] d_new = new double[n];         // Residual d(t+1)
            double r_new_norm = 0;                  // Norm ||r(t+1)||


            double[] r0 = vectorB;                  // Residual at step r(t=0) = r0
            double r0_norm = Matrices.VectorMethods.VectorEuclideanNorm(r0); // Norm of Residual ||r0|| 
            double[] d = new double[n];                   // Direction Vector


            for (int t = 0; t < maxIterations; t++)
            {

                double[] Ax = MatrixMethods.MatrixVectorMultiplicationWithSkyline(valuesMatrixA, diagonalOffsetsMatrixA, x);

                // r(t) = b - A * x(t)
                for (int i = 0; i < n; i++)
                {
                    r[i] = vectorB[i] - Ax[i];  // r = b - A * x
                    d[i] = r[i];  // Save for later d(t)
                }

                // a = (r^T * r) / (d^T * A * d)
                // Numerator: r^T * r
                // Denominator: r^T * A * r

                double[] Ad = MatrixMethods.MatrixVectorMultiplicationWithSkyline(valuesMatrixA, diagonalOffsetsMatrixA, d);
                double num = VectorMethods.VectorDotProduct(r, r); // r^T * r
                double denom = VectorMethods.VectorDotProduct(Ad, d);  // (Ad)^T * d 

                a = num / denom; // a = Numerator / Denominator

                for (int i = 0; i < n; i++)
                {
                    x_new[i] = x[i] + a * d[i];  // x(t+1) = x(t) + a * r(t)
                    r_new[i] = r[i] - a * Ad[i]; // r(t+1) = r(t) - a * A * d(t)
                }

                r_new_norm = VectorMethods.VectorEuclideanNorm(r_new);

                if (r_new_norm / r0_norm <= epsilon)  // If ||r(t+1)|| / ||r0|| < ε ---> Break Loop
                {
                    currentIteration = t;  // Save current Break Iteration
                    break;
                }

                // b = (r^T * r) / (r^T * r)
                num = VectorMethods.VectorDotProduct(r_new, r_new);
                denom = VectorMethods.VectorDotProduct(r, r);
                b = num / denom;

                for (int i = 0; i < n; i++)
                {
                    d_new[i] = r_new[i] + b * d[i];
                }

                x = x_new;
                r = r_new;
                d = d_new;
            }

            Console.WriteLine("Residual ||rt|| = " + r_new_norm);
            Console.WriteLine("Number of Iterations N = " + (currentIteration + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(x_new);

            return x_new;
        }
    }
}
