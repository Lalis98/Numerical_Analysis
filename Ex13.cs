using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex13
    {
        public static double[] GradientDescent(double[,] matrixA, double[] vectorB, int maxIterations, double epsilon)
        {
            int currentIteration = maxIterations;
            double a = 0;  // Step Size of Method
            int n = matrixA.GetLength(0); // Dimension of Matrix
            double[] x = new double[n];   // Vector x(t)
            double[] x_new = new double[n];  // Vector x(t+1)
            double[] r = new double[n];   // Residual r(t) = b - A * x(t)

            double[] r0 = vectorB;
            double r0_norm = Matrices.VectorMethods.VectorEuclideanNorm(r0);
            double rt_new_norm = 0;


            for (int t = 0; t < maxIterations; t++) // t --> Number of Iterations
            {

                // r(t) = b - A * x(t)
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;

                    for (int j = 0; j < n; j++)
                    {
                        sum = sum + matrixA[i, j] * x[j];  // A * x
                    }

                    r[i] = vectorB[i] - sum;  // r = b - A * x
                }

                // a = (r^T * r) / (r^T * A * r)
                // Numerator: r^T * r
                // Denominator: r^T * A * r

                double num = 0; // Numerator
                for (int i = 0; i < n; i++)
                {
                    num = num + r[i] * r[i]; // r^T * r
                }

                double denom = 0; // Denominator
                double[] q = new double[n];

                
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0;j < n; j++)
                    {
                        sum = sum + r[j] * matrixA[i, j];  // q = r^T * A
                    }

                    q[i] = sum;
                    
                    denom = denom + q[i] * r[i];  // q^T * r
                }

                a = num / denom; // a = Numerator / Denominator

                for (int i = 0; i < n ; i++)
                {
                    x_new[i] = x[i] + a * r[i];  // x(t+1) = x(t) + a * r(t)
                }


                double[] rt_new = new double[n];  // Residual ||r(t)||

                // Calculate the Residual r(t) = b - A * x(t)
                for (int row = 0; row < n; row++)
                {
                    double sum_row = 0;

                    for (int col = 0; col < n; col++)
                    {
                        sum_row += matrixA[row, col] * x_new[col];
                    }

                    rt_new[row] = vectorB[row] - sum_row;
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

        public static double[] GradientDescentWithCSR(double[] valuesArrayMatrixA, int[] colIndicesMatrixA, int[] rowOffsetsMatrixA, double[] vectorB, int maxIterations, double epsilon)
        {

            int currentIteration = maxIterations;  // Initilize Current Iterations of the for loop
            double a = 0;  // Step Size of Method
            int n = vectorB.Length; // Dimension of Matrix
            double[] x = new double[n];   // Vector x(t)
            double[] x_new = new double[n];  // Vector x(t+1)
            double[] r = new double[n];   // Residual r(t) = b - A * x(t)

            double[] r0 = vectorB;
            double r0_norm = Matrices.VectorMethods.VectorEuclideanNorm(r0);
            double rt_new_norm = 0;

            for (int t = 0; t < maxIterations ; t++)
            {
                // r(t) = b - A * x(t)

                double[] Ax = MatrixMethods.MatrixVectorMultiplicationWithCSR(valuesArrayMatrixA, colIndicesMatrixA, rowOffsetsMatrixA, x); // A * x
                for (int i = 0; i < n; i++)
                {
                    r[i] = vectorB[i] - Ax[i];  // r = b - A * x
                }

                // a = (r^T * r) / (r^T * A * r)
                // Numerator: r^T * r
                // Denominator: r^T * A * r = r^T * q

                double num = 0; // Numerator
                for (int i = 0; i < n; i++)
                {
                    num = num + r[i] * r[i]; // r^T * r
                }

                double denom = 0; // Denominator
                double[] q = MatrixMethods.MatrixVectorMultiplicationWithCSR(valuesArrayMatrixA, colIndicesMatrixA, rowOffsetsMatrixA, r); // q = A * r
                for (int i = 0; i < n; i++)
                {
                    denom = denom + q[i] * r[i];  // r^T * q
                }

                a = num / denom; // Step Size of Method

                for (int i = 0; i < n; i++)
                {
                    x_new[i] = x[i] + a * r[i];  // x(t+1) = x(t) + a * r(t)
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||
                double[] Ax_new = MatrixMethods.MatrixVectorMultiplicationWithCSR(valuesArrayMatrixA, colIndicesMatrixA, rowOffsetsMatrixA, x_new); // A * x 
                
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

        public static double[] ConjugateGradient(double[,] matrixA, double[] vectorB, int maxIterations, double epsilon)
        {
            int currentIteration = maxIterations;   // Current Number of Iterations of the for loop
            double a = 0;                           // Step Size of Method
            double b = 0;
            int n = matrixA.GetLength(0);           // Dimension of Matrix
            double[] x = new double[n];             // Vector x(t)
            double[] x_new = new double[n];         // Vector x(t+1)
            double[] r = new double[n];             // Residual r(t) = b - A * x(t)
            double[] r_new = new double[n];         // Residual r(t+1)
            double[] d_new = new double[n];         // Residual d(t+1)
            double r_new_norm = 0;                  // Norm ||r(t+1)||


            double[] r0 = vectorB;                  // Residual at step r(t=0) = r0
            double r0_norm = Matrices.VectorMethods.VectorEuclideanNorm(r0); // Norm of Residual ||r0|| 
            double rt_new_norm = 0;         // Initilization of norm ||r(t)||
            double[] d = new double[n];                   // Direction Vector

            for (int t = 0; t < maxIterations; t++)
            {
                // r(t) = b - A * x(t)
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;

                    for (int j = 0; j < n; j++)
                    {
                        sum = sum + matrixA[i, j] * x[j];  // A * x
                    }

                    r[i] = vectorB[i] - sum;  // r = b - A * x
                    d[i] = r[i];
                }

                // a = (r^T * r) / (d^T * A * d)
                // Numerator: r^T * r
                // Denominator: r^T * A * r

                double num = 0; // Numerator
                for (int i = 0; i < n; i++)
                {
                    num = num + r[i] * r[i]; // r^T * r
                }

                double denom = 0; // Denominator
                double[] Ad = new double[n];


                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        sum = sum + d[j] * matrixA[i, j];  // (Ad)^T = d^T * A
                    }

                    Ad[i] = sum;

                    denom = denom + d[i] * Ad[i];  // (Ad)^T * d 
                }

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


                num = 0;
                denom = 0;
                // b = (r^T * r) / (r^T * r)
                for (int i = 0; i < n; i++)
                {
                    num = num + r_new[i] * r_new[i];
                    denom = denom + r[i] * r[i];
                }

                b = num / denom;

                for (int i = 0;i < n; i++)
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
            Matrices.VectorMethods.PrintVector(x_new);

            return x_new;
        }

        public static double[] ConjugateGradientWithCSR(double[] valuesArrayMatrixA, int[] colIndicesMatrixA, int[] rowOffsetsMatrixA, double[] vectorB, int maxIterations, double epsilon)
        {
            int currentIteration = maxIterations;   // Current Number of Iterations of the for loop
            double a = 0;                           // Step Size of Method
            double b = 0;
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

                double[] Ax = MatrixMethods.MatrixVectorMultiplicationWithCSR(valuesArrayMatrixA, colIndicesMatrixA, rowOffsetsMatrixA, x);

                // r(t) = b - A * x(t)
                for (int i = 0; i < n; i++)
                {
                    r[i] = vectorB[i] - Ax[i];  // r = b - A * x
                    d[i] = r[i];  // Save for later d(t)
                }

                // a = (r^T * r) / (d^T * A * d)
                // Numerator: r^T * r
                // Denominator: r^T * A * r


                double num = VectorMethods.VectorDotProduct(r, r); // r^T * r

                double[] Ad = MatrixMethods.MatrixVectorMultiplicationWithCSR(valuesArrayMatrixA, colIndicesMatrixA, rowOffsetsMatrixA, d);

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
            Matrices.VectorMethods.PrintVector(x_new);

            return x_new;
        }

    }
}
