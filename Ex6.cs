using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex6
    {
        public static double[,] FactorizationCholeskyStoredInMatrixA(double[,] matrixA)
        {
            if (matrixA.GetLength(0) != matrixA.GetLength(1))
            {
                throw new Exception("Matrix should be square!");
            }

            int n = matrixA.GetLength(0);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double sum = 0.0;

                    for (int k = 0; k < j; k++)  // Ensure k starts from 0
                    {
                        sum = sum + matrixA[i, k] * matrixA[j, k];
                    }

                    if (i == j)
                    {
                        matrixA[i, j] = Math.Sqrt(matrixA[i, i] - sum);
                    }
                    else
                    {
                        matrixA[i, j] = (1.0 / matrixA[j, j]) * (matrixA[i, j] - sum);
                        matrixA[j, i] = matrixA[i, j]; // Simultaneously assign L^T into the lower triangular part
                    }
                }
            }


            Matrices.MatrixMethods.PrintMatrix(matrixA);

            return matrixA;
        }
    }
}
