using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex10
    {
        public static (double[], int[]) FactorizationCholeskySkyline(double[] valuesMatrixA, int[] diagonalOffsetsMatrixA)
        {
            double[] a = { 0.0 };

            int n = diagonalOffsetsMatrixA.Length - 1;  // Size of Positive Definite Matrix (nxn)

            double[] valuesMatrixL = new double[valuesMatrixA.Length];  // Values of L
            int[] diagonalOffsetsMatrixL = diagonalOffsetsMatrixA;  // Diagonal Offsets of L

            for (int i = 0; i < n; i++)
            {
                int hi = diagonalOffsetsMatrixA[i + 1] - diagonalOffsetsMatrixA[i] - 1;  // Height of Column i
                int mi = i - hi;  // Free Height of Column i

                for (int j = mi; j < i; j++)
                {

                    int hj = diagonalOffsetsMatrixA[j + 1] - diagonalOffsetsMatrixA[j] - 1;  // Height of Column j
                    int mj = j - hj;   // Free Height of Column j

                    double sum1 = 0.0;

                    // Calculate the sum for the elements of the lower triangular matrix
                    for (int k = Math.Max(mi, mj); k < j; k++)  // k < j
                    {
                        int L_ik_index = diagonalOffsetsMatrixA[i] + (i - k);
                        int L_jk_index = diagonalOffsetsMatrixA[j] + (j - k);

                        sum1 = sum1 + valuesMatrixL[L_ik_index] * valuesMatrixL[L_jk_index];
                    }

                    int L_ij_index = diagonalOffsetsMatrixA[i] + (i - j);
                    int L_jj_index = diagonalOffsetsMatrixA[j];

                    valuesMatrixL[L_ij_index] = (valuesMatrixA[L_ij_index] - sum1) / valuesMatrixL[L_jj_index];

                }

                double sum2 = 0;
                for (int k = mi; k < i; k++)  // k < j
                {
                    int L_ik_index = diagonalOffsetsMatrixA[i] + (i - k);

                    sum2 = sum2 + valuesMatrixL[L_ik_index] * valuesMatrixL[L_ik_index];
                }

                int L_ii_index = diagonalOffsetsMatrixA[i];
                valuesMatrixL[L_ii_index] = Math.Sqrt(valuesMatrixA[L_ii_index] - sum2);
            }

            Console.Write("L = ");
            VectorMethods.PrintVector(valuesMatrixL);

            return (valuesMatrixL, diagonalOffsetsMatrixL);

        }
    }
}
