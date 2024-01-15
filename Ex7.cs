using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex7
    {
        /// <summary>
        /// Solves the linear system represented by a sparse matrix A in Compressed Sparse Row (CSR) format,
        /// where A * x = b, and x is a vector.
        /// </summary>
        /// <param name="matrixA">The input sparse matrix A in dense format.</param>
        /// <param name="vectorX">The input vector x.</param>
        /// <returns>
        /// The resulting vector b after solving the linear system A * x = b.
        /// </returns>
        /// <exception cref="Exception">Thrown when the dimensions of the matrix and vector are incompatible for multiplication.</exception>
        public static double[] SolveLinearSystemMatrixVectorWithCSR(double[,] matrixA, double[] vectorX)
        {
            int numRowsVectorX = vectorX.Length;
            int numColumnsMatrixA = matrixA.GetLength(1);

            if (numRowsVectorX != numColumnsMatrixA) 
            {
                throw new Exception("Cannot Multiply Matrix A with Vector X with these dimensions!");
            }
            double[] resultVector = new double[matrixA.GetLength(0)];

            // Retrieve the CSR format arrays representing the transpose of matrixA.
            (double[] valuesArray, int[] colIndices, int[] rowOffsets) = Matrices.StoringMatrices.SaveMatrixCSR(matrixA);

            for (int i = 0; i < rowOffsets.Length - 1; i++)
            {
                resultVector[i] = 0.0;

                for (int j = rowOffsets[i]; j < rowOffsets[i + 1]; j++)
                {
                    resultVector[i] = resultVector[i] + valuesArray[j] * vectorX[colIndices[j]];
                }

            }

            Console.WriteLine("");
            Console.WriteLine("b = {" + string.Join(", ", resultVector) + "}");

            return resultVector;
        }

        /// <summary>
        /// Solves the linear system represented by the transpose of a sparse matrix A in Compressed Sparse Row (CSR) format,
        /// where A^T * x = b, and x is a vector.
        /// </summary>
        /// <param name="matrixA">The input matrix A in dense format.</param>
        /// <param name="vectorX">The input vector x.</param>
        /// <returns>
        /// The resulting vector b after solving the linear system A^T * x = b.
        /// </returns>
        /// <exception cref="Exception">Thrown when the dimensions of the matrix and vector are incompatible for multiplication.</exception>
        public static double[] SolveLinearSystemTransposeMatrixVectorWithCSR(double[,] matrixA, double[] vectorX)
        {

            int numRowsVectorX = vectorX.Length;
            int numColumnsMatrixA = matrixA.GetLength(0);
            int numRowsMatrixA = matrixA.GetLength(1);

            if (numRowsVectorX != numColumnsMatrixA)
            {
                throw new Exception("Cannot Multiply Matrix A with Vector X with these dimensions!");
            }
            double[] resultVector = new double[matrixA.GetLength(1)];

            // Retrieve the CSC format arrays representing the transpose of matrixAT.
            (double[] valuesArray, int[] colIndices, int[] rowOffsets) = Matrices.StoringMatrices.SaveMatrixCSR(matrixA);

            for (int i = 0; i < rowOffsets.Length - 1; i++)
            {
                for (int j = rowOffsets[i]; j < rowOffsets[i + 1]; j++)
                {
                    resultVector[colIndices[j]] = resultVector[colIndices[j]] + valuesArray[j] * vectorX[i];
                }

            }

            Console.WriteLine("");
            Console.WriteLine("b = {" + string.Join(", ", resultVector) + "}");

            return resultVector;
        }
    }
}
