using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex8
    {
        /// <summary>
        /// Solves a linear system represented by a sparse matrix A in COO (Coordinate List) format,
        /// where A * x = b. The matrix A is assumed to be stored in row-major or column-major order
        /// based on the specified storing method.
        /// </summary>
        /// <param name="matrixA">The input sparse matrix A in dense format.</param>
        /// <param name="vectorX">The input dense vector x.</param>
        /// <param name="storingMethodMajor">
        /// The storing method for the COO format. 
        /// Valid values are "Row" (default) or "Column".
        /// </param>
        /// <returns>
        /// The resulting dense vector b after solving the linear system A * x = b.
        /// </returns>
        /// <exception cref="Exception">Thrown when the dimensions of the matrix and vector are incompatible for multiplication.</exception>
        public static double[] SolveLinearSystemMatrixVectorWithCOO(double[,] matrixA, double[] vectorX, string storingMethodMajor = "Row")
        {
            int numRowsMatrixA = matrixA.GetLength(0);
            int numColsMatrixA = matrixA.GetLength(1);
            int numVectorX = vectorX.Length;

            if (numColsMatrixA != numVectorX)
            {
                throw new Exception("Cannot Multiply Matrix A with Vector X with these dimensions!");
            }

            // Row Major COO
            if (storingMethodMajor == "Row")
            {
                (double[] valuesArray, int[] rowsArray, int[] colsArray) = Matrices.StoringMatrices.SaveMatrixCOO(matrixA, storingMethodMajor);

                double[] resultVector = new double[numRowsMatrixA];

                for (int i = 0; i < valuesArray.Length; i++)
                {
                    double value = valuesArray[i];
                    int rowIndex = rowsArray[i];
                    int colIndex = colsArray[i];

                    resultVector[rowIndex] = resultVector[rowIndex] + value * vectorX[colIndex];
                }

/*                Console.WriteLine("");
                Console.WriteLine("b = {" + string.Join(", ", resultVector) + "}");*/

                return resultVector;
            }

            // Column Major COO
            else
            {
                (double[] valuesArray, int[] colsArray, int[] rowsArray) = Matrices.StoringMatrices.SaveMatrixCOO(matrixA, storingMethodMajor);

                double[] resultVector = new double[numRowsMatrixA];

                for (int i = 0; i < valuesArray.Length; i++)
                {
                    double value = valuesArray[i];
                    int colIndex = colsArray[i];
                    int rowIndex = rowsArray[i];

                    resultVector[rowIndex] = resultVector[rowIndex] + value * vectorX[colIndex];
                }

/*                Console.WriteLine("");
                Console.WriteLine("b = {" + string.Join(", ", resultVector) + "}");*/

                return resultVector;

            }

        }

        /// <summary>
        /// Solves a linear system represented by the transpose of a sparse matrix A in COO (Coordinate List) format,
        /// where A^T * x = b. The matrix A is assumed to be stored in row-major or column-major order
        /// based on the specified storing method.
        /// </summary>
        /// <param name="matrixA">The input sparse matrix A in dense format.</param>
        /// <param name="vectorX">The input dense vector x.</param>
        /// <param name="storingMethodMajor">
        /// The storing method for the COO format. 
        /// Valid values are "Row" (default) or "Column".
        /// </param>
        /// <returns>
        /// The resulting dense vector b after solving the linear system A^T * x = b.
        /// </returns>
        /// <exception cref="Exception">Thrown when the dimensions of the matrix and vector are incompatible for multiplication.</exception>
        public static double[] SolveLinearSystemTransposeMatrixVectorWithCOO(double[,] matrixA, double[] vectorX, string storingMethodMajor = "Row")
        {
            int numRowsMatrixA = matrixA.GetLength(0);
            int numColsMatrixA = matrixA.GetLength(1);
            int numVectorX = vectorX.Length;

            if (numRowsMatrixA != numVectorX)
            {
                throw new Exception("Cannot Multiply Matrix A with Vector X with these dimensions!");
            }

            // Row Major COO
            if (storingMethodMajor == "Row")
            {
                (double[] valuesArray, int[] rowsArray, int[] colsArray) = Matrices.StoringMatrices.SaveMatrixCOO(matrixA, storingMethodMajor);

                double[] resultVector = new double[numColsMatrixA];

                for (int i = 0; i < valuesArray.Length; i++)
                {
                    double value = valuesArray[i];
                    int rowIndex = rowsArray[i];
                    int colIndex = colsArray[i];

                    resultVector[colIndex] = resultVector[colIndex] + value * vectorX[rowIndex];
                }

/*                Console.WriteLine("");
                Console.WriteLine("Result Vector = {" + string.Join(", ", resultVector) + "}");*/

                return resultVector;
            }
            // Column Major COO
            else
            {
                (double[] valuesArray, int[] colsArray, int[] rowsArray) = Matrices.StoringMatrices.SaveMatrixCOO(matrixA, storingMethodMajor);

                double[] resultVector = new double[numColsMatrixA];

                for (int i = 0; i < valuesArray.Length; i++)
                {
                    double value = valuesArray[i];
                    int colIndex = colsArray[i];
                    int rowIndex = rowsArray[i];

                    resultVector[colIndex] = resultVector[colIndex] + value * vectorX[rowIndex];
                }

/*                Console.WriteLine("");
                Console.WriteLine("Result Vector = {" + string.Join(", ", resultVector) + "}");*/

                return resultVector;
            }

        }

        /// <summary>
        /// Test for solving a linear system A * x = b with the SolveLinearSystemMatrixVectorWithCOO function.
        /// </summary>
        /// <remarks>Uses a specific example and compares computed result with the expected result.</remarks>
        public static void TestSolveLinearSystemMatrixVectorWithCOO()
        {
            double[,] A =
            {
                { 1, 2, 3, 4,
                },
                { 1, 1, 1, 1,
                },
            };

            double[] x = { 1, 2, 3, 4 };

            double[] expected = { 30, 10 };

            double[] computed = SolveLinearSystemMatrixVectorWithCOO(A, x, "Row");

            Matrices.VectorMethods.CheckVector(expected, computed);
            Matrices.VectorMethods.PrintVector(expected);
            Matrices.VectorMethods.PrintVector(computed);

        }

        /// <summary>
        /// Test for solving a linear system A^T * x = b with the SolveLinearSystemTransposeMatrixVectorWithCOO function.
        /// </summary>
        /// <remarks>Uses a specific example and compares computed result with the expected result.</remarks>
        public static void TestSolveLinearSystemTransposeMatrixVectorWithCOO()
        {
            double[,] A =
            {
                { 1, 2, 3, 4,
                },
                { 1, 1, 1, 1,
                },
            };

            double[] x = { 1, 1 };

            double[] expected = { 2, 3, 4, 5 };

            double[] computed = SolveLinearSystemTransposeMatrixVectorWithCOO(A, x, "Row");

            Matrices.VectorMethods.CheckVector(expected, computed);
            Matrices.VectorMethods.PrintVector(expected);
            Matrices.VectorMethods.PrintVector(computed);

        }
    }
}
