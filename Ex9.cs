using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex9
    {
        /// <summary>
        /// Returns the value of an element in a banded matrix.
        /// </summary>
        /// <param name="bandedMatrix">The array representing the upper triangular elements of the banded matrix.</param>
        /// <param name="bandedMatrixDimension">The dimension of the square banded matrix.</param>
        /// <param name="matrixBandwidth">The bandwidth of the banded matrix.</param>
        /// <param name="i">The row index of the element.</param>
        /// <param name="j">The column index of the element.</param>
        /// <returns>The value of the element at the specified indices (i, j).</returns>
        /// <exception cref="Exception">Thrown when the provided indices are invalid for the banded matrix.</exception>
        public static double ReturnValueFromBandedMatrix(double[] bandedMatrix, int bandedMatrixDimension, int matrixBandwidth, int i, int j)
        {
            if (i < 0 || i >= bandedMatrixDimension || j < 0 || j >= bandedMatrixDimension)
            {
                throw new Exception("Invalid indices i, j for the Banded Matrix.");
            }

            // Check if the indices are within the bandwidth
            if (Math.Abs(i - j) > matrixBandwidth)
            {
                return 0.0; // Return 0 for elements outside the bandwidth
            }

            // Switch i and j if is on lower triangle
            if (i > j)
            {
                int temp = i;
                i = j;
                j = temp;
            }

            // Calculate the index for the upper triangular elements
            int index = i + (j - i) * bandedMatrixDimension;

            Console.WriteLine("Result Value:" + bandedMatrix[index]);

            return bandedMatrix[index];
        }

        public static void TestReturnValueFromBandedMatrix()
        {
            double[] banded_matrix = { 1.0, 4.0, 3.0, 5.0, 0, 5, 1, 0, 2 };
            int dimension = 4;
            int bandwidth = 2;
            int row = 1;
            int col = 2;

            double expected = 5.0 ;

            double computed = ReturnValueFromBandedMatrix(banded_matrix, dimension, bandwidth, row, col);

            if (computed != expected) 
            {
                Console.WriteLine("Expected Value is not equal with the Computed!");
            }

            Console.WriteLine("Expected Value:" + expected);
            Console.WriteLine("Computed Value:" + computed);
        }

    }
}
