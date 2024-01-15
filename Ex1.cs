using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MscNumericalLinearAlgebra.Matrices;

namespace MscNumericalLinearAlgebra.ExcerciseSeries2
{
    public static class Ex1
    {
        public static double[] MultiplyMatrixAVectorX(double[,] matrixA, double[] vectorX)
        {
            int lenRowMatrixA = matrixA.GetLength(0);
            int lenColMatrixA = matrixA.GetLength(1);
            int lenVectorX = vectorX.GetLength(0);

            if (lenColMatrixA != lenVectorX)
            {
                throw new Exception("The sizes of Matrix A and vector B are not equal!");
            }

            double[] resultVector = new double[lenRowMatrixA];

            for (int i = 0; i < lenRowMatrixA; i++)
            {
                for (int j = 0; j < lenColMatrixA; j++)
                {
                    resultVector[i] = resultVector[i] + matrixA[i,j] * vectorX[j];
                } 
            }

            Console.WriteLine("");
            Console.WriteLine("b = {" + string.Join(", ", resultVector) + "}");

            return resultVector;
        }

        public static void TestMultiplyMatrixAVectorX()
        {
            double[,] matrixA =
                {
                { 1.0, 0.0, 3.0 },
                { 0.0, 4.0, 0.0 }
            };

            double[] vectorX = { 1.0, 1.0, 1.0 };

            double[] expected = { 4.0, 4.0 };
            double[] computed = MultiplyMatrixAVectorX(matrixA, vectorX);

            VectorMethods.CheckVector(computed, expected);
                   
        }
    }
}
