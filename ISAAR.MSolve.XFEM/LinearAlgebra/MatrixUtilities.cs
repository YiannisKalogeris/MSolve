﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.LinearAlgebra
{
    /// <summary>
    /// TODO: Some of these could be extension methods, but effort should be spent into implementing them in the  
    /// Matrix classes directly. If they indeed only need IMatrix, then they should be implemented as extensions.
    /// </summary>
    static class MatrixUtilities
    {
        //public static void AddPartialToSymmetricTotalMatrix(Matrix partialMatrix,
        //    SymmetricMatrix totalMatrix)
        //{
        //    Debug.Assert(partialMatrix.NumRows == totalMatrix.NumRows, "Non matching rows");
        //    Debug.Assert(partialMatrix.NumColumns == totalMatrix.NumColumns, "Non matching columns");
        //    for (int row = 0; row < totalMatrix.NumRows; ++row)
        //    {
        //        for (int col = row; col < totalMatrix.NumColumns; ++col)
        //        {
        //            totalMatrix[row, col] += partialMatrix[row, col];
        //        }
        //    }
        //}

        public static void PrintUpperTriangleDense(IMatrixView matrix)
        {
            StringBuilder builder = new StringBuilder();
            for (int row = 0; row < matrix.NumRows; ++row)
            {
                for (int col = 0; col < matrix.NumColumns; ++col)
                {
                    if (col < row) builder.Append(0.0);
                    else builder.Append(matrix[row, col]);
                    builder.Append(' ');
                }
                builder.Append("\n");
            }
            builder.Append("\n");
            Console.Write(builder.ToString());
        }

        /// <summary>
        /// The output is a dense matrix
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="oldToNewIndices"></param>
        /// <returns></returns>
        public static Matrix Reorder(IIndexable2D matrix, int[] oldToNewIndices) // TODO: add this to Matrix or as an extension method
        {
            if (matrix.NumRows != matrix.NumColumns) throw new ArgumentException("The matrix must be square");
            if (matrix.NumRows != oldToNewIndices.Length) throw new NonMatchingDimensionsException(
                "Mismatch in the dimensions of the matrix and the permutation array");

            var newMatrix = Matrix.CreateZero(matrix.NumRows, matrix.NumRows);
            for (int row = 0; row < matrix.NumRows; ++row)
            {
                int newRow = oldToNewIndices[row];
                for (int col = 0; col < matrix.NumColumns; ++col)
                {
                    newMatrix[newRow, oldToNewIndices[col]] = matrix[row, col];
                }
            }
            return newMatrix;
        }

        public static double[] Substract(double[] vector1, double[] vector2)
        {
            if (vector1.Length != vector2.Length) throw new NonMatchingDimensionsException(
                "Vector 1 length  = " + vector1.Length + ", but vector 2 length = " + vector2.Length);

            double[] result = new double[vector1.Length];
            for (int i = 0; i < vector1.Length; ++i) result[i] = vector1[i] - vector2[i];
            return result;
        }
    }
}
