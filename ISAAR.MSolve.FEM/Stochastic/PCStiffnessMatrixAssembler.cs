using ISAAR.MSolve.PreProcessor.Stochastic;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Stochastic
{
    public static class PCStiffnessMatrixAssembler
    {
        public static double[,] PCGlobalStiffnessMatrix(int KLterms, int PCorder, double[,] Kmean,double[][,] Kstd, bool isGaussian ) {
            int PCterms = PolynomialChaosCoefficientsProvider.HermitePC(KLterms,PCorder).Item4;
            double[] PsiSquareNorm = PolynomialChaosCoefficientsProvider.HermitePC(KLterms, PCorder).Item3;
            double[][,] PCCoefficientsMatrices = PolynomialChaosCoefficientsProvider.PCCoefficientsCalculator(KLterms, PCorder, isGaussian);
            double[,] StiffnessMatrixSSFEM =new double[PCterms*Kmean.GetLength(0), PCterms * Kmean.GetLength(0)];

            double[][,] Ki = new double[KLterms + 1][,];
            for (int i = 0; i <= KLterms; i++) {
                if (i == 0)
                {
                    Ki[i] = Kmean;
                }
                else Ki[i] = Kstd[i + 1];
            }

            double[,][,] KK = new double[PCterms, PCterms][,];
            for (int j = 0; j < PCterms; j++)
            {
                for (int k = 0; k < PCterms; k++)
                {
                    KK[j, k] = new double[Kmean.GetLength(0), Kmean.GetLength(0)];
                    for (int i = 0; i <= KLterms; i++)
                    {
                        if (PCCoefficientsMatrices[i][j, k] != 0)
                        {
                            KK[j, k] = PCStiffnessMatrixAssembler.MatrixMatrixAddition(KK[j, k], PCStiffnessMatrixAssembler.MatrixByConstantMult(Ki[i], PCCoefficientsMatrices[i][j, k]));
                        }
                    }
                }
            }

            for (int i = 0; i < PCterms; i++)
            {
                for (int j = 0; j < PCterms; j++)
                {
                    for (int k1 = 0; k1 < Kmean.GetLength(0); k1++)
                    {
                        for (int k2 = 0; k2 < Kmean.GetLength(0); k2++)
                        {
                            StiffnessMatrixSSFEM[i * k1, j * k2] = KK[i,j][k1, k2];
                        }
                    }
                 }
            }

         return StiffnessMatrixSSFEM;
        }

        private static double[,] MatrixByConstantMult(double[,] Matrix, double c)
        {
            double[,] newMatrix = new double[Matrix.GetLength(0), Matrix.GetLength(1)];
            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(1); j++)
                {
                    newMatrix[i, j] = Matrix[i, j] * c;
                }
            }
            return newMatrix;
        }

        private static double[,] MatrixMatrixAddition(double[,] Matrix1, double[,] Matrix2)
        {
            double[,] newMatrix = new double[Matrix1.GetLength(0), Matrix1.GetLength(1)];
            for (int i = 0; i < Matrix1.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix1.GetLength(1); j++)
                {
                    newMatrix[i, j] = Matrix1[i, j] + Matrix2[i, j];
                }
            }
            return newMatrix;
        }
    }
}
