using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.PreProcessor.Stochastic
{
    public static class PolynomialChaosCoefficientsProvider
    {

        public static double[][,] PCCoefficientsCalculator(int KLterms, int PCorder, bool isGaussian)
        {
            //if (KLterms < 1) throw new ArgumentException("KLterms cannot be less than 1");
            //if (PCorder < 1) throw new ArgumentException("p cannot be less than 1");

            double[,] ALPHA = PolynomialChaosCoefficientsProvider.HermitePC(KLterms, PCorder).Item1;
            double[][] Psi = PolynomialChaosCoefficientsProvider.HermitePC(KLterms, PCorder).Item2;
            double[] PsiSquareNorm = PolynomialChaosCoefficientsProvider.HermitePC(KLterms, PCorder).Item3;
            int PCterms= PolynomialChaosCoefficientsProvider.HermitePC(KLterms, PCorder).Item4;

            double[][,] PCCoefficientsMatrices = new double[KLterms + 1][,];
            for (int i = 0; i <= KLterms; i++)
            {
                PCCoefficientsMatrices[i] = new double[PCterms, PCterms];
            }
            for (int i = 0; i < PCterms; i++)
            {
                PCCoefficientsMatrices[0][i, i] = PsiSquareNorm[i];
            }
            for (int j = 0; j < PCterms; j++)
            {
                double[] ALPHAp = new double[ALPHA.GetLength(1)];
                for (int j1 = 0; j1 < ALPHA.GetLength(1); j1++)
                {
                    ALPHAp[j1] = ALPHA[j, j1];
                }
                for (int k = 0; k < PCterms; k++)
                {
                    double[] ALPHAq = new double[ALPHA.GetLength(1)];
                    for (int k1 = 0; k1 < ALPHA.GetLength(1); k1++)
                    {
                        ALPHAq[k1] = ALPHA[k, k1];
                    }
                    for (int i = 0; i < KLterms; i++)
                    {
                       int[] idx = PolynomialChaosCoefficientsProvider.IDX(i, KLterms);
                        bool Condition = true;
                        for (int k2 = 0; k2 < idx.Length; k2++)
                        {
                            if (ALPHAp[idx[k2]] != ALPHAq[idx[k2]])
                            {
                                Condition = false;
                            }
                        }
                        if (Condition == true)
                        {
                            int dummy1 = 0;
                            int dummy2 = 0;
                            if (ALPHAp[i] - 1 == ALPHAq[i])
                            {
                                dummy1 = 1;
                            }
                            if (ALPHAp[i] == ALPHAq[i] - 1)
                            {
                                dummy2 = 1;
                            }
                            PCCoefficientsMatrices[i + 1][j, k] = (PolynomialChaosCoefficientsProvider.Factorial((int)ALPHAp[i]) * dummy1
                                + PolynomialChaosCoefficientsProvider.Factorial((int)ALPHAq[i]) * dummy2) * PsiSquareNorm[j] / PolynomialChaosCoefficientsProvider.Factorial((int)ALPHAp[i]);
                        }
                    }

              }
         }

         return PCCoefficientsMatrices;
        }


        //auxillary methods required
        public static int[] IDX(int i, int KLterms)
        {
            int[] a = new int[i];
            int L1 = a.Length;
            int[] b = new int[KLterms - i-1];
            int L2 = b.Length;
            int L = L1 + L2;

            int[] idx = new int[L];
            if (L1 != 0 && L2 !=0)
            {
                for (int j = 0; j < i-1; j++)
                {
                    a[j] = j;
                }
                for (int j = i+1; j < KLterms; j++)
                {
                    b[j-(i+1)] = j;
                }

             Array.Copy(a, idx, L1);
             Array.Copy(b, 0, idx, L1, L2);
            }
            if (L1 == 0 && L2 != 0)
            {
                for (int j = i + 1; j < KLterms; j++)
                {
                    b[j - (i + 1)] = j;
                }
                Array.Copy(b, idx, L2);
            }
            if (L1 != 0 && L2 == 0)
            {
                for (int j = 0; j < i - 1; j++)
                {
                    a[j] = j;
                }
                Array.Copy(a, idx, L1);
            }

            return idx;
            
        }

        private static Tuple<double[,], double[][], double[], int> HermitePC(int KLterms, int PCorder)
        //private static (double[,] ALPHA, double[][] Psi, double[] PsiSquareNorm, int PCterms) HermitePC(int KLterms, int PCorder)
        {
            int PCterms = 0;
            for (int k = 0; k <= PCorder; k++)
            {
                PCterms += (int)(PolynomialChaosCoefficientsProvider.Factorial(KLterms + k - 1) / (PolynomialChaosCoefficientsProvider.Factorial(k) * PolynomialChaosCoefficientsProvider.Factorial(KLterms - 1)));
            }

            double[][] He = new double[PCorder + 2][];
            He[0] = new double[1] { 1 };
            He[1] = new double[2] { 1, 0 };

            for (int i = 2; i <= PCorder + 1; i++)
            {
                He[i] = new double[i + 1];
                double[] temp1 = new double[He[i - 1].Length + 1];
                for (int j = 0; j < He[i - 1].Length; j++)
                {
                    temp1[j] = He[i - 1][j];
                }
                double[] temp2 = new double[He[i - 2].Length + 2];
                for (int j = 2; j < He[i - 2].Length + 2; j++)
                {
                    temp2[j] = He[i - 2][j - 2] * (i - 1);
                }
                for (int i1 = 0; i1 < temp2.Length; i1++)
                {
                    He[i][i1] = temp1[i1] - temp2[i1];
                }
            }

            double[,][] H = new double[PCorder + 1, KLterms][];
            for (int j = 0; j <= KLterms - 1; j++)
            {
                for (int i = 0; i <= PCorder; i++)
                {
                    H[i, j] = new double[He[i].Length];
                    for (int k = 0; k <= He[i].Length - 1; k++)
                    {
                        H[i, j][k] = He[i][k];
                    }
                }
            }

            double[][] Psi = new double[PCterms][];
            double[,] ALPHA = PolynomialChaosCoefficientsProvider.MultiIndex(KLterms, PCorder);
            for (int i = 0; i < PCterms; i++)
            {
                double[][] mult = new double[KLterms+1][];
                mult[0] = new double[1] { 1 };
                for (int j = 0; j < KLterms; j++)
                {
                    mult[j+1] = PolynomialChaosCoefficientsProvider.Convolute(mult[j], H[(int)(ALPHA[i, j] ), j]);
                }
                Psi[i] = mult[KLterms];
            }

            double[,] FactorialALPHA = new double[ALPHA.GetLength(0), ALPHA.GetLength(1)];
            for (int i = 0; i <= ALPHA.GetLength(0) - 1; i++)
            {
                for (int j = 0; j <= ALPHA.GetLength(1) - 1; j++)
                {
                    FactorialALPHA[i, j] = Factorial((int)ALPHA[i, j]);
                }
            }

            double[] PsiSquareNorm = new double[FactorialALPHA.GetLength(0)];
            for (int i = 0; i < FactorialALPHA.GetLength(0); i++)
            {
                double temp = 1;
                for (int j = 0; j < FactorialALPHA.GetLength(1); j++)
                {
                    temp*= FactorialALPHA[i, j];
                }
                PsiSquareNorm[i] = temp;
            }


        return new Tuple<double[,], double[][], double[], int>(ALPHA, Psi, PsiSquareNorm,PCterms);
        }

        public static double Factorial(int f)
        {
            double result = 1;
            for (int i = 2; i <= f; i++)
            {
                result *= i;
            }
            return result;
        }

        public static double[,] MultiIndex(int KLterms, int PCorder)
        {
            double[][,] alpha = new double[PCorder + 1][,];
            double[] dimAlpha = new double[PCorder + 1];
            alpha[0] = new double[1, KLterms];
            dimAlpha[0] = alpha[0].GetLength(0);
            if (KLterms == 1)
            {
                for (int i = 1; i <= PCorder; i++) //?????????? -1??
                {
                    alpha[i] = new double[1, 1];
                    alpha[i][0, 0] = i;
                    dimAlpha[i] = alpha[i].GetLength(0);
                }
            }
            if (KLterms > 1)
            {
                for (int i = 1; i <= PCorder; i++)
                {
                    double[,] s = NChooseK(KLterms - 1, i + KLterms - 1);
                    double[] s1 = new double[s.GetLength(0)];
                    double[] s2 = new double[s.GetLength(0)];
                    for (int j = 0; j <= s2.Length - 1; j++)
                    {
                        s2[j] = KLterms + i;
                    }
                    double[,] S = new double[s1.Length, s.GetLength(1) + 2];
                    for (int j = 0; j <= s2.Length - 1; j++)
                    {
                        for (int k = 0; k <= s.GetLength(1)+1; k++)
                        {
                            if (k == 0)
                            {
                                S[j, k] = s1[j];
                            }
                            else if (k == s.GetLength(1)+1 )
                            {
                                S[j, k] = s2[j];
                            }
                            else
                                S[j, k] = s[j ,k-1];
                        }

                    }
                    double[,] DS = Diff(S);
                    double[,] flippedDS = FlipUD(DS);
                    alpha[i] = new double[flippedDS.GetLength(0), flippedDS.GetLength(1)];
                    for (int j = 0; j <= flippedDS.GetLength(0) - 1; j++)
                    {
                        for (int k = 0; k <= flippedDS.GetLength(1) - 1; k++)
                        {
                            alpha[i][j, k] = flippedDS[j, k]-1;
                        }
                    }
                    dimAlpha[i] = alpha[i].GetLength(0);

                }
            }
            // convert cell alpha into matrix ALPHA
            int totalDimAlpha = 0;
            for (int i = 0; i <= PCorder; i++)
            {
                totalDimAlpha +=(int) dimAlpha[i];
            }
            int count = 0;
            double[,] ALPHA = new double[totalDimAlpha, KLterms];
            for (int i = 0; i <= PCorder; i++) {
                for (int j = 0; j <= dimAlpha[i]-1; j++)
                {
                    for (int k = 0; k <= KLterms-1; k++)
                    {
                        ALPHA[count,k]= alpha[i][j,k];
                    }
                    count += 1;
                }
            }
                return ALPHA;
        }

        private static double[,] NChooseK(int k, int n)
        {
            int NoOfCombinations = 0;
            NoOfCombinations = (int)(PolynomialChaosCoefficientsProvider.Factorial(n) / (PolynomialChaosCoefficientsProvider.Factorial(k) * PolynomialChaosCoefficientsProvider.Factorial(n - k)));

            var Combinations = new double[NoOfCombinations, k];
            var result = new int[k];
            var stack = new Stack<int>();
            stack.Push(1);
            int j = 0;

            while (stack.Count > 0)
            {
                var index = stack.Count - 1;
                var value = stack.Pop();

                while (value <= n)
                {
                    result[index++] = value++;
                    stack.Push(value);
                    if (index == k)
                    {
                        for (int i = 0; i <= result.Length-1 ; i++)
                        {
                            Combinations[j, i] = result[i];
                        }
                        j += 1;
                        break;
                        
                    }
                    
                }
                
            }
            return Combinations;
        }

        private static double[,] Diff(double[,] M)
        {
            double[,] DM = new double[M.GetLength(0) , M.GetLength(1)-1];
            for (int i = 0; i <= M.GetLength(0) - 1; i++)
            {
                for (int j = 0; j <= M.GetLength(1) - 2; j++)
                {
                    DM[i, j] = M[i , j+1] - M[i, j];
                }
            }
            return DM;
        }

        private static double[,] FlipUD(double[,] M)
        {
            double[,] FlippedM = new double[M.GetLength(0), M.GetLength(1)];
            for (int i = 0; i <= M.GetLength(0) - 1; i++)
            {
                for (int j = 0; j <= M.GetLength(1) - 1; j++)
                {
                    FlippedM[i, j] = M[M.GetLength(0) - 1 - i, j];
                }
            }
            return FlippedM;
        }

        private static double[] Convolute(double[] v1, double[] v2)
        {
            int vecSize1 = v1.Length;
            int vecSize2 = v2.Length;
            int size = vecSize1 + vecSize2 - 1;
            var result = new double[size];
            var v1Ext = new double[size];
            var v2Ext = new double[size];

            var v2Reverse = v2.Reverse().ToArray();
            Array.Clear(v2Ext, 0, size);
            Array.Copy(v2Reverse, 0, v2Ext, size - vecSize2, vecSize2);
            int index = size - 1;
            for (int i = 0; i < size - vecSize1 + 1; i++)
            {
                Array.Clear(v1Ext, 0, size);
                Array.Copy(v1, 0, v1Ext, i, vecSize1);
                for (int j = 0; j < size; j++)
                    result[index] += v1Ext[j] * v2Ext[j];
                index--;
            }

            for (int i = size - vecSize2 - 1; i >= 0; i--)
            {
                Array.Clear(v2Ext, 0, size);
                Array.Copy(v2Reverse, 0, v2Ext, i, vecSize2);
                for (int j = 0; j < size; j++)
                    result[index] += v1Ext[j] * v2Ext[j];
                index--;
            }
            return result;
        }
    }
}
