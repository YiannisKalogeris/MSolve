using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Interpolation.Inverse
{
	/// <summary>
	/// Inverse mapping of the isoparametric interpolation of a hexahedral finite element with 8 nodes. Since the original
	/// mapping is linear, there are analytic formulas, which are presented in 
	/// "The inverse mapping and distortion measures for 8-node hexahedral isoparametric elements", K. -Y. YuanY. -S. HuangH. -T. YangT. H. H. Pian, 1994
	/// https://link.springer.com/article/10.1007/BF00350284
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InverseInterpolationHexa8 : IInverseInterpolation3D
    {
        private readonly IReadOnlyList<Node_v2> nodes;
        private readonly InterpolationHexa8 directMapping;

        public InverseInterpolationHexa8(IReadOnlyList<Node_v2> nodes, InterpolationHexa8 directMapping)
	    {
            this.nodes = nodes;
            this.directMapping = directMapping;
	    }

		public NaturalPoint3D TransformPointCartesianToNatural(CartesianPoint3D point)
		{
            // Adapted from the old Hexa8. Uses an iterative (NR I think) method.

            double[] mins = new double[] {nodes[0].X,nodes[0].Y,nodes[0].Z };
            double[] maxes = new double[] {nodes[0].X,nodes[0].Y,nodes[0].Z };
            for (int i = 0; i <nodes.Count; i++)
            {
                mins[0] = mins[0] > nodes[i].X ? nodes[i].X : mins[0];
                mins[1] = mins[1] > nodes[i].Y ? nodes[i].Y : mins[1];
                mins[2] = mins[2] > nodes[i].Z ? nodes[i].Z : mins[2];
                maxes[0] = maxes[0] < nodes[i].X ? nodes[i].X : maxes[0];
                maxes[1] = maxes[1] < nodes[i].Y ? nodes[i].Y : maxes[1];
                maxes[2] = maxes[2] < nodes[i].Z ? nodes[i].Z : maxes[2];
            }
            //return new double[] { (node.X - mins[0]) / ((maxes[0] - mins[0]) / 2) - 1,
            //    (node.Y - mins[1]) / ((maxes[1] - mins[1]) / 2) - 1,
            //    (node.Z - mins[2]) / ((maxes[2] - mins[2]) / 2) - 1 };

            bool maybeInsideElement = point.X <= maxes[0] && point.X >= mins[0] &&
                point.Y <= maxes[1] && point.Y >= mins[1] &&
                point.Z <= maxes[2] && point.Z >= mins[2];
            if (maybeInsideElement == false) return null;

            const int jacobianSize = 3;
            const int maxIterations = 1000;
            const double tolerance = 1e-10;
            int iterations = 0;
            double deltaNaturalCoordinatesNormSquare = 100;
            double[] naturalCoordinates = new double[] { 0, 0, 0 };
            const double toleranceSquare = tolerance * tolerance;

            while (deltaNaturalCoordinatesNormSquare > toleranceSquare && iterations < maxIterations)
            {
                iterations++;
                var shapeFunctions = directMapping.EvaluateFunctionsAt(new NaturalPoint3D(naturalCoordinates));
                //var shapeFunctions = CalcH8Shape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
                double[] coordinateDifferences = new double[] { 0, 0, 0 };
                for (int i = 0; i < shapeFunctions.Length; i++)
                {
                    coordinateDifferences[0] += shapeFunctions[i] *nodes[i].X;
                    coordinateDifferences[1] += shapeFunctions[i] *nodes[i].Y;
                    coordinateDifferences[2] += shapeFunctions[i] *nodes[i].Z;
                }
                coordinateDifferences[0] = point.X - coordinateDifferences[0];
                coordinateDifferences[1] = point.Y - coordinateDifferences[1];
                coordinateDifferences[2] = point.Z - coordinateDifferences[2];

                //double[,] faXYZ = GetCoordinatesTranspose(element);
                //double[] nablaShapeFunctions = CalcH8NablaShape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
                //var inverseJacobian = CalcH8JDetJ(faXYZ, nablaShapeFunctions).Item2;

                Matrix naturalShapeFunctionDerivatives = directMapping.EvaluateNaturalGradientsAt(new NaturalPoint3D(naturalCoordinates));
                var jacobian = new IsoparametricJacobian3D(nodes, naturalShapeFunctionDerivatives);
                Matrix inverseJacobian = jacobian.InverseMatrix;

                double[] deltaNaturalCoordinates = new double[] { 0, 0, 0 };
                for (int i = 0; i < jacobianSize; i++)
                    for (int j = 0; j < jacobianSize; j++)
                        deltaNaturalCoordinates[i] += inverseJacobian[j, i] * coordinateDifferences[j];
                for (int i = 0; i < 3; i++)
                    naturalCoordinates[i] += deltaNaturalCoordinates[i];

                deltaNaturalCoordinatesNormSquare = 0;
                for (int i = 0; i < 3; i++)
                    deltaNaturalCoordinatesNormSquare += deltaNaturalCoordinates[i] * deltaNaturalCoordinates[i];
                //deltaNaturalCoordinatesNormSquare = Math.Sqrt(deltaNaturalCoordinatesNormSquare);
            }

            return new NaturalPoint3D(naturalCoordinates.Count(x => Math.Abs(x) - 1.0 > tolerance) > 0 ? new double[0] : naturalCoordinates);
        }
	}
}
