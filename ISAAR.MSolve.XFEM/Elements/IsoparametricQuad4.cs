﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Interpolation.ShapeFunctions;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Elements
{
    class IsoparametricQuad4: IFiniteElement2D
    {
        private const int DOFS_COUNT = 8;

        private readonly IReadOnlyList<Node2D> nodes;
        private readonly IsoparametricInterpolation2D interpolation;

        public IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> MaterialsOfGaussPoints { get; }

        public static IsoparametricQuad4 CreateHomogeneous(Node2D[] nodes, IFiniteElementMaterial2D material)
        {
            var nodesCopy = new Node2D[nodes.Length];
            nodes.CopyTo(nodesCopy, 0);

            var gpToMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in GaussQuadrature2D.Order2x2.GenerateIntegrationPoints()) // TODO: remove the integration logic from the element class
            {
                gpToMaterials[point] = material.Clone();
            }

            return new IsoparametricQuad4(nodesCopy, gpToMaterials);
        }

        /// <summary>
        /// Parameters passed to this constructor will not be copied. The static factory methods are more robust.
        /// </summary>
        /// <param name="nodes">Is not deep copied.</param>
        /// <param name="materialsOfGaussPoints">Is not deep copied</param>
        public IsoparametricQuad4(IReadOnlyList<Node2D> nodes, 
            IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> materialsOfGaussPoints)
        {
            // TODO: Add checks here: order of nodes and suitability of gauss points. 
            // Or add checks in the callers (in this class and in enriched elements)?
            this.nodes = nodes;
            this.interpolation = new IsoparametricInterpolation2D(this.nodes, NaturalShapeFunctions2D.Quad4);
            this.MaterialsOfGaussPoints = materialsOfGaussPoints;
        }

        public SymmetricMatrix2D<double> BuildStiffnessMatrix()
        {
            var stiffness = new SymmetricMatrix2D<double>(DOFS_COUNT);
            foreach (var entry in MaterialsOfGaussPoints)
            {
                GaussPoint2D gaussPoint = entry.Key;
                IFiniteElementMaterial2D material = entry.Value;

                // Calculate the necessary quantities for the integration
                ShapeFunctionDerivatives2D interpolationDerivatives = 
                    this.interpolation.EvaluateDerivativesAt(gaussPoint.X, gaussPoint.Y);
                Matrix2D<double> deformation = CalculateDeformationMatrix(interpolationDerivatives);
                Matrix2D<double> constitutive = material.CalculateConstitutiveMatrix();
                double thickness = material.Thickness;

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D<double> partial = (deformation.Transpose() * constitutive) * deformation; // Perhaps this could be done in a faster way taking advantage of symmetry.
                partial.Scale(material.Thickness * interpolationDerivatives.Jacobian.Determinant * gaussPoint.Weight);
                Debug.Assert(partial.Rows == DOFS_COUNT);
                Debug.Assert(partial.Columns == DOFS_COUNT);
                MatrixUtilities.AddPartialToSymmetricTotalMatrix(partial, stiffness);
            }
            return stiffness;
        }

        /// <summary>
        /// Calculates the deformation matrix B. Dimensions = 3x8.
        /// B is a linear transformation FROM the nodal values of the displacement field TO the the derivatives of
        /// the displacement field in respect to the cartesian axes (i.e. the stresses): {dU/dX} = [B] * {d} => 
        /// {u,x v,y u,y, v,x} = [... Bk ...] * {u1 v1 u2 v2 u3 v3 u4 v4}, where k = 1, ... nodesCount is a node and
        /// Bk = [dNk/dx 0; 0 dNk/dY; dNk/dy dNk/dx] (3x2)
        /// </summary>
        /// <param name="shapeFunctionDerivatives">The shape function derivatives calculated at a specific 
        ///     integration point</param>
        /// <returns></returns>
        public Matrix2D<double> CalculateDeformationMatrix(ShapeFunctionDerivatives2D shapeFunctionDerivatives)
        {
            var deformationMatrix = new Matrix2D<double>(3, DOFS_COUNT);
            for (int nodeIndex = 0; nodeIndex < 4; ++nodeIndex)
            {
                int col1 = 2 * nodeIndex;
                int col2 = 2 * nodeIndex + 1;
                Tuple<double, double> dNdX = shapeFunctionDerivatives.CartesianDerivativesOfNode(nodeIndex);

                deformationMatrix[0, col1] = dNdX.Item1;
                deformationMatrix[1, col2] = dNdX.Item2;
                deformationMatrix[2, col1] = dNdX.Item2;
                deformationMatrix[2, col2] = dNdX.Item1;
            }
            return deformationMatrix;
        }
    }
}
