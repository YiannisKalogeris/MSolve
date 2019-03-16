using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: How can we ensure that the model has the correct shape and discretization?
namespace ISAAR.MSolve.Analyzers.Multiscale
{
    /// <summary>
    /// This works only for cubic (3D) RVEs, centered at (0,0,0), with 1 thermal dof per node and linear boundary conditions.
    /// </summary>
    public class ThermalCubeRve : IReferenceVolumeElement
    {
        private const int numDimensions = 3;

        private readonly Vector3 macroscopicTemperatureGradient;
        private readonly double xMin, yMin, zMin, xMax, yMax, zMax;
        private readonly IEnumerable<INode> nonEmbeddedNodes;
        private readonly HashSet<INode> nodesXmin, nodesXmax, nodesYmin, nodesYmax, nodesZmin, nodesZmax;

        /// <summary>
        /// </summary>
        /// <param name="model"></param>
        /// <param name="xyzMinCoords"></param>
        /// <param name="xyzMaxCoords"></param>
        /// <param name="macroscopicTemperatureGradient"></param>
        /// <param name="meshTolerance">The default is 1E-10 * min(|xMax-xMin|, |yMax-yMin|)</param>
        public ThermalCubeRve(IEnumerable<INode> nonEmbeddedNodes, Vector3 xyzMinCoords, Vector3 xyzMaxCoords,
            Vector3 macroscopicTemperatureGradient, double meshTolerance)
        {
            this.nonEmbeddedNodes = nonEmbeddedNodes;
            this.xMin = xyzMinCoords[0];
            this.yMin = xyzMinCoords[1];
            this.zMin = xyzMinCoords[2];
            this.xMax = xyzMaxCoords[0];
            this.yMax = xyzMaxCoords[1];
            this.zMax = xyzMaxCoords[2];
            this.macroscopicTemperatureGradient = macroscopicTemperatureGradient;

            // Find the nodes of each edge
            nodesXmin = new HashSet<INode>();
            nodesXmax = new HashSet<INode>();
            nodesYmin = new HashSet<INode>();
            nodesYmax = new HashSet<INode>();
            nodesZmin = new HashSet<INode>();
            nodesZmax = new HashSet<INode>();
            foreach (INode node in nonEmbeddedNodes)
            {
                // Some edges are prioritized for corner nodes. //TODO: should the corner nodes be handled differently?
                if (Math.Abs(node.Z - zMax) <= meshTolerance) nodesZmax.Add(node);
                else if (Math.Abs(node.Y - yMax) <= meshTolerance) nodesYmax.Add(node);
                else if (Math.Abs(node.X - xMax) <= meshTolerance) nodesXmax.Add(node);
                else if (Math.Abs(node.Z - zMin) <= meshTolerance) nodesZmin.Add(node);
                else if (Math.Abs(node.Y - yMin) <= meshTolerance) nodesYmin.Add(node);
                else if (Math.Abs(node.X - xMin) <= meshTolerance) nodesXmin.Add(node);
            }
        }

        public ThermalCubeRve(IEnumerable<INode> nonEmbeddedNodes, Vector3 bottomLeftCoords, Vector3 topRightCoords,
            Vector3 macroscopicTemperatureGradient) : 
            this(nonEmbeddedNodes, bottomLeftCoords, topRightCoords, macroscopicTemperatureGradient, 
                1E-10 * topRightCoords.Subtract(bottomLeftCoords).MinAbsolute())
        {
        }

        public void ApplyBoundaryConditions()
        {
            double dTdx = macroscopicTemperatureGradient[0];
            double dTdy = macroscopicTemperatureGradient[1];
            double dTdz = macroscopicTemperatureGradient[2];
            foreach (var node in nodesXmin) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 0.0 });
            foreach (var node in nodesYmin) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 0.0 });
            foreach (var node in nodesZmin) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 0.0 });
            foreach (var node in nodesXmax) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = dTdx });
            foreach (var node in nodesYmax) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = dTdy });
            foreach (var node in nodesZmax) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = dTdz });
        }

        public double CalculateRveVolume() => (xMax - xMin) * (yMax - yMin) * (zMax - zMin);

        public IMatrixView CalculateKinematicRelationsMatrix(ISubdomain_v2 subdomain)
        {
            ISubdomainConstrainedDofOrdering constrainedDofOrdering = subdomain.ConstrainedDofOrdering;
            var kinematicRelations = Matrix.CreateZero(numDimensions, constrainedDofOrdering.NumConstrainedDofs);
            CalculateKinematicsOfEdge(constrainedDofOrdering, nodesXmin, kinematicRelations);
            CalculateKinematicsOfEdge(constrainedDofOrdering, nodesYmin, kinematicRelations);
            CalculateKinematicsOfEdge(constrainedDofOrdering, nodesZmin, kinematicRelations);
            CalculateKinematicsOfEdge(constrainedDofOrdering, nodesXmax, kinematicRelations);
            CalculateKinematicsOfEdge(constrainedDofOrdering, nodesYmax, kinematicRelations);
            CalculateKinematicsOfEdge(constrainedDofOrdering, nodesZmax, kinematicRelations);
            return kinematicRelations;
        }

        private void CalculateKinematicsOfEdge(ISubdomainConstrainedDofOrdering constrainedDofOrdering,
            IEnumerable<INode> edgeNodes, Matrix kinematicRelations)
        {
            foreach (INode node in edgeNodes)
            {
                int dofIdx = constrainedDofOrdering.ConstrainedDofs[node, DOFType.Temperature];
                kinematicRelations[0, dofIdx] = node.X;
                kinematicRelations[1, dofIdx] = node.Y;
                kinematicRelations[2, dofIdx] = node.Z;
            }
        }
    }
}
