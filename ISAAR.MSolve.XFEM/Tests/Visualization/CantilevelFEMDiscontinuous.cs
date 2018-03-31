﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Analysis;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Output;
using ISAAR.MSolve.XFEM.Output.VTK;
using ISAAR.MSolve.XFEM.Tensors;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Visualization
{
    class CantilevelFEMDiscontinuous
    {
        private static readonly double DIM_X = 100, DIM_Y = 20;
        private static readonly int ELEMENTS_PER_X = 50, ELEMENTS_PER_Y = 20;
        private static readonly double E = 2e6, v = 0.3;
        private static readonly double LOAD = 1e5;
        private static readonly string OUTPUT_FILE = "CantileverFEMDiscontinuous";

        public static void Main()
        {
            Model2D model = CreateModel();
            VectorMKL solution = Solve(model);
            WriteOutput(model, solution);
        }

        private static Model2D CreateModel()
        {
            var model = new Model2D();

            // Mesh
            var meshGenerator = new UniformRectilinearMeshGenerator(DIM_X, DIM_Y, ELEMENTS_PER_X, ELEMENTS_PER_Y);
            Tuple<XNode2D[], List<XNode2D[]>> meshEntities = meshGenerator.CreateMesh();

            // Nodes
            foreach (XNode2D node in meshEntities.Item1) model.AddNode(node);

            // Elements
            var integration = new XSimpleIntegration2D();
            //var integration = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
            //        new RectangularSubgridIntegration2D<XContinuumElement2D>(1, GaussLegendre2D.Order2x2));

            foreach (XNode2D[] elementNodes in meshEntities.Item2)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStrain(E, v);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                    elementNodes, materialField, integration));
            }

            // Constraints
            var finder = new EntityFinder(model, 1e-6);
            foreach (var node in finder.FindNodesWithX(0.0))
            {
                model.AddConstraint(node, StandardDOFType.X, 0.0);
                model.AddConstraint(node, StandardDOFType.Y, 0.0);
            }

            // Loads
            //var loadedNodes = finder.FindNodesWithX(DIM_X);
            //double loadPerNode = - LOAD / loadedNodes.Count;
            //foreach (var node in loadedNodes) model.AddNodalLoad(node, StandardDOFType.Y, loadPerNode);
            XNode2D topRightNode = finder.FindNodeWith(DIM_X, DIM_Y);
            model.AddNodalLoad(topRightNode, StandardDOFType.Y, -LOAD);

            return model;
        }

        private static VectorMKL Solve(Model2D model)
        {
            model.EnumerateDofs();
            var analysis = new LinearStaticAnalysisSkyline(model);
            analysis.Solve();
            return analysis.Solution;
        }

        private static void WriteOutput(Model2D model, VectorMKL solution)
        {
            var writer = new DiscontinuousMeshVTKWriter(model);
            writer.InitializeFile(OUTPUT_FILE);

            double[] solutionArry = solution.CopyToArray(); //TODO: change the framework to work with Vector instead of double[]

            var displacementsOut = new DisplacementOutput(model);

            IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<double[]>> elementDisplacements =
               displacementsOut.FindElementWiseDisplacements(solutionArry);
            writer.WriteVector2DField("displacements", elementDisplacements);
            
            var stressRecovery = new StressRecovery(model);
            IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>> elementStresses = 
                stressRecovery.ComputeElementWiseStresses(solutionArry);
            writer.WriteTensor2DField("stresses", elementStresses);

            writer.CloseCurrentFile();
        }
    }
}
