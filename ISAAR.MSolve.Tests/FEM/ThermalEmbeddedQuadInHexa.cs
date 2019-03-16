using ISAAR.MSolve.FEM.Embedding;
using Xunit;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Analyzers.Multiscale;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Problems;
using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Geometry.Shapes;
using System.Diagnostics;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ThermalEmbeddedQuadInHexa
    {
        private const int embeddedNode1ID = 1000;
        private const int embeddedNode2ID = 1001;
        private const int embeddedNode3ID = 1002;
        private const int embeddedNode4ID = 1003;
        private const int subdomainID = 0;
        private const int hostElementsIDStart = 0;
        private const int embeddedElementsIDStart = 1;
        private const double minX = -1, minY = -1, minZ = -1, maxX = 1, maxY = 1, maxZ = 1;
        private const double embeddedThickness = 0.10;
        private const int numElementsX = 2, numElementsY = 2, numElementsZ = 2;
        private static readonly Vector3 temperatureGradient = Vector3.Create(100.0, 0, 0);
        private const double conductivityMatrix = 1.0, conductivitySheet = 1000.0;

        [Fact]
        public static void ThermalEmbeddedElementExample()
        {
            Model_v2 model = new Model_v2();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Choose model
            ThermalExampleWithEmbedded(model);

            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var provider = new ProblemThermal_v2(model, solver);
            var rve = new ThermalCubeRve(model.NodesDictionary.Where(x => x.Key < embeddedNode1ID).Select(kv => kv.Value),
                Vector3.Create(minX, minY, minZ), Vector3.Create(maxX, maxY, maxZ), temperatureGradient);
            var homogenization = new HomogenizationAnalyzer(model, solver, provider, rve);

            homogenization.Initialize();
            homogenization.Solve();

            IMatrix conductivity = homogenization.EffectiveConstitutiveTensors[subdomainID];
            Debug.WriteLine($"C = [ {conductivity[0, 0]} {conductivity[0, 1]} {conductivity[0, 2]};");
            Debug.WriteLine($"      {conductivity[1, 0]} {conductivity[1, 1]} {conductivity[1, 2]};");
            Debug.WriteLine($"      {conductivity[2, 0]} {conductivity[2, 1]} {conductivity[2, 2]}]");
        }

        private static void AddHostElements(Model_v2 model)
        {
            // Material
            double density = 1.0;
            double k = 1.0;
            double c = 1.0;

            // Generate mesh
            var meshGenerator = new UniformMeshGenerator3D(minX, minY, minZ, maxX, maxY, maxZ, numElementsX, numElementsY, numElementsZ);
            (IReadOnlyList<Node_v2> vertices, IReadOnlyList<CellConnectivity_v2> cells) = meshGenerator.CreateMesh();

            // Add nodes to the model
            for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

            // Add the elements to the model
            var elementFactory = new ThermalElement3DFactory(new ThermalMaterial(density, c, conductivityMatrix));
            for (int e = 0; e < cells.Count; ++e)
            {
                ThermalElement3D element = elementFactory.CreateElement(cells[e].CellType, cells[e].Vertices);
                var elementWrapper = new Element_v2() { ID = e + hostElementsIDStart, ElementType = element };
                foreach (Node_v2 node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(elementWrapper.ID, elementWrapper);
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }
        }

        private static void AddEmbeddedElements(Model_v2 model)
        {
            // Material
            double density = 1.0;
            double c = 1.0;
            var embeddedMaterial = new ThermalMaterial(density, c, conductivitySheet);

            // Nodes
            int numNonEmbeddedNodes = model.NodesDictionary.Count;
            //TODO: if we define the element on the XZ plane, it will not work. The element will only access node.X, node.Y.
            var embeddedNodes = new Node_v2[]
            {
                new Node_v2() { ID = embeddedNode1ID, X = minX, Y = minY, Z = minZ },
                new Node_v2() { ID = embeddedNode2ID, X = maxX, Y = minY, Z = minZ },
                new Node_v2() { ID = embeddedNode3ID, X = maxX, Y = maxY, Z = minZ  },
                new Node_v2() { ID = embeddedNode4ID, X = minX, Y = maxY, Z = minZ  }
            };
            foreach (var node in embeddedNodes) model.NodesDictionary.Add(node.ID, node);

            // Elements
            var elementFactory = new ThermalElement2DFactory(embeddedThickness, embeddedMaterial);
            var elementType = elementFactory.CreateElement(CellType.Quad4, embeddedNodes);
            int numNonEmbeddedElements = model.ElementsDictionary.Count();
            int embeddedElementID = hostElementsIDStart + numNonEmbeddedElements;
            var elementWrapper = new Element_v2() { ID = embeddedElementID, ElementType = elementType };
            foreach (var node in elementType.Nodes) elementWrapper.AddNode(node);
            model.ElementsDictionary[elementWrapper.ID] = elementWrapper;
            model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);

            // Apply embedding
            var embeddedGrouping = new ThermalEmbeddedGrouping(model,
                model.ElementsDictionary.Where(x => x.Key < numNonEmbeddedElements).Select(kv => kv.Value),
                model.ElementsDictionary.Where(x => x.Key >= numNonEmbeddedElements).Select(kv => kv.Value),
                ThermalHostElementTransformationVector.CreateHost3D());
            embeddedGrouping.ApplyEmbedding();
        }


        private static void ThermalExampleWithEmbedded(Model_v2 model)
        {
            AddHostElements(model);
            AddEmbeddedElements(model);
        }
    }
}
