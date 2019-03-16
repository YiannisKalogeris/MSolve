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
        private const double thickness = 1.0;
        private const int numElementsX = 2, numElementsY = 2;
        private static readonly Vector2 temperatureGradient = Vector2.Create(100.0, 0);
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
            var rve = new ThermalSquareRve(model.NodesDictionary.Where(x => x.Key < embeddedNode1ID).Select(kv => kv.Value),
                Vector2.Create(minX, minY), Vector2.Create(maxX, maxY), thickness, temperatureGradient);
            var homogenization = new HomogenizationAnalyzer(model, solver, provider, rve);

            homogenization.Initialize();
            homogenization.Solve();

            IMatrix conductivity = homogenization.EffectiveConstitutiveTensors[subdomainID];
            Console.WriteLine();
        }

        private static void AddHostElements(Model_v2 model)
        {
            // Material
            double density = 1.0;
            double k = 1.0;
            double c = 1.0;

            // Nodes
            int numNodes = 27;
            var nodes = new Node_v2[numNodes];
            nodes[0] = new Node_v2 { ID = 0,   X = maxX,              Y = maxY,              Z = maxZ              };
            nodes[1] = new Node_v2 { ID = 1,   X = maxX,              Y = (maxY + minY) / 2, Z = maxZ              };
            nodes[2] = new Node_v2 { ID = 2,   X = maxX,              Y = minY,              Z = maxZ              };
            nodes[3] = new Node_v2 { ID = 3,   X = maxX,              Y = maxY,              Z = (maxZ + minZ) / 2 };
            nodes[4] = new Node_v2 { ID = 4,   X = maxX,              Y = (maxY + minY) / 2, Z = (maxZ + minZ) / 2 };
            nodes[5] = new Node_v2 { ID = 5,   X = maxX,              Y = minY,              Z = (maxZ + minZ) / 2 };
            nodes[6] = new Node_v2 { ID = 6,   X = maxX,              Y = maxY,              Z = minZ              };
            nodes[7] = new Node_v2 { ID = 7,   X = maxX,              Y = (maxY + minY) / 2, Z = minZ              };
            nodes[8] = new Node_v2 { ID = 8,   X = maxX,              Y = minY,              Z = minZ              };
            nodes[9] = new Node_v2 { ID = 9,   X = (maxX + minX) / 2, Y = maxY,              Z = maxZ              };
            nodes[10] = new Node_v2 { ID = 10, X = (maxX + minX) / 2, Y = (maxY + minY) / 2, Z = maxZ              };
            nodes[11] = new Node_v2 { ID = 11, X = (maxX + minX) / 2, Y = minY,              Z = maxZ              };
            nodes[12] = new Node_v2 { ID = 12, X = (maxX + minX) / 2, Y = maxY,              Z = (maxZ + minZ) / 2 };
            nodes[13] = new Node_v2 { ID = 13, X = (maxX + minX) / 2, Y = (maxY + minY) / 2, Z = (maxZ + minZ) / 2 };
            nodes[14] = new Node_v2 { ID = 14, X = (maxX + minX) / 2, Y = minY,              Z = (maxZ + minZ) / 2 };
            nodes[15] = new Node_v2 { ID = 15, X = (maxX + minX) / 2, Y = maxY,              Z = minZ              };
            nodes[16] = new Node_v2 { ID = 16, X = (maxX + minX) / 2, Y = (maxY + minY) / 2, Z = minZ              };
            nodes[17] = new Node_v2 { ID = 17, X = (maxX + minX) / 2, Y = minY,              Z = minZ              };
            nodes[18] = new Node_v2 { ID = 18, X = minX,              Y = maxY,              Z = maxZ              };
            nodes[19] = new Node_v2 { ID = 19, X = minX,              Y = (maxY + minY) / 2, Z = maxZ              };
            nodes[20] = new Node_v2 { ID = 20, X = minX,              Y = minY,              Z = maxZ              };
            nodes[21] = new Node_v2 { ID = 21, X = minX,              Y = maxY,              Z = (maxZ + minZ) / 2 };
            nodes[22] = new Node_v2 { ID = 22, X = minX,              Y = (maxY + minY) / 2, Z = (maxZ + minZ) / 2 };
            nodes[23] = new Node_v2 { ID = 23, X = minX,              Y = minY,              Z = (maxZ + minZ) / 2 };
            nodes[24] = new Node_v2 { ID = 24, X = minX,              Y = maxY,              Z = minZ              };
            nodes[25] = new Node_v2 { ID = 25, X = minX,              Y = (maxY + minY) / 2, Z = minZ              };
            nodes[26] = new Node_v2 { ID = 26, X = minX,              Y = minY,              Z = minZ              };
            for (int i = 0; i < numNodes; ++i) model.NodesDictionary[i] = nodes[i];

            // Elements
            int numElements = 8;
            var elementFactory = new ThermalElement3DFactory(new ThermalMaterial(density, c, conductivityMatrix));
            var elements = new ThermalElement3D[8];
            elements[0] = elementFactory.CreateElement(CellType.Hexa8, new Node_v2[] { nodes[13], nodes[4], nodes[3], nodes[12], nodes[10], nodes[1], nodes[0], nodes[9] });
            elements[1] = elementFactory.CreateElement(CellType.Hexa8, new Node_v2[] { nodes[14], nodes[5], nodes[4], nodes[13], nodes[11], nodes[2], nodes[1], nodes[10] });
            elements[2] = elementFactory.CreateElement(CellType.Hexa8, new Node_v2[] { nodes[16], nodes[7], nodes[6], nodes[15], nodes[13], nodes[4], nodes[3], nodes[12] });
            elements[3] = elementFactory.CreateElement(CellType.Hexa8, new Node_v2[] { nodes[17], nodes[8], nodes[7], nodes[16], nodes[14], nodes[5], nodes[4], nodes[13] });
            elements[4] = elementFactory.CreateElement(CellType.Hexa8, new Node_v2[] { nodes[22], nodes[13], nodes[12], nodes[21], nodes[19], nodes[10], nodes[9], nodes[18] });
            elements[5] = elementFactory.CreateElement(CellType.Hexa8, new Node_v2[] { nodes[23], nodes[14], nodes[13], nodes[22], nodes[20], nodes[11], nodes[10], nodes[19] });
            elements[6] = elementFactory.CreateElement(CellType.Hexa8, new Node_v2[] { nodes[25], nodes[16], nodes[15], nodes[24], nodes[22], nodes[13], nodes[12], nodes[21] });
            elements[7] = elementFactory.CreateElement(CellType.Hexa8, new Node_v2[] { nodes[26], nodes[17], nodes[16], nodes[25], nodes[23], nodes[14], nodes[13], nodes[22] });

            for (int i = 0; i < numElements; ++i)
            {
                var elementWrapper = new Element_v2() { ID = i, ElementType = elements[i] };
                foreach (var node in elements[i].Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary[i] = elementWrapper;
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

            model.NodesDictionary.Add(embeddedNode1ID, new Node_v2() { ID = embeddedNode1ID, X = minX, Y = minY });
            model.NodesDictionary.Add(embeddedNode2ID, new Node_v2() { ID = embeddedNode2ID, X = minX, Y = maxY });

            // Elements
            Node_v2[] startEndNodes = { model.NodesDictionary[embeddedNode1ID], model.NodesDictionary[embeddedNode2ID] };
            var elementType = new ThermalRod(startEndNodes, crossSectionArea, embeddedMaterial);
            int numNonEmbeddedElements = model.ElementsDictionary.Count();
            int embeddedElementID = hostElementsIDStart + numNonEmbeddedElements;
            var elementWrapper = new Element_v2() { ID = embeddedElementID, ElementType = elementType };
            foreach (var node in startEndNodes) elementWrapper.AddNode(node);
            model.ElementsDictionary[elementWrapper.ID] = elementWrapper;
            model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);

            // Apply embedding
            var embeddedGrouping = new ThermalEmbeddedGrouping(model,
                model.ElementsDictionary.Where(x => x.Key < numNonEmbeddedElements).Select(kv => kv.Value),
                model.ElementsDictionary.Where(x => x.Key >= numNonEmbeddedElements).Select(kv => kv.Value),
                new ThermalElementTransformationVector());
            embeddedGrouping.ApplyEmbedding();
        }


        private static void ThermalExampleWithEmbedded(Model_v2 model)
        {
            AddHostElements(model);
            AddEmbeddedElements(model);
        }
    }
}
