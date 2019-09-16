using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
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
using System.Diagnostics;


namespace ISAAR.MSolve.Tests.FEM
{
    public class NCandlesProblem
    {
        private const int subdomainID = 0;
        private const double minX = 0, minY = 0, maxX = 1, maxY = 1;
        private const double thickness = 1.0;
        private const int numElementsX = 100, numElementsY = 100;
        private const double conductivityMatrix = 1;

        [Fact]
        public static void RunExample()
        {
            Model_v2 model = new Model_v2();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Choose model
            PlateThermalExample(model);

            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var provider = new ProblemThermal_v2(model, solver);

            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

        }

        private static void AddElements(Model_v2 model)
        {
            // Material
            double density = 1.0;
            double c = 1.0;

            // Generate mesh
            var meshGenerator = new UniformMeshGenerator2D_v2(minX, minY, maxX, maxY, numElementsX, numElementsY);
            (IReadOnlyList<Node_v2> vertices, IReadOnlyList<CellConnectivity_v2> cells) = meshGenerator.CreateMesh();

            // Add nodes to the model
            for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

            // Add the elements to the model
            var elementFactory = new ThermalElement2DFactory(1.0, new ThermalMaterial(density, c, conductivityMatrix));
            for (int e = 0; e < cells.Count; ++e)
            {
                ThermalElement2D element = elementFactory.CreateElement(cells[e].CellType, cells[e].Vertices);
                var elementWrapper = new Element_v2() { ID = e, ElementType = element };
                foreach (Node_v2 node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(elementWrapper.ID, elementWrapper);
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }
        }


        private static void PlateThermalExample(Model_v2 model)
        {
            AddElements(model);
            // Dirichlet BC
            model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[6].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });

            // Neumann BC
            double q = 50.0;
            model.Loads.Add(new Load_v2() { Amount = q / 2.0, Node = model.NodesDictionary[2], DOF = DOFType.Temperature });
            model.Loads.Add(new Load_v2() { Amount = q, Node = model.NodesDictionary[5], DOF = DOFType.Temperature });
            model.Loads.Add(new Load_v2() { Amount = q / 2.0, Node = model.NodesDictionary[8], DOF = DOFType.Temperature });
        }
    }
}
