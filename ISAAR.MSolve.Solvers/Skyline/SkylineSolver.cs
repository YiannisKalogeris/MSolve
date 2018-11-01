﻿using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: factorization should be done during Solve() to time correctly. Alternatively, an observer should record the durations.
//TODO: investigate if it is possible to avoid casting the matrix provided by the analyzer/assembler into skyline. Perhaps the
//      the matrix could be obtained directly from the assembler and the analyzer could provide delegates with the operations it 
//      wants done on the matrix, instead of doing them itself.
//TODO: try to abstract the subdomain logic from ther analyzers. I doubt it is possible though.
//TODO: directly pass the single linear system instead of a list that must be checked. The same holds for all solvers and 
//      assemblers.
namespace ISAAR.MSolve.Solvers.Skyline
{
    /// <summary>
    /// Direct solver for models with only 1 subdomain. Uses Cholesky factorization on symmetric positive definite matrices 
    /// stored in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineSolver: ISolver_v2
    {
        private const string name = "SkylineSolver"; // for error messages
        private readonly SkylineAssembler assembler = new SkylineAssembler();
        private readonly FreeDofOrderer dofOrderer; //TODO: this should probably be accessed from the subdomain
        private readonly double factorizationPivotTolerance;
        private readonly LinearSystem_v2<SkylineMatrix, Vector> linearSystem;
        private CholeskySkyline factorizedMatrix;
        private LinearSystem_v2<SkylineMatrix, Vector> linearSystem_v2;

        public SkylineSolver(int subdomainID, double factorizationPivotTolerance = 1E-15) //TODO: subdomainID should not be provided by the user or needed at all
        {
            this.linearSystem = new LinearSystem_v2<SkylineMatrix, Vector>(subdomainID);
            this.LinearSystems = new Dictionary<int, ILinearSystem_v2>(1) { { subdomainID, linearSystem } };
            this.factorizationPivotTolerance = factorizationPivotTolerance;
        }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider elementMatrixProvider)
        {
            return assembler.BuildGlobalMatrix(subdomain.ΙElementsDictionary.Values, dofOrderer, elementMatrixProvider);
        }

        public void Initialize()
        {
            //TODO: perhaps I should order the dofs here.
        }

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public void Solve()
        {
            if (linearSystem.IsMatrixModified)
            {
                factorizedMatrix = linearSystem.Matrix.FactorCholesky(true, factorizationPivotTolerance);
                linearSystem.IsMatrixModified = false; //TODO: this is bad, since someone else might see it as unchanged. Better use observers.
                linearSystem.IsMatrixFactorized = true;
            }
            linearSystem.Solution = factorizedMatrix.SolveLinearSystem(linearSystem.RhsVector);
        }
    }
}
