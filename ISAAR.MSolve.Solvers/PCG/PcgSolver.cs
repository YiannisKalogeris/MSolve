﻿using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.CG;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: perhaps the user should choose the PCG settings himself and pass it. In this case, this should be named IterativeSolver.
//TODO: the maxIterations of PCG should be able to use the order of the matrix as a default value.
//TODO: IIndexable2D is not a good choice if all solvers must cast it to the matrix types the operate on.
namespace ISAAR.MSolve.Solvers.PCG
{
    /// <summary>
    /// Iterative solver for models with only 1 subdomain. Uses the Proconditioned Conjugate Gradient algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgSolver: ISolver_v2
    {
        private const string name = "PcgSolver"; // for error messages
        private readonly IGlobalMatrixAssembler<IMatrix> assembler;
        private readonly FreeDofOrderer dofOrderer; //TODO: this should probably be accessed from the subdomain
        private readonly LinearSystem_v2<IMatrix, IVector> linearSystem;
        private readonly PreconditionedConjugateGradient pcgAlgorithm;
        private readonly IPreconditionerBuilder preconditionerBuilder;
        private IPreconditioner preconditioner;

        public PcgSolver(LinearSystem_v2<IMatrix, IVector> linearSystem,
            int maxIterations, double residualTolerance, IPreconditionerBuilder preconditionerBuilder,
            IGlobalMatrixAssembler<IMatrix> globalMatrixAssembler)
        {
            this.linearSystem = linearSystem;
            this.LinearSystems = new Dictionary<int, ILinearSystem_v2>(1) { { linearSystem.ID, linearSystem } };
            this.pcgAlgorithm = new PreconditionedConjugateGradient(maxIterations, residualTolerance);
            this.preconditionerBuilder = preconditionerBuilder;
            this.assembler = globalMatrixAssembler;
        }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider elementMatrixProvider)
        {
            return assembler.BuildGlobalMatrix(subdomain.ΙElementsDictionary.Values, dofOrderer, elementMatrixProvider);
        }

        public void Initialize()
        {
        }

        /// <summary>
        /// Solves the linear system with PCG method. If the matrix has been modified, a new preconditioner will be computed.
        /// </summary>
        public void Solve()
        {
            if (linearSystem.IsMatrixModified)
            {
                preconditioner = preconditionerBuilder.BuildPreconditioner(linearSystem.Matrix);
                linearSystem.IsMatrixModified = false;
            }
            (IVector solution, CGStatistics stats) = 
                pcgAlgorithm.Solve(linearSystem.Matrix, linearSystem.RhsVector, preconditioner);
            linearSystem.Solution = solution;
        }
    }
}
