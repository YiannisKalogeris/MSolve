﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.FETI
{
    interface IFetiPreconditioner
    {
        void SolveLinearSystem(Vector rhs, Vector lhs);
    }
}