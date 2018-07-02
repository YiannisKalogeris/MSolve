﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;

namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
    /// <summary>
    /// Collection of integration points that are generated by a traditional 1D quadrature rule, independent of the 
    /// element type. Thus these points can be cached as static fields of each enum class, so that accessing them is fast 
    /// and storing them is done only once for all elements. All integration points are with respect to a natural 
    /// coordinate system.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IQuadrature1D
    {
        /// <summary>
        /// The integration points are sorted based on an order strictly defined for each quadrature.
        /// </summary>
        IReadOnlyList<GaussPoint1D> IntegrationPoints { get; }
    }
}
