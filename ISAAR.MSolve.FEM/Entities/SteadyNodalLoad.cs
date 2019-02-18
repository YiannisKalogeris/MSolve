﻿using ISAAR.MSolve.Discretization.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: This is probably covered by Load.cs
namespace ISAAR.MSolve.FEM.Entities
{
    /// <summary>
    /// You should use <see cref="Load"/> instead.
    /// </summary>
    public class SteadyNodalLoad: ITimeDependentNodalLoad
    {
        private readonly double constantloadAmount;

        public SteadyNodalLoad(double constantloadAmount)
        {
            this.constantloadAmount = constantloadAmount;
        }

        public Node Node { get; set; }
        public DOFType DOF { get; set; }

        public double GetLoadAmount(int timeStep)
        {
            return constantloadAmount;
        }
    }
}