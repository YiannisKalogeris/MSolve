using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Embedding;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEmbeddedElement_v2
    {
        IList<EmbeddedNode_v2> EmbeddedNodes { get; }

        /// <summary>
        /// For the provided <paramref name="node"/>, find and return the dofs of this element and their indices in the local matrices and vectors.
        /// </summary>
        /// <param name="element"></param>
        /// <param name="node"></param>
        /// <exception cref="ArgumentException">Thrown if <paramref name="element"/> does not contain <paramref name="node"/>.</exception>
        Dictionary<DOFType, int> GetInternalNodalDOFs(Element_v2 element, Node_v2 node);


        double[] GetLocalDOFValues(Element_v2 hostElement, double[] hostDOFValues);
    }
}
