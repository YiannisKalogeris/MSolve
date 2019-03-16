using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Embedding;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEmbeddedHostElement_v2
    {
        /// <summary>
        /// Maps the global cartesian coordinates of <paramref name="node"/> to the local natural (e.g. isoparametric) system of the element.
        /// </summary>
        /// <param name="element"></param>
        /// <param name="node"></param>
        /// <param name="transformationVector"></param>
        EmbeddedNode_v2 BuildHostElementEmbeddedNode(Element_v2 element, Node_v2 node,
            IEmbeddedDOFInHostTransformationVector_v2 transformationVector);

        /// <summary>
        /// Returns the shape functions of the element evaluated at the natural coordinates defined by <paramref name="node"/>. 
        /// Optionally in the same array there could be the shape functions derivatives and the inverse jacobian, evaluated at the same point.
        /// </summary>
        /// <param name="element"></param>
        /// <param name="node">
        /// Its coordinates are defined in the natural (e.g. isoparametric) system of the element.
        /// </param>
        double[] GetShapeFunctionsForNode(Element_v2 element, EmbeddedNode_v2 node); 
    }
}
