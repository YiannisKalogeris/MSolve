using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Elements;

namespace ISAAR.MSolve.FEM.Embedding
{
    public abstract class ThermalHostElementTransformationVector : IEmbeddedDOFInHostTransformationVector_v2
    {
        private const int commonDofsPerNode = 1;
        private const int hostDofsPerNode = 1;

        private ThermalHostElementTransformationVector() { }

        public static ThermalHostElementTransformationVector CreateHost2D() => new ThermalHost2DTransformationVector();

        public static ThermalHostElementTransformationVector CreateHost3D() => new ThermalHost3DTransformationVector();

        private readonly DOFType[] thermalDOFTypes = new DOFType[] { DOFType.Temperature };

        public IList<DOFType> GetDependentDOFTypes { get { return thermalDOFTypes; } }

        public IList<IList<DOFType>> GetDOFTypesOfHost(EmbeddedNode_v2 node)
        {
            return node.EmbeddedInElement.ElementType.GetElementDOFTypes(node.EmbeddedInElement);
        }

        public double[][] GetTransformationVector(EmbeddedNode_v2 node)
        {
            //CheckElementType(node.EmbeddedInElement.ElementType);

            int hostShapeFunctionLength = this.NumShapeFunctionsOfElement(node.EmbeddedInElement.ElementType);
            //const int hostShapeFunctionLength = 4; //TODO: Use the interpolation for this. Probably for the next line too.
            double[] hostShapeFunctions = ((IEmbeddedHostElement_v2)node.EmbeddedInElement.ElementType).GetShapeFunctionsForNode(node.EmbeddedInElement, node);

            var transformation = new double[commonDofsPerNode][];
            for (int j = 0; j < commonDofsPerNode; j++)
            {
                transformation[j] = new double[hostShapeFunctionLength * hostDofsPerNode];
                for (int k = 0; k < hostShapeFunctionLength; k++) transformation[j][hostDofsPerNode * k + j] = hostShapeFunctions[k];
            }

            return transformation;
        }

        protected abstract int NumShapeFunctionsOfElement(IFiniteElement_v2 element);

        private class ThermalHost2DTransformationVector : ThermalHostElementTransformationVector
        {
            protected override int NumShapeFunctionsOfElement(IFiniteElement_v2 element)
                => ((ThermalElement2D)element).Interpolation.NumFunctions;
        }

        private class ThermalHost3DTransformationVector : ThermalHostElementTransformationVector
        {
            protected override int NumShapeFunctionsOfElement(IFiniteElement_v2 element)
                => ((ThermalElement3D)element).Interpolation.NumFunctions;
        }
       
    }
}

