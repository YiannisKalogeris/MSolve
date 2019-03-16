using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: This could be done simultaneously with ordering the free dofs, to improve performance.
namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Orders the constrained dofs of a subdomain, independendtly from the free ones.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class ConstrainedDofOrderingStrategy
    {
        /// <summary>
        /// Orders the constrained freedom degrees of one of the model's subdomains.
        /// </summary>
        /// <param name="subdomain">A subdomain of the whole model.</param>
        internal (int numSubdomainConstrainedDofs, DofTable subdomainConstrainedDofs) OrderSubdomainDofs(ISubdomain_v2 subdomain)
        {
            //TODO: This should be done simultaneously with ordering the free dofs.

            int totalDOFs = 0;
            Dictionary<int, List<DOFType>> nodalDOFTypesDictionary = new Dictionary<int, List<DOFType>>(); //TODO: use Set instead of List
            foreach (IElement_v2 element in subdomain.Elements)
            {
                for (int i = 0; i < element.Nodes.Count; i++)
                {
                    if (!nodalDOFTypesDictionary.ContainsKey(element.Nodes[i].ID))
                        nodalDOFTypesDictionary.Add(element.Nodes[i].ID, new List<DOFType>());
                    nodalDOFTypesDictionary[element.Nodes[i].ID].AddRange(element.ElementType.DofEnumerator.GetDOFTypesForDOFEnumeration(element)[i]);
                }
            }

            var constrainedDofs = new DofTable();
            foreach (INode node in subdomain.Nodes)
            {
                //List<DOFType> dofTypes = new List<DOFType>();
                //foreach (Element element in node.ElementsDictionary.Values)
                //{
                //    if (elementsDictionary.ContainsKey(element.ID))
                //    {
                //        foreach (DOFType dof in element.ElementType.DOFTypes)
                //            dofTypes.Add(dof);
                //    }
                //}

                Dictionary<DOFType, int> dofsDictionary = new Dictionary<DOFType, int>();
                foreach (DOFType dofType in nodalDOFTypesDictionary[node.ID].Distinct<DOFType>())
                {
                    foreach (var constraint in node.Constraints) //TODO: access the constraints from the subdomain
                    {
                        if (constraint.DOF == dofType)
                        {
                            constrainedDofs[node, dofType] = totalDOFs++;
                            break;
                        }
                    }
                }
            }

            return (totalDOFs, constrainedDofs);
        }

        ///// <summary>
        ///// Orders the constrained freedom degrees of one of the model's subdomains.
        ///// </summary>
        ///// <param name="subdomain">A subdomain of the whole model.</param>
        //internal (int numSubdomainConstrainedDofs, DofTable subdomainConstrainedDofs) OrderSubdomainDofs(ISubdomain_v2 subdomain)
        //{
        //    //TODO: This doesn't work if an embedded node has been constrained. (Yes it happens)
        //    var constrainedDofs = new DofTable();
        //    int dofCounter = 0;
        //    foreach (INode node in subdomain.Nodes)
        //    {
        //        foreach (Constraint constraint in node.Constraints) constrainedDofs[node, constraint.DOF] = dofCounter++;
        //    }
        //    return (dofCounter, constrainedDofs);
        //}
    }
}
