package org.openscience.cdk.isomorphism;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.chesmapper.map.main.Settings;
import org.mg.javalib.util.ArrayUtil;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.isomorphism.mcss.RGraph;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.isomorphism.mcss.RNode;
import org.openscience.cdk.tools.manipulator.BondManipulator;

public class MyUniversalIsomorphismTester
{

	final static int ID1 = 0;
	final static int ID2 = 1;
	private static long start;
	private static long timeout = -1;

	///////////////////////////////////////////////////////////////////////////
	//                            Query Methods
	//
	// This methods are simple applications of the RGraph model on atom containers
	// using different constrains and search options. They give an exemple of the
	// most common queries but of course it is possible to define other type of
	// queries exploiting the constrain and option combinations
	//

	////
	// Isomorphism search

	/**
	 * Tests if g1 and g2 are isomorph.
	 *
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  second molecule. May be an {@link IQueryAtomContainer}.
	 * @return     true if the 2 molecule are isomorph
	 * @throws     CDKException if the first molecule is an instance of IQueryAtomContainer
	 */
	public static boolean isIsomorph(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		if (g1 instanceof IQueryAtomContainer)
			throw new CDKException("The first IAtomContainer must not be an IQueryAtomContainer");

		if (g2.getAtomCount() != g1.getAtomCount())
			return false;
		// check single atom case
		if (g2.getAtomCount() == 1)
		{
			IAtom atom = g1.getAtom(0);
			IAtom atom2 = g2.getAtom(0);
			if (atom instanceof IQueryAtom)
			{
				IQueryAtom qAtom = (IQueryAtom) atom;
				return qAtom.matches(g2.getAtom(0));
			}
			else if (atom2 instanceof IQueryAtom)
			{
				IQueryAtom qAtom = (IQueryAtom) atom2;
				return qAtom.matches(g1.getAtom(0));
			}
			else
			{
				String atomSymbol = atom.getSymbol();
				return g1.getAtom(0).getSymbol().equals(atomSymbol);
			}
		}
		return (getIsomorphMap(g1, g2) != null);
	}

	/**
	 * Returns the first isomorph mapping found or null.
	 *
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  second molecule. May be an {@link IQueryAtomContainer}.
	 * @return     the first isomorph mapping found projected of g1. This is a List of RMap objects containing Ids of matching bonds.
	 */
	public static List<RMap> getIsomorphMap(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		if (g1 instanceof IQueryAtomContainer)
			throw new CDKException("The first IAtomContainer must not be an IQueryAtomContainer");

		List<RMap> result = null;

		List<List<RMap>> rMapsList = search(g1, g2, getBitSet(g1), getBitSet(g2), false, false);

		if (!rMapsList.isEmpty())
		{
			result = rMapsList.get(0);
		}

		return result;
	}

	/**
	 * Returns the first isomorph 'atom mapping' found for g2 in g1.
	 *
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  second molecule. May be an {@link IQueryAtomContainer}.
	 * @return     the first isomorph atom mapping found projected on g1.
	 * This is a List of RMap objects containing Ids of matching atoms.
	 * @throws CDKException if the first molecules is not an instance of {@link IQueryAtomContainer}
	 */
	public static List<RMap> getIsomorphAtomsMap(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		if (g1 instanceof IQueryAtomContainer)
			throw new CDKException("The first IAtomContainer must not be an IQueryAtomContainer");

		List<RMap> list = checkSingleAtomCases(g1, g2);
		if (list == null)
		{
			return makeAtomsMapOfBondsMap(MyUniversalIsomorphismTester.getIsomorphMap(g1, g2), g1, g2);
		}
		else if (list.isEmpty())
		{
			return null;
		}
		else
		{
			return list;
		}
	}

	/**
	 * Returns all the isomorph 'mappings' found between two
	 * atom containers.
	 *
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  second molecule. May be an {@link IQueryAtomContainer}.
	 * @return     the list of all the 'mappings'
	 */
	public static List<List<RMap>> getIsomorphMaps(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		return search(g1, g2, getBitSet(g1), getBitSet(g2), true, true);
	}

	public static BitSet getBitSet(IAtomContainer ac, IBond bonds[])
	{
		BitSet bs;
		int n = ac.getBondCount();
		if (n != 0)
		{
			bs = new BitSet(n);
			for (int i = 0; i < n; i++)
				if (ArrayUtil.indexOf(bonds, ac.getBond(i)) != -1)
					bs.set(i);
		}
		else
			bs = new BitSet();
		return bs;
	}

	public static List<List<RMap>> getIsomorphMaps(IAtomContainer g1, IBond b1[], IAtomContainer g2, IBond b2[])
			throws CDKException
	{
		if (g1 instanceof IQueryAtomContainer)
			throw new CDKException("The first IAtomContainer must not be an IQueryAtomContainer");
		BitSet c1 = getBitSet(g1, b1);
		BitSet c2 = getBitSet(g2, b2);
		return MyUniversalIsomorphismTester.search(g1, g2, c1, c2, true, true, c1, c2);
	}

	/////
	// Subgraph search

	/**
	 * Returns all the subgraph 'bond mappings' found for g2 in g1.
	 * This is an {@link List} of {@link List}s of {@link RMap} objects.
	 *
	 * Note that if the query molecule is a single atom, then bond mappings
	 * cannot be defined. In such a case, the {@link RMap} object refers directly to
	 * atom - atom mappings. Thus RMap.id1 is the index of the target atom
	 * and RMap.id2 is the index of the matching query atom (in this case,
	 * it will always be 0). Note that in such a case, there is no need
	 * to call {@link #makeAtomsMapsOfBondsMaps(List, IAtomContainer, IAtomContainer)},
	 * though if it is called, then the
	 * return value is simply the same as the return value of this method.
	 *
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  second molecule. May be an {@link IQueryAtomContainer}.
	 * @return     the list of all the 'mappings' found projected of g1
	 *
	 * @see #makeAtomsMapsOfBondsMaps(List, IAtomContainer, IAtomContainer)
	 */
	public static List<List<RMap>> getSubgraphMaps(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		return search(g1, g2, new BitSet(), getBitSet(g2), true, true);
	}

	/**
	 * Returns the first subgraph 'bond mapping' found for g2 in g1.
	 *
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  second molecule. May be an {@link IQueryAtomContainer}.
	 * @return     the first subgraph bond mapping found projected on g1. This is a {@link List} of
	 *             {@link RMap} objects containing Ids of matching bonds.
	 */
	public static List<RMap> getSubgraphMap(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		List<RMap> result = null;
		List<List<RMap>> rMapsList = search(g1, g2, new BitSet(), getBitSet(g2), false, false);

		if (!rMapsList.isEmpty())
		{
			result = rMapsList.get(0);
		}

		return result;
	}

	/**
	 * Returns all subgraph 'atom mappings' found for g2 in g1, where g2 must be a substructure
	 * of g1. If it is not a substructure, null will be returned.
	 * This is an {@link List} of {@link List}s of {@link RMap} objects.
	 *
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  substructure to be mapped. May be an {@link IQueryAtomContainer}.
	 * @return     all subgraph atom mappings found projected on g1. This is a
	 *             {@link List} of {@link RMap} objects containing Ids of matching atoms.
	 */
	public static List<List<RMap>> getSubgraphAtomsMaps(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		List<RMap> list = checkSingleAtomCases(g1, g2);
		if (list == null)
		{
			return makeAtomsMapsOfBondsMaps(MyUniversalIsomorphismTester.getSubgraphMaps(g1, g2), g1, g2);
		}
		else
		{
			List<List<RMap>> atomsMap = new ArrayList<List<RMap>>();
			atomsMap.add(list);
			return atomsMap;
		}
	}

	/**
	 * Returns the first subgraph 'atom mapping' found for g2 in g1, where g2 must be a substructure
	 * of g1. If it is not a substructure, null will be returned.
	 *
	 * @param  g1 first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2 substructure to be mapped. May be an {@link IQueryAtomContainer}.
	 * @return    the first subgraph atom mapping found projected on g1.
	 *            This is a {@link List} of {@link RMap} objects containing Ids of matching atoms.
	 */
	public static List<RMap> getSubgraphAtomsMap(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		List<RMap> list = checkSingleAtomCases(g1, g2);
		if (list == null)
		{
			return makeAtomsMapOfBondsMap(MyUniversalIsomorphismTester.getSubgraphMap(g1, g2), g1, g2);
		}
		else if (list.isEmpty())
		{
			return null;
		}
		else
		{
			return list;
		}
	}

	/**
	 * Tests if g2 a subgraph of g1.
	 *
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  second molecule. May be an {@link IQueryAtomContainer}.
	 * @return     true if g2 a subgraph on g1
	 */
	public static boolean isSubgraph(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		if (g1 instanceof IQueryAtomContainer)
			throw new CDKException("The first IAtomContainer must not be an IQueryAtomContainer");

		if (g2.getAtomCount() > g1.getAtomCount())
			return false;
		// test for single atom case
		if (g2.getAtomCount() == 1)
		{
			IAtom atom = g2.getAtom(0);
			for (int i = 0; i < g1.getAtomCount(); i++)
			{
				IAtom atom2 = g1.getAtom(i);
				if (atom instanceof IQueryAtom)
				{
					IQueryAtom qAtom = (IQueryAtom) atom;
					if (qAtom.matches(atom2))
						return true;
				}
				else if (atom2 instanceof IQueryAtom)
				{
					IQueryAtom qAtom = (IQueryAtom) atom2;
					if (qAtom.matches(atom))
						return true;
				}
				else
				{
					if (atom2.getSymbol().equals(atom.getSymbol()))
						return true;
				}
			}
			return false;
		}
		if (!testSubgraphHeuristics(g1, g2))
			return false;
		return (getSubgraphMap(g1, g2) != null);
	}

	////
	// Maximum common substructure search

	/**
	 * Returns all the maximal common substructure between twp atom containers.
	 *
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  second molecule. May be an {@link IQueryAtomContainer}.
	 * @return     the list of all the maximal common substructure
	 *             found projected of g1 (list of AtomContainer )
	 */
	public static List<IAtomContainer> getOverlaps(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		start = System.currentTimeMillis();
		List<List<RMap>> rMapsList = search(g1, g2, new BitSet(), new BitSet(), true, false);

		// projection on G1
		List<IAtomContainer> graphList = projectList(rMapsList, g1, ID1);

		// reduction of set of solution (isomorphism and substructure
		// with different 'mappings'

		return getMaximum(graphList);
	}

	/**
	 * Transforms an AtomContainer into a {@link BitSet} (which's size = number of bond
	 * in the atomContainer, all the bit are set to true).
	 *
	 * @param  ac  {@link IAtomContainer} to transform
	 * @return     The bitSet
	 */
	public static BitSet getBitSet(IAtomContainer ac)
	{
		BitSet bs;
		int n = ac.getBondCount();

		if (n != 0)
		{
			bs = new BitSet(n);
			for (int i = 0; i < n; i++)
			{
				bs.set(i);
			}
		}
		else
		{
			bs = new BitSet();
		}

		return bs;
	}

	//////////////////////////////////////////////////
	//          Internal methods

	/**
	 * Builds the {@link RGraph} ( resolution graph ), from two atomContainer
	 * (description of the two molecules to compare)
	 * This is the interface point between the CDK model and
	 * the generic MCSS algorithm based on the RGRaph.
	 *
	 * @param  g1  Description of the first molecule
	 * @param  g2  Description of the second molecule
	 * @return     the rGraph
	 */
	public static RGraph buildRGraph(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		return buildRGraph(g1, g2, null, null);
	}

	public static RGraph buildRGraph(IAtomContainer g1, IAtomContainer g2, BitSet c1, BitSet c2) throws CDKException
	{
		RGraph rGraph = new RGraph();
		nodeConstructor(rGraph, g1, g2, c1, c2);
		arcConstructor(rGraph, g1, g2);
		return rGraph;
	}

	/**
	 * General {@link RGraph} parsing method (usually not used directly)
	 * This method is the entry point for the recursive search
	 * adapted to the atom container input.
	 *
	 * @param  g1                first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2                second molecule. May be an {@link IQueryAtomContainer}.
	 * @param  c1                initial condition ( bonds from g1 that
	 *                           must be contains in the solution )
	 * @param  c2                initial condition ( bonds from g2 that
	 *                           must be contains in the solution )
	 * @param  findAllStructure  if false stop at the first structure found
	 * @param  findAllMap        if true search all the 'mappings' for one same
	 *                           structure
	 * @return                   a List of Lists of {@link RMap} objects that represent the search solutions
	 */
	public static List<List<RMap>> search(IAtomContainer g1, IAtomContainer g2, BitSet c1, BitSet c2,
			boolean findAllStructure, boolean findAllMap) throws CDKException
	{
		return search(g1, g2, c1, c2, findAllStructure, findAllMap, null, null);
	}

	public static List<List<RMap>> search(IAtomContainer g1, IAtomContainer g2, BitSet c1, BitSet c2,
			boolean findAllStructure, boolean findAllMap, BitSet c1exclusive, BitSet c2exclusive) throws CDKException
	{
		// remember start time
		start = System.currentTimeMillis();

		// handle single query atom case separately
		if (g2.getAtomCount() == 1)
		{
			List<List<RMap>> matches = new ArrayList<List<RMap>>();
			IAtom queryAtom = g2.getAtom(0);

			// we can have a IQueryAtomContainer *or* an IAtomContainer
			if (queryAtom instanceof IQueryAtom)
			{
				IQueryAtom qAtom = (IQueryAtom) queryAtom;
				for (IAtom atom : g1.atoms())
				{
					if (qAtom.matches(atom))
					{
						List<RMap> lmap = new ArrayList<RMap>();
						lmap.add(new RMap(g1.getAtomNumber(atom), 0));
						matches.add(lmap);
					}
				}
			}
			else
			{
				for (IAtom atom : g1.atoms())
				{
					if (queryAtom.getSymbol().equals(atom.getSymbol()))
					{
						List<RMap> lmap = new ArrayList<RMap>();
						lmap.add(new RMap(g1.getAtomNumber(atom), 0));
						matches.add(lmap);
					}
				}
			}
			return matches;
		}

		// reset result
		List<List<RMap>> rMapsList = new ArrayList<List<RMap>>();

		// build the RGraph corresponding to this problem
		RGraph rGraph = buildRGraph(g1, g2, c1exclusive, c2exclusive);
		// Set time data
		rGraph.setTimeout(MyUniversalIsomorphismTester.timeout);
		rGraph.setStart(MyUniversalIsomorphismTester.start);
		// parse the RGraph with the given constrains and options
		rGraph.parse(c1, c2, findAllStructure, findAllMap);
		List<BitSet> solutionList = rGraph.getSolutions();

		// conversions of RGraph's internal solutions to G1/G2 mappings
		for (BitSet set : solutionList)
		{
			rMapsList.add(rGraph.bitSetToRMap(set));
		}

		return rMapsList;
	}

	//////////////////////////////////////
	//    Manipulation tools

	/**
	 * Projects a list of {@link RMap} on a molecule.
	 *
	 * @param  rMapList  the list to project
	 * @param  g         the molecule on which project
	 * @param  id        the id in the {@link RMap} of the molecule g
	 * @return           an AtomContainer
	 */
	public static IAtomContainer project(List<RMap> rMapList, IAtomContainer g, int id)
	{
		IAtomContainer ac = g.getBuilder().newInstance(IAtomContainer.class);

		Map<IAtom, IAtom> table = new HashMap<IAtom, IAtom>();
		IAtom a1;
		IAtom a2;
		IAtom a;
		IBond bond;

		for (Iterator<RMap> i = rMapList.iterator(); i.hasNext();)
		{
			RMap rMap = i.next();
			if (id == MyUniversalIsomorphismTester.ID1)
			{
				bond = g.getBond(rMap.getId1());
			}
			else
			{
				bond = g.getBond(rMap.getId2());
			}

			a = bond.getAtom(0);
			a1 = (IAtom) table.get(a);

			if (a1 == null)
			{
				try
				{
					a1 = (IAtom) a.clone();
				}
				catch (CloneNotSupportedException e)
				{
					Settings.LOGGER.error(e);
				}
				ac.addAtom(a1);
				table.put(a, a1);
			}

			a = bond.getAtom(1);
			a2 = table.get(a);

			if (a2 == null)
			{
				try
				{
					a2 = (IAtom) a.clone();
				}
				catch (CloneNotSupportedException e)
				{
					Settings.LOGGER.error(e);
				}
				ac.addAtom(a2);
				table.put(a, a2);
			}
			IBond newBond = g.getBuilder().newInstance(IBond.class, a1, a2, bond.getOrder());
			newBond.setFlag(CDKConstants.ISAROMATIC, bond.getFlag(CDKConstants.ISAROMATIC));
			ac.addBond(newBond);
		}
		return ac;
	}

	/**
	 * Projects a list of RMapsList on a molecule.
	 *
	 * @param  rMapsList  list of RMapsList to project
	 * @param  g          the molecule on which project
	 * @param  id         the id in the RMap of the molecule g
	 * @return            a list of AtomContainer
	 */
	public static List<IAtomContainer> projectList(List<List<RMap>> rMapsList, IAtomContainer g, int id)
	{
		List<IAtomContainer> graphList = new ArrayList<IAtomContainer>();

		for (List<RMap> rMapList : rMapsList)
		{
			IAtomContainer ac = project(rMapList, g, id);
			graphList.add(ac);
		}
		return graphList;
	}

	/**
	 * Removes all redundant solution.
	 *
	 * @param  graphList  the list of structure to clean
	 * @return            the list cleaned
	 * @throws CDKException if there is a problem in obtaining subgraphs
	 */
	private static List<IAtomContainer> getMaximum(List<IAtomContainer> graphList) throws CDKException
	{
		List<IAtomContainer> reducedGraphList = new ArrayList<IAtomContainer>();
		reducedGraphList.addAll(graphList);

		for (int i = 0; i < graphList.size(); i++)
		{
			IAtomContainer gi = graphList.get(i);

			for (int j = i + 1; j < graphList.size(); j++)
			{
				IAtomContainer gj = graphList.get(j);

				// Gi included in Gj or Gj included in Gi then
				// reduce the irrelevant solution
				if (isSubgraph(gj, gi))
				{
					reducedGraphList.remove(gi);
				}
				else if (isSubgraph(gi, gj))
				{
					reducedGraphList.remove(gj);
				}
			}
		}
		return reducedGraphList;
	}

	/**
	 *  Checks for single atom cases before doing subgraph/isomorphism search.
	 *
	 * @param  g1  AtomContainer to match on. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  AtomContainer as query. May be an {@link IQueryAtomContainer}.
	 * @return     {@link List} of {@link List} of {@link RMap} objects for the Atoms (not Bonds!), null if no single atom case
	 * @throws     CDKException if the first molecule is an instance of IQueryAtomContainer
	*/
	public static List<RMap> checkSingleAtomCases(IAtomContainer g1, IAtomContainer g2) throws CDKException
	{
		if (g1 instanceof IQueryAtomContainer)
			throw new CDKException("The first IAtomContainer must not be an IQueryAtomContainer");

		if (g2.getAtomCount() == 1)
		{
			List<RMap> arrayList = new ArrayList<RMap>();
			IAtom atom = g2.getAtom(0);
			if (atom instanceof IQueryAtom)
			{
				IQueryAtom qAtom = (IQueryAtom) atom;
				for (int i = 0; i < g1.getAtomCount(); i++)
				{
					if (qAtom.matches(g1.getAtom(i)))
						arrayList.add(new RMap(i, 0));
				}
			}
			else
			{
				String atomSymbol = atom.getSymbol();
				for (int i = 0; i < g1.getAtomCount(); i++)
				{
					if (g1.getAtom(i).getSymbol().equals(atomSymbol))
						arrayList.add(new RMap(i, 0));
				}
			}
			return arrayList;
		}
		else if (g1.getAtomCount() == 1)
		{
			List<RMap> arrayList = new ArrayList<RMap>();
			IAtom atom = g1.getAtom(0);
			for (int i = 0; i < g2.getAtomCount(); i++)
			{
				IAtom atom2 = g2.getAtom(i);
				if (atom2 instanceof IQueryAtom)
				{
					IQueryAtom qAtom = (IQueryAtom) atom2;
					if (qAtom.matches(atom))
						arrayList.add(new RMap(0, i));
				}
				else
				{
					if (atom2.getSymbol().equals(atom.getSymbol()))
						arrayList.add(new RMap(0, i));
				}
			}
			return arrayList;
		}
		else
		{
			return null;
		}
	}

	/**
	 *  This makes maps of matching atoms out of a maps of matching bonds as produced by the
	 *  get(Subgraph|Ismorphism)Maps methods.
	 *
	 * @param  l   The list produced by the getMap method.
	 * @param  g1  The first atom container. Must not be a {@link IQueryAtomContainer}.
	 * @param  g2  The second one (first and second as in getMap). May be an {@link IQueryAtomContainer}.
	 * @return     A List of {@link List}s of {@link RMap} objects of matching Atoms.
	 */
	public static List<List<RMap>> makeAtomsMapsOfBondsMaps(List<List<RMap>> l, IAtomContainer g1, IAtomContainer g2)
	{
		if (l == null)
		{
			return l;
		}
		if (g2.getAtomCount() == 1)
			return l; // since the RMap is already an atom-atom mapping
		List<List<RMap>> result = new ArrayList<List<RMap>>();
		for (List<RMap> l2 : l)
		{
			result.add(makeAtomsMapOfBondsMap(l2, g1, g2));
		}
		return result;
	}

	/**
	 *  This makes a map of matching atoms out of a map of matching bonds as produced by the
	 *  get(Subgraph|Ismorphism)Map methods.
	 *
	 * @param  l   The list produced by the getMap method.
	 * @param  g1  first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  g2  second molecule. May be an {@link IQueryAtomContainer}.
	 * @return     The mapping found projected on g1. This is a {@link List} of {@link RMap} objects
	 *             containing Ids of matching atoms.
	 */
	public static List<RMap> makeAtomsMapOfBondsMap(List<RMap> l, IAtomContainer g1, IAtomContainer g2)
	{
		if (l == null)
			return (l);
		List<RMap> result = new ArrayList<RMap>();
		for (int i = 0; i < l.size(); i++)
		{
			IBond bond1 = g1.getBond(l.get(i).getId1());
			IBond bond2 = g2.getBond(l.get(i).getId2());
			IAtom[] atom1 = BondManipulator.getAtomArray(bond1);
			IAtom[] atom2 = BondManipulator.getAtomArray(bond2);
			for (int j = 0; j < 2; j++)
			{
				List<IBond> bondsConnectedToAtom1j = g1.getConnectedBondsList(atom1[j]);
				for (int k = 0; k < bondsConnectedToAtom1j.size(); k++)
				{
					if (bondsConnectedToAtom1j.get(k) != bond1)
					{
						IBond testBond = (IBond) bondsConnectedToAtom1j.get(k);
						for (int m = 0; m < l.size(); m++)
						{
							IBond testBond2;
							if (((RMap) l.get(m)).getId1() == g1.getBondNumber(testBond))
							{
								testBond2 = g2.getBond(((RMap) l.get(m)).getId2());
								for (int n = 0; n < 2; n++)
								{
									List<IBond> bondsToTest = g2.getConnectedBondsList(atom2[n]);
									if (bondsToTest.contains(testBond2))
									{
										RMap map;
										if (j == n)
										{
											map = new RMap(g1.getAtomNumber(atom1[0]), g2.getAtomNumber(atom2[0]));
										}
										else
										{
											map = new RMap(g1.getAtomNumber(atom1[1]), g2.getAtomNumber(atom2[0]));
										}
										if (!result.contains(map))
										{
											result.add(map);
										}
										RMap map2;
										if (j == n)
										{
											map2 = new RMap(g1.getAtomNumber(atom1[1]), g2.getAtomNumber(atom2[1]));
										}
										else
										{
											map2 = new RMap(g1.getAtomNumber(atom1[0]), g2.getAtomNumber(atom2[1]));
										}
										if (!result.contains(map2))
										{
											result.add(map2);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		return result;
	}

	/**
	 *  Builds  the nodes of the {@link RGraph} ( resolution graph ), from
	 *  two atom containers (description of the two molecules to compare)
	 *
	 * @param  gr   the target RGraph
	 * @param  ac1   first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  ac2   second molecule. May be an {@link IQueryAtomContainer}.
	 * @throws CDKException if it takes too long to identify overlaps
	 */
	private static void nodeConstructor(RGraph gr, IAtomContainer ac1, IAtomContainer ac2, BitSet c1, BitSet c2)
			throws CDKException
	{
		if (ac1 instanceof IQueryAtomContainer)
			throw new CDKException("The first IAtomContainer must not be an IQueryAtomContainer");

		// resets the target graph.
		gr.clear();

		// compares each bond of G1 to each bond of G2
		for (int i = 0; i < ac1.getBondCount(); i++)
		{
			for (int j = 0; j < ac2.getBondCount(); j++)
			{
				IBond bondA2 = ac2.getBond(j);
				if (bondA2 instanceof IQueryBond)
				{
					IQueryBond queryBond = (IQueryBond) bondA2;
					IQueryAtom atom1 = (IQueryAtom) (bondA2.getAtom(0));
					IQueryAtom atom2 = (IQueryAtom) (bondA2.getAtom(1));
					IBond bond = ac1.getBond(i);
					if (queryBond.matches(bond))
					{
						// ok, bonds match
						if (atom1.matches(bond.getAtom(0)) && atom2.matches(bond.getAtom(1))
								|| atom1.matches(bond.getAtom(1)) && atom2.matches(bond.getAtom(0)))
						{
							// ok, atoms match in either order
							gr.addNode(new RNode(i, j));
						}
					}
				}
				else
				{
					// if both bonds are compatible then create an association node
					// in the resolution graph
					if (((c1 == null || c1.get(i)) && (c2 == null || c2.get(j))) //
							&& ( // bond type conditions
							( // same bond order and same aromaticity flag (either both on or off)
							ac1.getBond(i).getOrder() == ac2.getBond(j).getOrder() && ac1.getBond(i).getFlag(
									CDKConstants.ISAROMATIC) == ac2.getBond(j).getFlag(CDKConstants.ISAROMATIC))// 
							|| ( // both bond are aromatic
							ac1.getBond(i).getFlag(CDKConstants.ISAROMATIC) && ac2.getBond(j).getFlag(
									CDKConstants.ISAROMATIC))) && ( // atome type conditions
							( // a1 = a2 && b1 = b2
							atomEquals(ac1.getBond(i).getAtom(0), ac2.getBond(j).getAtom(0)) && //
									atomEquals(ac1.getBond(i).getAtom(1), ac2.getBond(j).getAtom(1)) //
							|| ( // a1 = b2 && b1 = a2
							atomEquals(ac1.getBond(i).getAtom(0), ac2.getBond(j).getAtom(1)) && //
							atomEquals(ac1.getBond(i).getAtom(1), ac2.getBond(j).getAtom(0))))))
					{
						gr.addNode(new RNode(i, j));
					}
				}
			}
		}
	}

	private static boolean atomEquals(IAtom a1, IAtom a2)
	{
		return a1.getSymbol().equals(a2.getSymbol())
				&& a1.getFlag(CDKConstants.ISAROMATIC) == a2.getFlag(CDKConstants.ISAROMATIC);
	}

	/**
	 *  Build edges of the {@link RGraph}s.
	 *  This method create the edge of the RGraph and
	 *  calculates the incompatibility and neighborhood
	 *  relationships between RGraph nodes.
	 *
	 * @param  gr   the rGraph
	 * @param  ac1   first molecule. Must not be an {@link IQueryAtomContainer}.
	 * @param  ac2   second molecule. May be an {@link IQueryAtomContainer}.
	 * @throws CDKException if it takes too long to get the overlaps
	 */
	private static void arcConstructor(RGraph gr, IAtomContainer ac1, IAtomContainer ac2) throws CDKException
	{
		// each node is incompatible with himself
		for (int i = 0; i < gr.getGraph().size(); i++)
		{
			RNode x = (RNode) gr.getGraph().get(i);
			x.getForbidden().set(i);
		}

		IBond a1;
		IBond a2;
		IBond b1;
		IBond b2;

		gr.setFirstGraphSize(ac1.getBondCount());
		gr.setSecondGraphSize(ac2.getBondCount());

		for (int i = 0; i < gr.getGraph().size(); i++)
		{
			RNode x = gr.getGraph().get(i);

			// two nodes are neighbors if their adjacency
			// relationship in are equivalent in G1 and G2
			// else they are incompatible.
			for (int j = i + 1; j < gr.getGraph().size(); j++)
			{
				RNode y = gr.getGraph().get(j);

				a1 = ac1.getBond(gr.getGraph().get(i).getRMap().getId1());
				a2 = ac2.getBond(gr.getGraph().get(i).getRMap().getId2());
				b1 = ac1.getBond(gr.getGraph().get(j).getRMap().getId1());
				b2 = ac2.getBond(gr.getGraph().get(j).getRMap().getId2());

				if (a2 instanceof IQueryBond)
				{
					if (a1.equals(b1) || a2.equals(b2) || !queryAdjacencyAndOrder(a1, b1, a2, b2))
					{
						x.getForbidden().set(j);
						y.getForbidden().set(i);
					}
					else if (hasCommonAtom(a1, b1))
					{
						x.getExtension().set(j);
						y.getExtension().set(i);
					}
				}
				else
				{
					if (a1.equals(b1) || a2.equals(b2) || (!getCommonSymbol(a1, b1).equals(getCommonSymbol(a2, b2))))
					{
						x.getForbidden().set(j);
						y.getForbidden().set(i);
					}
					else if (hasCommonAtom(a1, b1))
					{
						x.getExtension().set(j);
						y.getExtension().set(i);
					}
				}
			}
		}
	}

	/**
	 * Determines if two bonds have at least one atom in common.
	 *
	 * @param  a  first bond
	 * @param  b  second bond
	 * @return    the symbol of the common atom or "" if
	 *            the 2 bonds have no common atom
	 */
	private static boolean hasCommonAtom(IBond a, IBond b)
	{
		return a.contains(b.getAtom(0)) || a.contains(b.getAtom(1));
	}

	/**
	 *  Determines if 2 bond have 1 atom in common and returns the common symbol.
	 *
	 * @param  a  first bond
	 * @param  b  second bond
	 * @return    the symbol of the common atom or "" if
	 *            the 2 bonds have no common atom
	 */
	private static String getCommonSymbol(IBond a, IBond b)
	{
		String symbol = "";

		if (a.contains(b.getAtom(0)))
		{
			symbol = b.getAtom(0).getSymbol();
		}
		else if (a.contains(b.getAtom(1)))
		{
			symbol = b.getAtom(1).getSymbol();
		}

		return symbol;
	}

	//	/**
	//	*  Determines if 2 bond have 1 atom in common if second is a query AtomContainer.
	//	*
	//	* @param  a1  first bond
	//	* @param  b1  second bond
	//	* @return    the symbol of the common atom or "" if
	//	*            the 2 bonds have no common atom
	//	*/
	//	private static boolean queryAdjacency(IBond a1, IBond b1, IBond a2, IBond b2)
	//	{
	//
	//		IAtom atom1 = null;
	//		IAtom atom2 = null;
	//
	//		if (a1.contains(b1.getAtom(0)))
	//		{
	//			atom1 = b1.getAtom(0);
	//		}
	//		else if (a1.contains(b1.getAtom(1)))
	//		{
	//			atom1 = b1.getAtom(1);
	//		}
	//
	//		if (a2.contains(b2.getAtom(0)))
	//		{
	//			atom2 = b2.getAtom(0);
	//		}
	//		else if (a2.contains(b2.getAtom(1)))
	//		{
	//			atom2 = b2.getAtom(1);
	//		}
	//
	//		if (atom1 != null && atom2 != null)
	//		{
	//			// well, this looks fishy: the atom2 is not always a IQueryAtom !
	//			return ((IQueryAtom) atom2).matches(atom1);
	//		}
	//		else
	//			return atom1 == null && atom2 == null;
	//
	//	}

	/**
	 *  Determines if 2 bond have 1 atom in common if second is a query AtomContainer
	 *  and whether the order of the atoms is correct (atoms match).
	 *
	 * @param  bond1  first bond
	 * @param  bond2  second bond
	 * @param queryBond1 first query bond
	 * @param queryBond2 second query bond
	 * @return    the symbol of the common atom or "" if the 2 bonds have no common atom
	 */
	private static boolean queryAdjacencyAndOrder(IBond bond1, IBond bond2, IBond queryBond1, IBond queryBond2)
	{

		IAtom centralAtom = null;
		IAtom centralQueryAtom = null;

		if (bond1.contains(bond2.getAtom(0)))
		{
			centralAtom = bond2.getAtom(0);
		}
		else if (bond1.contains(bond2.getAtom(1)))
		{
			centralAtom = bond2.getAtom(1);
		}

		if (queryBond1.contains(queryBond2.getAtom(0)))
		{
			centralQueryAtom = queryBond2.getAtom(0);
		}
		else if (queryBond1.contains(queryBond2.getAtom(1)))
		{
			centralQueryAtom = queryBond2.getAtom(1);
		}

		if (centralAtom != null && centralQueryAtom != null && ((IQueryAtom) centralQueryAtom).matches(centralAtom))
		{
			IQueryAtom queryAtom1 = (IQueryAtom) queryBond1.getConnectedAtom(centralQueryAtom);
			IQueryAtom queryAtom2 = (IQueryAtom) queryBond2.getConnectedAtom(centralQueryAtom);
			IAtom atom1 = bond1.getConnectedAtom(centralAtom);
			IAtom atom2 = bond2.getConnectedAtom(centralAtom);
			if (queryAtom1.matches(atom1) && queryAtom2.matches(atom2) || queryAtom1.matches(atom2)
					&& queryAtom2.matches(atom1))
			{
				return true;
			}
			else
				return false;
		}
		else
			return centralAtom == null && centralQueryAtom == null;

	}

	/**
	 *  Checks some simple heuristics for whether the subgraph query can
	 *  realistically be a subgraph of the supergraph. If, for example, the
	 *  number of nitrogen atoms in the query is larger than that of the supergraph
	 *  it cannot be part of it.
	 *
	 * @param  ac1  the supergraph to be checked. Must not be an {@link IQueryAtomContainer}.
	 * @param  ac2  the subgraph to be tested for. May be an {@link IQueryAtomContainer}.
	 * @return    true if the subgraph ac2 has a chance to be a subgraph of ac1
	 * @throws CDKException if the first molecule is an instance of {@link IQueryAtomContainer}
	 */
	private static boolean testSubgraphHeuristics(IAtomContainer ac1, IAtomContainer ac2) throws CDKException
	{
		if (ac1 instanceof IQueryAtomContainer)
			throw new CDKException("The first IAtomContainer must not be an IQueryAtomContainer");

		int ac1SingleBondCount = 0;
		int ac1DoubleBondCount = 0;
		int ac1TripleBondCount = 0;
		int ac1AromaticBondCount = 0;
		int ac2SingleBondCount = 0;
		int ac2DoubleBondCount = 0;
		int ac2TripleBondCount = 0;
		int ac2AromaticBondCount = 0;
		int ac1SCount = 0;
		int ac1OCount = 0;
		int ac1NCount = 0;
		int ac1FCount = 0;
		int ac1ClCount = 0;
		int ac1BrCount = 0;
		int ac1ICount = 0;
		int ac1CCount = 0;

		int ac2SCount = 0;
		int ac2OCount = 0;
		int ac2NCount = 0;
		int ac2FCount = 0;
		int ac2ClCount = 0;
		int ac2BrCount = 0;
		int ac2ICount = 0;
		int ac2CCount = 0;

		IBond bond;
		IAtom atom;
		for (int i = 0; i < ac1.getBondCount(); i++)
		{
			bond = ac1.getBond(i);
			if (bond.getFlag(CDKConstants.ISAROMATIC))
				ac1AromaticBondCount++;
			else if (bond.getOrder() == IBond.Order.SINGLE)
				ac1SingleBondCount++;
			else if (bond.getOrder() == IBond.Order.DOUBLE)
				ac1DoubleBondCount++;
			else if (bond.getOrder() == IBond.Order.TRIPLE)
				ac1TripleBondCount++;
		}
		for (int i = 0; i < ac2.getBondCount(); i++)
		{
			bond = ac2.getBond(i);
			if (bond instanceof IQueryBond)
				continue;
			if (bond.getFlag(CDKConstants.ISAROMATIC))
				ac2AromaticBondCount++;
			else if (bond.getOrder() == IBond.Order.SINGLE)
				ac2SingleBondCount++;
			else if (bond.getOrder() == IBond.Order.DOUBLE)
				ac2DoubleBondCount++;
			else if (bond.getOrder() == IBond.Order.TRIPLE)
				ac2TripleBondCount++;
		}

		if (ac2SingleBondCount > ac1SingleBondCount)
			return false;
		if (ac2AromaticBondCount > ac1AromaticBondCount)
			return false;
		if (ac2DoubleBondCount > ac1DoubleBondCount)
			return false;
		if (ac2TripleBondCount > ac1TripleBondCount)
			return false;

		for (int i = 0; i < ac1.getAtomCount(); i++)
		{
			atom = ac1.getAtom(i);
			if (atom.getSymbol().equals("S"))
				ac1SCount++;
			else if (atom.getSymbol().equals("N"))
				ac1NCount++;
			else if (atom.getSymbol().equals("O"))
				ac1OCount++;
			else if (atom.getSymbol().equals("F"))
				ac1FCount++;
			else if (atom.getSymbol().equals("Cl"))
				ac1ClCount++;
			else if (atom.getSymbol().equals("Br"))
				ac1BrCount++;
			else if (atom.getSymbol().equals("I"))
				ac1ICount++;
			else if (atom.getSymbol().equals("C"))
				ac1CCount++;
		}
		for (int i = 0; i < ac2.getAtomCount(); i++)
		{
			atom = ac2.getAtom(i);
			if (atom instanceof IQueryAtom)
				continue;
			if (atom.getSymbol().equals("S"))
				ac2SCount++;
			else if (atom.getSymbol().equals("N"))
				ac2NCount++;
			else if (atom.getSymbol().equals("O"))
				ac2OCount++;
			else if (atom.getSymbol().equals("F"))
				ac2FCount++;
			else if (atom.getSymbol().equals("Cl"))
				ac2ClCount++;
			else if (atom.getSymbol().equals("Br"))
				ac2BrCount++;
			else if (atom.getSymbol().equals("I"))
				ac2ICount++;
			else if (atom.getSymbol().equals("C"))
				ac2CCount++;
		}

		if (ac1SCount < ac2SCount)
			return false;
		if (ac1NCount < ac2NCount)
			return false;
		if (ac1OCount < ac2OCount)
			return false;
		if (ac1FCount < ac2FCount)
			return false;
		if (ac1ClCount < ac2ClCount)
			return false;
		if (ac1BrCount < ac2BrCount)
			return false;
		if (ac1ICount < ac2ICount)
			return false;
		return ac1CCount >= ac2CCount;

	}

	/**
	 * Sets the time in milliseconds until the substructure search will be breaked.
	 * @param timeout
	 * Time in milliseconds. -1 to ignore the timeout.
	 */
	public static void setTimeout(long timeout)
	{
		MyUniversalIsomorphismTester.timeout = timeout;
	}

}
