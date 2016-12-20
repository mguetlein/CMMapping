package org.chesmapper.map.property;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.chesmapper.map.data.DatasetFile;
import org.chesmapper.map.data.fragments.MatchEngine;
import org.chesmapper.map.dataInterface.DefaultFragmentProperty;
import org.chesmapper.map.dataInterface.FragmentPropertySet;
import org.chesmapper.map.main.BinHandler;
import org.chesmapper.map.main.Settings;
import org.chesmapper.map.main.TaskProvider;
import org.chesmapper.map.util.ExternalToolUtil;
import org.mg.javalib.util.ArrayUtil;
import org.mg.javalib.util.FileUtil;
import org.openscience.cdk.fingerprint.IBitFingerprint;

public class OBMolPrintFingerprintSet extends FragmentPropertySet // AbstractPropertySet
{
	//	int minF = 3;

	public static OBMolPrintFingerprintSet[] FINGERPRINTS = null;

	static
	{
		try
		{
			FINGERPRINTS = new OBMolPrintFingerprintSet[] { new OBMolPrintFingerprintSet() };
		}
		catch (Exception e)
		{
			Settings.LOGGER.error(e);
		}
	}

	//	HashMap<DatasetFile, List<DefaultFragmentProperty>> props = new HashMap<DatasetFile, List<DefaultFragmentProperty>>();

	//	private OBMolPrintFingerprintSet()
	//	{
	//	}

	//	public static OBMolPrintFingerprintSet fromString(String string)
	//	{
	//		for (OBMolPrintFingerprintSet set : FINGERPRINTS)
	//			if (set.toString().equals(string))
	//				return set;
	//		return null;
	//	}
	//
	//	@Override
	//	public boolean equals(Object o)
	//	{
	//		return (o instanceof OBMolPrintFingerprintSet)
	//				&& ((OBMolPrintFingerprintSet) o).toString().equals(toString());
	//	}

	public String toString()
	{
		return Settings.text("features.ob.molprint2d");
	}

	@Override
	public String getDescription()
	{
		return Settings.text("features.ob.molprint2d.desc", Settings.OPENBABEL_STRING);
	}

	@Override
	public boolean compute(DatasetFile dataset)
	{
		try
		{
			List<IBitFingerprint> fingerprints = new ArrayList<IBitFingerprint>();
			TaskProvider.debug("Computing MolPrint2D with OpenBabel");

			//		String smiles = "";
			//		for (String smi : dataset.getSmiles())
			//			smiles += smi+"\n";
			//		File smi = File.createTempFile(dataset.getShortName(), "smi");
			//		FileUtil.writeStringToFile(smi.getAbsolutePath(), smiles);
			File tmp = File.createTempFile(dataset.getShortName(), "OBfingerprint");
			String cmd[] = { BinHandler.BABEL_BINARY.getLocation(), "-isdf", dataset.getSDF(),
					"-ompd" };
			TaskProvider.debug("Running babel: " + ArrayUtil.toString(cmd, " ", "", ""));
			ExternalToolUtil.run("ob-fingerprints", cmd, tmp);

			BufferedReader reader = new BufferedReader(new FileReader(tmp));
			String line = null;
			int molIdx = 0;
			HashMap<String, Set<Integer>> map = new HashMap<>();
			while ((line = reader.readLine()) != null)
			{
				String s[] = line.split("\t");
				if (!s[0].endsWith((molIdx + 1) + ""))
					throw new IllegalStateException(line);
				for (int featureIdx = 1; featureIdx < s.length; featureIdx++)
				{
					String f = s[featureIdx];
					if (!map.containsKey(f))
						map.put(f, new HashSet<Integer>());
					map.get(f).add(molIdx);
				}
				molIdx++;
			}
			reader.close();

			props.put(dataset, new ArrayList<DefaultFragmentProperty>());
			for (String f : map.keySet())
			{
				//				if (map.get(f).size() >= minF)
				//				{
				DefaultFragmentProperty p = new DefaultFragmentProperty(this, renameFeature(f), f,
						"", MatchEngine.OpenBabel);
				String v[] = new String[dataset.numCompounds()];
				for (int i = 0; i < v.length; i++)
					v[i] = map.get(f).contains(i) ? "1" : "0";
				p.setStringValues(v);
				p.setFrequency(map.get(f).size());
				props.get(dataset).add(p);
				//				}
			}
			updateFragments();
			return true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return false;
		}
	}

	@Override
	public boolean hasFixedMatchEngine()
	{
		return true;
	}

	@Override
	public boolean isComputationSlow()
	{
		return false;
	}

	@Override
	public boolean isSizeDynamicHigh(DatasetFile dataset)
	{
		return true;
	}

	private static String molToSym[];

	private static String renameFeature(String f)
	{
		for (int i = f.length() - 1; i >= 0; i--)
		{
			if (f.charAt(i) == ';')
			{
				int begin = -1;
				if (i == 1 || i == 2)
					begin = 0;
				else if (f.charAt(i - 2) == '-')
					begin = i - 1;
				else if (f.charAt(i - 3) == '-')
					begin = i - 2;
				//				System.out.println(begin);
				//				System.out.println(i);
				int mol = Integer.valueOf(f.substring(begin, i));
				//				System.out.println(mol);
				String sym = molToSym[mol];
				//				System.out.println(mcml);
				f = f.substring(0, begin) + sym + f.substring(i);
				//				System.out.println(f);
				int shortenedBy = (i - begin) - sym.length();
				if (shortenedBy > 0)
					i -= shortenedBy;
			}
		}
		f = f.replaceFirst(";1-", " 1:-");
		f = f.replaceFirst(";2-", " 2:-");
		f = f.replaceAll(";1-", ",-");
		f = f.replaceAll(";2-", ",-");
		f = f.replaceAll("-1-", "");
		f = f.replaceAll("-([0-9]+)-", "$1x");
		f = f.replaceAll(";$", "");
		return f;
	}

	static
	{
		String s[] = FileUtil.readStringFromFile("/usr/share/openbabel/2.3.2/types.txt")
				.split("\n");
		int numLines = 0;
		int numCols = 0;
		String h[] = null;
		HashMap<String, String[]> cols = new HashMap<>();
		int lineCount = 0;
		for (String line : s)
		{
			if (line.startsWith("#"))
				continue;
			else if (numLines == 0)
			{
				String l[] = line.split("[\\s]+");
				if (l.length != 2)
					throw new IllegalStateException(ArrayUtil.toString(l));
				numLines = Integer.valueOf(l[0]) + 1;
				numCols = Integer.valueOf(l[1]);
			}
			else if (h == null)
			{
				h = line.split("\\s");
				if (h.length != numCols)
					throw new IllegalStateException();
				for (String col : h)
					cols.put(col, new String[numLines]);
			}
			else
			{
				//				System.out.println(lineCount + ": " + line);
				int col = 0;
				for (String v : line.split("[\\s]+"))
				{
					cols.get(h[col])[lineCount] = v;
					col++;
				}
				lineCount++;
			}
		}
		//		System.out.println(numLines + " x " + numCols);
		//		System.out.println(ArrayUtil.toString(h));
		//		System.out.println(ArrayUtil.toString(cols.get(h[0])));

		int maxMol = 0;
		for (String mol : cols.get("MOL"))
			maxMol = Math.max(maxMol, Integer.valueOf(mol));

		molToSym = new String[maxMol + 1];

		int row = 0;
		for (String mol : cols.get("MOL"))
		{
			int m = Integer.valueOf(mol);
			if (m != 0)
			{
				String sym = cols.get("SYB")[row];
				//				if (m == 1)
				//					System.out.println("[" + row + "] " + m + ": " + mcml);

				String sy = molToSym[m];
				if (sy == null)
					molToSym[m] = sym;
				else if (!sy.equals(sym) && !sy.endsWith(sym) && !sy.contains(sym + "/"))
					molToSym[m] += "/" + sym;
			}
			row++;
		}

		for (int i = 1; i < molToSym.length; i++)
		{
			if (molToSym[i] == null)
				continue;
			if (molToSym[i].equals("C.3/C/C."))
				molToSym[i] = "C";
			if (molToSym[i].equals("N.2/N.4/N.pl3/N"))
				molToSym[i] = "N";
			if (molToSym[i].equals("N.2/N.pl3/N.pl3N"))
				molToSym[i] = "N";
			if (molToSym[i].equals("O.2/O.co2/O.co2O/O."))
				molToSym[i] = "O";
		}

		for (int i = 1; i < molToSym.length; i++)
			System.out.println(i + ": " + molToSym[i]);
	}

	public static void main(String[] args)
	{
		System.out.println(renameFeature("11;1-1-1;1-2-3;2-1-1;2-1-3;2-1-8;2-1-9;2-1-11;"));
	}

}