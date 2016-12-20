package org.chesmapper.map.weka;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.chesmapper.map.data.DatasetFile;
import org.chesmapper.map.dataInterface.CompoundData;
import org.chesmapper.map.dataInterface.CompoundProperty;
import org.chesmapper.map.dataInterface.NominalProperty;
import org.chesmapper.map.dataInterface.NumericProperty;
import org.chesmapper.map.main.Settings;
import org.chesmapper.map.main.TaskProvider;
import org.mg.wekalib.data.ArffWritable;
import org.mg.wekalib.data.ArffWriter;

public class CompoundArffWriter implements ArffWritable
{
	public static File writeArffFile(DatasetFile dataset, List<CompoundData> compounds,
			List<CompoundProperty> features)
	{
		return writeArffFile(dataset.getFeatureTableFilePath("arff"), compounds, features);
	}

	public static File writeArffFile(String arffFile, List<CompoundData> compounds,
			List<CompoundProperty> features)
	{
		File file = new File(arffFile);
		if (!Settings.CACHING_ENABLED || !file.exists())
		{
			TaskProvider.debug("writing arff file: " + arffFile);
			try
			{
				ArffWriter.writeToArffFile(file, new CompoundArffWriter(compounds, features));
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}
		}
		else
			TaskProvider.debug("arff file already exists: " + arffFile);
		return file;
	}

	List<CompoundData> compounds;
	List<CompoundProperty> features;
	boolean sparse = true;

	private CompoundArffWriter(List<CompoundData> compounds, List<CompoundProperty> features)
	{
		this.compounds = compounds;
		this.features = features;

		for (CompoundProperty p : features)
			if (p instanceof NumericProperty)
			{
				sparse = false;
				break;
			}
		if (features.size() < 100)
			sparse = false;
	}

	@Override
	public String getRelationName()
	{
		return "CheS-Mapper-Dataset";
	}

	@Override
	public List<String> getAdditionalInfo()
	{
		return null;
	}

	@Override
	public int getNumAttributes()
	{
		return features.size();
	}

	@Override
	public String getAttributeName(int attribute)
	{
		//return features.get(attribute).getUniqueName();
		return features.get(attribute).getName() + "_" + features.get(attribute).hashCode();
	}

	@Override
	public String[] getAttributeDomain(int attribute)
	{
		if (features.get(attribute) instanceof NumericProperty)
			return null;
		else
			return ((NominalProperty) features.get(attribute)).getDomain();
	}

	@Override
	public int getNumInstances()
	{
		return compounds.size();
	}

	Map<Integer, Map<String, Integer>> nominalFeatureMap = new HashMap<>();

	@Override
	public double getAttributeValueAsDouble(int instance, int attribute) throws Exception
	{
		if (features.get(attribute) instanceof NumericProperty)
		{
			Double v = compounds.get(instance)
					.getNormalizedValueCompleteDataset((NumericProperty) features.get(attribute));
			if (v == null)
				return Double.NaN;
			else
				return v;
		}
		else
		{
			String s = compounds.get(instance)
					.getStringValue((NominalProperty) features.get(attribute));
			if (s == null)
				return Double.NaN;
			if (!nominalFeatureMap.containsKey(attribute))
			{
				Map<String, Integer> map = new HashMap<>();
				int i = 0;
				for (String v : getAttributeDomain(attribute))
					map.put(v, i++);
				nominalFeatureMap.put(attribute, map);
			}
			return nominalFeatureMap.get(attribute).get(s);
		}
	}

	@Override
	public String getAttributeValue(int instance, int attribute)
	{
		if (features.get(attribute) instanceof NumericProperty)
		{
			Double v = compounds.get(instance)
					.getNormalizedValueCompleteDataset((NumericProperty) features.get(attribute));
			if (v == null)
				return "?";
			else
				return v.toString();
		}
		else
		{
			String s = compounds.get(instance)
					.getStringValue((NominalProperty) features.get(attribute));
			if (s == null)
				return "?";
			else if (s.length() > 1)
				return "\"" + s + "\"";
			else
				return s;
		}
	}

	@Override
	public boolean isSparse()
	{
		return sparse;
	}

	@Override
	public String getMissingValue(int attribute)
	{
		return "?";
	}

}
