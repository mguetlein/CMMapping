package org.chesmapper.map.alg.build3d;

import java.io.FileOutputStream;

import org.chesmapper.map.data.DatasetFile;
import org.chesmapper.map.data.FeatureService;
import org.chesmapper.map.main.Settings;
import org.mg.javalib.gui.property.Property;
import org.mg.javalib.gui.property.SelectProperty;
import org.mg.javalib.util.StringUtil;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.modeling.builder3d.TemplateHandler3D;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class CDK3DBuilder extends AbstractReal3DBuilder
{
	public static final String[] FORCEFIELDS = { "mm2", "mmff94" };
	public static final CDK3DBuilder INSTANCE = new CDK3DBuilder();

	private CDK3DBuilder()
	{
	}

	SelectProperty forcefield = new SelectProperty("forcefield", FORCEFIELDS, FORCEFIELDS[0]);

	@Override
	public Property[] getProperties()
	{
		return new Property[] { forcefield };
	}

	@Override
	public boolean[] build3D(DatasetFile dataset, String outfile)
	{
		return FeatureService.generateCDK3D(dataset, outfile, forcefield.getValue().toString());
	}

	@Override
	public String getName()
	{
		return Settings.text("build3d.cdk");
	}

	@Override
	public String getDescription()
	{
		return Settings.text("build3d.cdk.desc", Settings.CDK_STRING);
	}

	public static void main(String[] args) throws Exception
	{
		//CN(C)C(=S)S[Zn]SC(=S)N(C)C
		String smiles = "Br(=O)(=O)[O-]";
		IAtomContainer mol = new SmilesParser(SilentChemObjectBuilder.getInstance())
				.parseSmiles(smiles);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		ModelBuilder3D mb3d = ModelBuilder3D.getInstance(TemplateHandler3D.getInstance(), "mmff94",
				DefaultChemObjectBuilder.getInstance());
		IAtomContainer mol3d = mb3d.generate3DCoordinates(mol, true);
		SDFWriter writer = new SDFWriter(new FileOutputStream(
				"/home/martin/.ches-mapper/babel3d/2.3.2/smi/" + StringUtil.getMD5(smiles), true));
		writer.write(mol3d);
		writer.close();
	}
}
