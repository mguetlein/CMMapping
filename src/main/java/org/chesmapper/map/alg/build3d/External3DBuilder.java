package org.chesmapper.map.alg.build3d;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;

import org.chesmapper.map.data.DatasetFile;
import org.chesmapper.map.main.Settings;
import org.chesmapper.map.main.TaskProvider;
import org.mg.javalib.io.External3DComputer;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

public class External3DBuilder extends AbstractReal3DBuilder
{
	public static final External3DBuilder INSTANCE = new External3DBuilder();

	@Override
	public String getName()
	{
		return External3DBuilder.class.getName();
	}

	@Override
	public String getDescription()
	{
		return "using external script";
	}

	@Override
	public boolean[] build3D(DatasetFile dataset, String outFile) throws Exception
	{
		boolean valid[] = new boolean[dataset.getCompounds().length];
		try
		{
			if (new File(outFile).exists())
				new File(outFile).delete();

			IAtomContainer mols[] = dataset.getCompounds();

			int count = 0;
			for (IAtomContainer iMolecule : mols)
			{
				IAtomContainer molecule = iMolecule;

				String smiles = dataset.getSmiles()[count];
				System.out.println(smiles);

				String molFormat = External3DComputer.get3D(smiles);

				ISimpleChemObjectReader reader = new MDLV2000Reader(
						new ByteArrayInputStream(molFormat.getBytes(StandardCharsets.UTF_8)));
				IChemFile content = (IChemFile) reader.read((IChemObject) new ChemFile());
				IAtomContainer mol = (IAtomContainer) ChemFileManipulator
						.getAllAtomContainers(content).get(0);
				reader.close();
				mol = (IAtomContainer) AtomContainerManipulator.removeHydrogens(mol);

				mol.setProperties(molecule.getProperties());

				SDFWriter writer = new SDFWriter(new FileOutputStream(outFile, true));
				writer.write(molecule);
				writer.close();
				count++;

				if (!TaskProvider.isRunning())
					return null;
			}

		}
		catch (CDKException e)
		{
			Settings.LOGGER.error(e);
		}
		catch (IOException e)
		{
			Settings.LOGGER.error(e);
		}
		return valid;
	}

}
