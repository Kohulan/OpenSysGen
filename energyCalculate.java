//Use this code in order to check the energy of a small molecule using different forcefields

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import org.openbabel.OBConversion;
import org.openbabel.OBForceField;
import org.openbabel.OBMol;

public class energyCalculate {
	//use for .pdb files as default
	public static int x; // file counter

	public static void main(String[] args) throws IOException {
		System.loadLibrary("openbabel_java");

		File folder = new File("input folder path");
		File[] listOfFiles = folder.listFiles();
		File fileout = new File("output folder path/out.csv");
		if (!fileout.exists()) {
			fileout.createNewFile();
		}

		for (x = 0; x < listOfFiles.length; x++) {
			File file = listOfFiles[x];
			if (file.isFile() && file.getName().endsWith(".pdb")) { //change .pdb if you are using a different file format
				InputStream is = new FileInputStream(file);
				@SuppressWarnings("resource")
				BufferedReader buf = new BufferedReader(new InputStreamReader(is));
				String line = buf.readLine();
				StringBuilder sb = new StringBuilder();
				while (line != null) {
					sb.append(line).append("\n");
					line = buf.readLine();
				}
				String content = sb.toString();

				OBConversion conv = new OBConversion();
				OBMol mol = new OBMol();
				conv.SetInFormat("pdb"); //change pdb if you are using a different file format, check OpenBabel documentation for supported file formats
				conv.ReadString(mol, content);
				OBForceField force = OBForceField.FindForceField("mmff94"); //ghemical,gaff and uff also works. For more information check OpenBabel
				long rotnum = mol.NumRotors();
				force.Setup(mol);
				double energy = force.Energy();
				String energycontent = file.getName()+","+rotnum + "," + energy + "\n";

				// if file doesnt exists, then create it

				System.out.println(energycontent);

				FileWriter fw = new FileWriter(fileout.getAbsoluteFile(), true);
				BufferedWriter bufferWriter = new BufferedWriter(fw);
				bufferWriter.write(energycontent);
				bufferWriter.close();

			}

		}

	}

}

