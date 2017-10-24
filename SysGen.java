import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Scanner;

import org.openbabel.OBAlign;
import org.openbabel.OBAtom;
import org.openbabel.OBBond;
import org.openbabel.OBConversion;
import org.openbabel.OBForceField;
import org.openbabel.OBMol;
import org.openbabel.OBMolBondIter;

public class SysGen {

	

	public static void main(String[] args) throws IOException {
		System.out.println("********************************************************************************\n"
				+ "Welcome to OpenSysGen : An Opensource tool for Systematic Conformer Generation\n"
				+ "                 Written By : Kohulan Rajan & Pranav Khade\n"
				+ "                        Published on : 2017/06/12\n"
				+ "                         Version : V 1.1 (stable)\n"
				+"        ©Bioinformatics Centre,Savithribai Phule Pune University,India\n"
				+ "********************************************************************************");
		Scanner reader = new Scanner(System.in);
		System.out.println("Enter the file type (Use only: pdb,mol2,sdf,mol): ");
		String type = reader.next(); // Takes user input file type
		System.out.println("Enter the required file output type (Use only: pdb,mol2,sdf,mol): ");
		String outtype = reader.next();
		reader.close();
		int x;
		File report= new File("Final_Report.tsv");
		if (!report.exists()) {
			report.createNewFile();
		}
		FileWriter fre = new FileWriter(report.getAbsoluteFile(), true);
		BufferedWriter bufferfre = new BufferedWriter(fre);
		bufferfre.write("File name\tRotatable Bond Count\tTotal Conformers\tTotal conformers <1.0 RMSD\n");
		for (x = 0; x < args.length; x++) {
			long startTime = System.nanoTime(); // Setting start time for each file
			File fileout = null;
			String filename = args[x];
			filename = filename.substring(0, filename.length() - 4);

			System.loadLibrary("openbabel_java");
			// Converting molecule to OpenBabel readable format
			OBConversion conv = new OBConversion();
			OBMol mol = new OBMol();
			OBMol mol_ref = new OBMol();
			ArrayList<OBAtom> rotbonds = new ArrayList<OBAtom>();
			conv.SetInFormat(type);
			conv.ReadFile(mol_ref, args[x]);
			conv.ReadFile(mol, args[x]);

			// Extracting rotatable bonds

			for (OBBond bond : new OBMolBondIter(mol)) {
				if (bond.IsRotor()) {
					OBAtom a2 = bond.GetBeginAtom();
					OBAtom a3 = bond.GetEndAtom();
					OBAtom a1 = a2.NextNbrAtom(a2.BeginBonds());
					OBAtom a4 = a3.NextNbrAtom(a3.BeginBonds());
					rotbonds.add(a1);
					rotbonds.add(a2);
					rotbonds.add(a3);
					rotbonds.add(a4);
				}
			}

			long rotnum = mol.NumRotors(); // Getting the number of available rotatable bonds
			if(rotnum <= 10){
				if (type.equals("pdb")) { // Checking file type pdb or not
					fileout = new File(filename + "__RMSD_Report" + ".csv");
					
					if (!fileout.exists()) {
						fileout.createNewFile();
						
					}
				} else if (type.equals("sdf") || type.equals("mol") || type.equals("mol2")) {
					fileout = new File(filename + "__RMSD_Report" + ".csv");
					
					if (!fileout.exists() && !report.exists()) {
						fileout.createNewFile();
						
					}
				} else {
					System.out.println("Error! Please Check File Type");
					bufferfre.write("Error! Please Check File Type\n");
					break;
				}
			}
			// Rotatable Bonds more than 10
			else{
				System.out.println("Sorry!!!\n The molecule you have used contains " + rotnum
						+ " rotatable bonds, Please use small molecules contains less than 10 Rotatable Bonds\n\n");
				bufferfre.write("Sorry!!!\n The molecule you have used contains " + rotnum
						+ " rotatable bonds, Please use small molecules contains less than 10 Rotatable Bonds\n\n");
				break;
			}
			
			OBAlign align = new OBAlign();
			OBForceField force = OBForceField.FindForceField("mmff94"); // Setting up the force field type
			conv.SetOutFormat(outtype); // Setting output type
			@SuppressWarnings("unused")
			boolean append = conv.WriteFile(mol, filename + "__conformer" + "." + outtype); // Creating empty conformer file for appending

			FileWriter fw = new FileWriter(fileout.getAbsoluteFile(), true);
			BufferedWriter bufferWriter = new BufferedWriter(fw);
			
			if(type.equals("pdb")){
			bufferWriter.write("Conformer No ,RMSD,Energy\n");
			
			}else{
				bufferWriter.write("Conformer No ,RMSD\n");
				
			}
			int counter = 1;
			int count_rmsd=0;
			// Rotatable Bonds 1
			if (rotnum == 1) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bond");
				
				for (double i = 10; i <= 360.0; i = i + 10.0) {
					mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2), rotbonds.get(3),
							i * (3.14159265358979323846 / 180));
					align.SetRefMol(mol_ref);
					align.SetTargetMol(mol);
					align.Align();
					force.Setup(mol);
					double energy = force.Energy();
					double rmsd = align.GetRMSD();
					if (type.equals("pdb")) {
						bufferWriter.write( counter + "," + rmsd + "," + energy + "\n");
						append = conv.Write(mol);
					} else {
						bufferWriter.write( counter + "," + rmsd + "\n");
						append = conv.Write(mol);
					}
					if(rmsd<=1.0){
						count_rmsd++;
					}
					counter++;
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
			}
			
			// Rotatable Bonds 2
			else if (rotnum == 2) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bonds");
				
				for (double i = 10; i <= 360.0; i = i + 10.0) {
					for (double j = 10; j <= 360.0; j = j + 10.0) {

						mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2), rotbonds.get(3),
								i * (3.14159265358979323846 / 180));
						mol.SetTorsion(rotbonds.get(4), rotbonds.get(5), rotbonds.get(6), rotbonds.get(7),
								j * (3.14159265358979323846 / 180));
						align.SetRefMol(mol_ref);
						align.SetTargetMol(mol);
						align.Align();
						force.Setup(mol);
						double energy = force.Energy();
						double rmsd = align.GetRMSD();
						if (type.equals("pdb")) {
							bufferWriter.write(counter + "," + rmsd + "," + energy + "\n");
							append = conv.Write(mol);
						} else {
							bufferWriter.write(counter + "," + rmsd + "\n");
							append = conv.Write(mol);
						}
						if(rmsd<=1.0){
							count_rmsd++;
						}
						counter++;
					}
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
				}
			
			// Rotatable Bonds 3
			else if (rotnum == 3) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bonds");
				
				for (double i = 10; i <= 360.0; i = i + 10.0) {
					for (double j = 10; j <= 360.0; j = j + 10.0) {
						for (double k = 10; k <= 360.0; k = k + 10.0) {

							mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2), rotbonds.get(3),
									i * (3.14159265358979323846 / 180));
							mol.SetTorsion(rotbonds.get(4), rotbonds.get(5), rotbonds.get(6), rotbonds.get(7),
									j * (3.14159265358979323846 / 180));
							mol.SetTorsion(rotbonds.get(8), rotbonds.get(9), rotbonds.get(10), rotbonds.get(11),
									k * (3.14159265358979323846 / 180));
							align.SetRefMol(mol_ref);
							align.SetTargetMol(mol);
							align.Align();
							force.Setup(mol);
							double energy = force.Energy();
							double rmsd = align.GetRMSD();
							if (type.equals("pdb")) {
								bufferWriter.write(counter + "," + rmsd + "," + energy + "\n");
								append = conv.Write(mol);
							} else {
								bufferWriter.write(counter + "," + rmsd + "\n");
								append = conv.Write(mol);
							}
							if(rmsd<=1.0){
								count_rmsd++;
							}
							counter++;
						}
					}
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
			}
			
			// Rotatable Bonds 4
			else if (rotnum == 4) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bonds");
				
				for (double i = 30; i <= 360.0; i = i + 30.0) {
					for (double j = 30; j <= 360.0; j = j + 30.0) {
						for (double k = 30; k <= 360.0; k = k + 30.0) {
							for (double l = 30; l <= 360.0; l = l + 30.0) {
								mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2), rotbonds.get(3),
										i * (3.14159265358979323846 / 180));
								mol.SetTorsion(rotbonds.get(4), rotbonds.get(5), rotbonds.get(6), rotbonds.get(7),
										j * (3.14159265358979323846 / 180));
								mol.SetTorsion(rotbonds.get(8), rotbonds.get(9), rotbonds.get(10), rotbonds.get(11),
										k * (3.14159265358979323846 / 180));
								mol.SetTorsion(rotbonds.get(12), rotbonds.get(13), rotbonds.get(14), rotbonds.get(15),
										l * (3.14159265358979323846 / 180));
								align.SetRefMol(mol_ref);
								align.SetTargetMol(mol);
								align.Align();
								force.Setup(mol);
								double energy = force.Energy();
								double rmsd = align.GetRMSD();
								if (type.equals("pdb")) {
									bufferWriter.write(counter + "," + rmsd + "," + energy + "\n");
									append = conv.Write(mol);
								} else {
									bufferWriter.write(counter + "," + rmsd + "\n");
									append = conv.Write(mol);
								}
								if(rmsd<=1.0){
									count_rmsd++;
								}
								counter++;
							}
						}
					}
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
			}
			
			// Rotatable Bonds 5
			else if (rotnum == 5) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bonds");
				
				for (double i = 30; i <= 360.0; i = i + 30.0) {
					for (double j = 30; j <= 360.0; j = j + 30.0) {
						for (double k = 30; k <= 360.0; k = k + 30.0) {
							for (double l = 30; l <= 360.0; l = l + 30.0) {
								for (double m = 30; m <= 360.0; m = m + 30.0) {
									mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2), rotbonds.get(3),
											i * (3.14159265358979323846 / 180));
									mol.SetTorsion(rotbonds.get(4), rotbonds.get(5), rotbonds.get(6), rotbonds.get(7),
											j * (3.14159265358979323846 / 180));
									mol.SetTorsion(rotbonds.get(8), rotbonds.get(9), rotbonds.get(10), rotbonds.get(11),
											k * (3.14159265358979323846 / 180));
									mol.SetTorsion(rotbonds.get(12), rotbonds.get(13), rotbonds.get(14), rotbonds.get(15),
											l * (3.14159265358979323846 / 180));
									mol.SetTorsion(rotbonds.get(16), rotbonds.get(17), rotbonds.get(18), rotbonds.get(19),
											m * (3.14159265358979323846 / 180));
									align.SetTargetMol(mol);
									align.SetRefMol(mol_ref);
									align.SetTargetMol(mol);
									align.Align();
									force.Setup(mol);
									double energy = force.Energy();
									double rmsd = align.GetRMSD();
									if (type.equals("pdb")) {
										bufferWriter
												.write(counter + "," + rmsd + "," + energy + "\n");
										append = conv.Write(mol);
									} else {
										bufferWriter.write(counter + "," + rmsd + "\n");
										append = conv.Write(mol);
									}
									if(rmsd<=1.0){
										count_rmsd++;
									}
									counter++;
								}
							}
						}
					}
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
			}
			
			// Rotatable Bonds 6
			else if (rotnum == 6) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bonds");
				
				for (double i = 45; i <= 360.0; i = i + 45.0) {
					for (double j = 45; j <= 360.0; j = j + 45.0) {
						for (double k = 45; k <= 360.0; k = k + 45.0) {
							for (double l = 45; l <= 360.0; l = l + 45.0) {
								for (double m = 45; m <= 360.0; m = m + 45.0) {
									for (double n = 45; n <= 360.0; n = n + 45.0) {
										mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2),
												rotbonds.get(3), i * (3.14159265358979323846 / 180));
										mol.SetTorsion(rotbonds.get(4), rotbonds.get(5), rotbonds.get(6),
												rotbonds.get(7), j * (3.14159265358979323846 / 180));
										mol.SetTorsion(rotbonds.get(8), rotbonds.get(9), rotbonds.get(10),
												rotbonds.get(11), k * (3.14159265358979323846 / 180));
										mol.SetTorsion(rotbonds.get(12), rotbonds.get(13), rotbonds.get(14),
												rotbonds.get(15), l * (3.14159265358979323846 / 180));
										mol.SetTorsion(rotbonds.get(16), rotbonds.get(17), rotbonds.get(18),
												rotbonds.get(19), m * (3.14159265358979323846 / 180));
										mol.SetTorsion(rotbonds.get(20), rotbonds.get(21), rotbonds.get(22),
												rotbonds.get(23), n * (3.14159265358979323846 / 180));
										align.SetTargetMol(mol);
										align.SetRefMol(mol_ref);
										align.SetTargetMol(mol);
										align.Align();
										force.Setup(mol);
										double energy = force.Energy();
										double rmsd = align.GetRMSD();
										if (type.equals("pdb")) {
											bufferWriter.write(
													counter + "," + rmsd + "," + energy + "\n");
											append = conv.Write(mol);
										} else {
											bufferWriter.write(counter + "," + rmsd + "\n");
											append = conv.Write(mol);
										}
										if(rmsd<=1.0){
											count_rmsd++;
										}
										counter++;
									}
								}
							}
						}
					}
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
			}
			
			// Rotatable Bonds 7
			else if (rotnum == 7) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bonds");
				
				for (double i = 60; i <= 360.0; i = i + 60.0) {
					for (double j = 60; j <= 360.0; j = j + 60.0) {
						for (double k = 60; k <= 360.0; k = k + 60.0) {
							for (double l = 60; l <= 360.0; l = l + 60.0) {
								for (double m = 60; m <= 360.0; m = m + 60.0) {
									for (double n = 60; n <= 360.0; n = n + 60.0) {
										for (double o = 60; o <= 360.0; o = o + 60.0) {
											mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2),
													rotbonds.get(3), i * (3.14159265358979323846 / 180));
											mol.SetTorsion(rotbonds.get(4), rotbonds.get(5), rotbonds.get(6),
													rotbonds.get(7), j * (3.14159265358979323846 / 180));
											mol.SetTorsion(rotbonds.get(8), rotbonds.get(9), rotbonds.get(10),
													rotbonds.get(11), k * (3.14159265358979323846 / 180));
											mol.SetTorsion(rotbonds.get(12), rotbonds.get(13), rotbonds.get(14),
													rotbonds.get(15), l * (3.14159265358979323846 / 180));
											mol.SetTorsion(rotbonds.get(16), rotbonds.get(17), rotbonds.get(18),
													rotbonds.get(19), m * (3.14159265358979323846 / 180));
											mol.SetTorsion(rotbonds.get(20), rotbonds.get(21), rotbonds.get(22),
													rotbonds.get(23), n * (3.14159265358979323846 / 180));
											mol.SetTorsion(rotbonds.get(24), rotbonds.get(25), rotbonds.get(26),
													rotbonds.get(27), o * (3.14159265358979323846 / 180));
											align.SetTargetMol(mol);
											align.SetRefMol(mol_ref);
											align.SetTargetMol(mol);
											align.Align();
											force.Setup(mol);
											double energy = force.Energy();
											double rmsd = align.GetRMSD();
											if (type.equals("pdb")) {
												bufferWriter.write(
														counter + "," + rmsd + "," + energy + "\n");
												append = conv.Write(mol);
											} else {
												bufferWriter.write(counter + "," + rmsd + "\n");
												append = conv.Write(mol);
											}
											if(rmsd<=1.0){
												count_rmsd++;
											}
											counter++;
										}
									}
								}
							}
						}
					}
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
			}
			
			// Rotatable Bonds 8
			else if (rotnum == 8) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bonds");
				
				for (double i = 60; i <= 360.0; i = i + 60.0) {
					for (double j = 60; j <= 360.0; j = j + 60.0) {
						for (double k = 60; k <= 360.0; k = k + 60.0) {
							for (double l = 60; l <= 360.0; l = l + 60.0) {
								for (double m = 60; m <= 360.0; m = m + 60.0) {
									for (double n = 60; n <= 360.0; n = n + 60.0) {
										for (double o = 60; o <= 360.0; o = o + 60.0) {
											for (double p = 60; p <= 360.0; p = p + 60.0) {
												mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2),
														rotbonds.get(3), i * (3.14159265358979323846 / 180));
												mol.SetTorsion(rotbonds.get(4), rotbonds.get(5), rotbonds.get(6),
														rotbonds.get(7), j * (3.14159265358979323846 / 180));
												mol.SetTorsion(rotbonds.get(8), rotbonds.get(9), rotbonds.get(10),
														rotbonds.get(11), k * (3.14159265358979323846 / 180));
												mol.SetTorsion(rotbonds.get(12), rotbonds.get(13), rotbonds.get(14),
														rotbonds.get(15), l * (3.14159265358979323846 / 180));
												mol.SetTorsion(rotbonds.get(16), rotbonds.get(17), rotbonds.get(18),
														rotbonds.get(19), m * (3.14159265358979323846 / 180));
												mol.SetTorsion(rotbonds.get(20), rotbonds.get(21), rotbonds.get(22),
														rotbonds.get(23), n * (3.14159265358979323846 / 180));
												mol.SetTorsion(rotbonds.get(24), rotbonds.get(25), rotbonds.get(26),
														rotbonds.get(27), o * (3.14159265358979323846 / 180));
												mol.SetTorsion(rotbonds.get(28), rotbonds.get(29), rotbonds.get(30),
														rotbonds.get(31), p * (3.14159265358979323846 / 180));
												align.SetTargetMol(mol);
												align.SetRefMol(mol_ref);
												align.SetTargetMol(mol);
												align.Align();
												force.Setup(mol);
												double energy = force.Energy();
												double rmsd = align.GetRMSD();
												if (type.equals("pdb")) {
													bufferWriter.write(counter + "," + rmsd + ","
															+ energy + "\n");
													append = conv.Write(mol);
												} else {
													bufferWriter.write(counter + "," + rmsd + "\n");
													append = conv.Write(mol);
												}
												if(rmsd<=1.0){
													count_rmsd++;
												}
												counter++;
											}
										}
									}
								}
							}
						}
					}
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
			}
			
			// Rotatable Bonds 9
			else if (rotnum == 9) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bonds");
				
				for (double i = 72; i <= 360.0; i = i + 72.0) {
					for (double j = 72; j <= 360.0; j = j + 72.0) {
						for (double k = 72; k <= 360.0; k = k + 72.0) {
							for (double l = 72; l <= 360.0; l = l + 72.0) {
								for (double m = 72; m <= 360.0; m = m + 72.0) {
									for (double n = 72; n <= 360.0; n = n + 72.0) {
										for (double o = 72; o <= 360.0; o = o + 72.0) {
											for (double p = 72; p <= 360.0; p = p + 72.0) {
												for (double q = 72; q <= 360.0; q = q + 72.0) {
													mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2),
															rotbonds.get(3), i * (3.14159265358979323846 / 180));
													mol.SetTorsion(rotbonds.get(4), rotbonds.get(5), rotbonds.get(6),
															rotbonds.get(7), j * (3.14159265358979323846 / 180));
													mol.SetTorsion(rotbonds.get(8), rotbonds.get(9), rotbonds.get(10),
															rotbonds.get(11), k * (3.14159265358979323846 / 180));
													mol.SetTorsion(rotbonds.get(12), rotbonds.get(13), rotbonds.get(14),
															rotbonds.get(15), l * (3.14159265358979323846 / 180));
													mol.SetTorsion(rotbonds.get(16), rotbonds.get(17), rotbonds.get(18),
															rotbonds.get(19), m * (3.14159265358979323846 / 180));
													mol.SetTorsion(rotbonds.get(20), rotbonds.get(21), rotbonds.get(22),
															rotbonds.get(23), n * (3.14159265358979323846 / 180));
													mol.SetTorsion(rotbonds.get(24), rotbonds.get(25), rotbonds.get(26),
															rotbonds.get(27), o * (3.14159265358979323846 / 180));
													mol.SetTorsion(rotbonds.get(28), rotbonds.get(29), rotbonds.get(30),
															rotbonds.get(31), p * (3.14159265358979323846 / 180));
													mol.SetTorsion(rotbonds.get(32), rotbonds.get(33), rotbonds.get(34),
															rotbonds.get(35), q * (3.14159265358979323846 / 180));
													align.SetTargetMol(mol);
													align.SetRefMol(mol_ref);
													align.SetTargetMol(mol);
													align.Align();
													force.Setup(mol);
													double energy = force.Energy();
													double rmsd = align.GetRMSD();
													if (type.equals("pdb")) {
														bufferWriter.write(counter + "," + rmsd
																+ "," + energy + "\n");
														append = conv.Write(mol);
													} else {
														bufferWriter.write(counter + "," + rmsd + "\n");
														append = conv.Write(mol);
													}
													if(rmsd<=1.0){
														count_rmsd++;
													}
													counter++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
			}
			
			// Rotatable Bonds 10
			else if (rotnum == 10) {
				System.out.println("Conformer Generation Started Successfully!!\n\n");
				System.out.println("The Small molecule " + filename + " contains: " + rotnum + " rotatable bonds");
				
				for (double i = 90; i <= 360.0; i = i + 90.0) {
					for (double j = 90; j <= 360.0; j = j + 90.0) {
						for (double k = 90; k <= 360.0; k = k + 90.0) {
							for (double l = 90; l <= 360.0; l = l + 90.0) {
								for (double m = 90; m <= 360.0; m = m + 90.0) {
									for (double n = 90; n <= 360.0; n = n + 90.0) {
										for (double o = 90; o <= 360.0; o = o + 90.0) {
											for (double p = 90; p <= 360.0; p = p + 90.0) {
												for (double q = 90; q <= 360.0; q = q + 90.0) {
													for (double r = 90; r <= 360.0; r = r + 90.0) {
														mol.SetTorsion(rotbonds.get(0), rotbonds.get(1), rotbonds.get(2),
																rotbonds.get(3), i * (3.14159265358979323846 / 180));
														mol.SetTorsion(rotbonds.get(4), rotbonds.get(5), rotbonds.get(6),
																rotbonds.get(7), j * (3.14159265358979323846 / 180));
														mol.SetTorsion(rotbonds.get(8), rotbonds.get(9), rotbonds.get(10),
																rotbonds.get(11), k * (3.14159265358979323846 / 180));
														mol.SetTorsion(rotbonds.get(12), rotbonds.get(13), rotbonds.get(14),
																rotbonds.get(15), l * (3.14159265358979323846 / 180));
														mol.SetTorsion(rotbonds.get(16), rotbonds.get(17), rotbonds.get(18),
																rotbonds.get(19), m * (3.14159265358979323846 / 180));
														mol.SetTorsion(rotbonds.get(20), rotbonds.get(21), rotbonds.get(22),
																rotbonds.get(23), n * (3.14159265358979323846 / 180));
														mol.SetTorsion(rotbonds.get(24), rotbonds.get(25), rotbonds.get(26),
																rotbonds.get(27), o * (3.14159265358979323846 / 180));
														mol.SetTorsion(rotbonds.get(28), rotbonds.get(29), rotbonds.get(30),
																rotbonds.get(31), p * (3.14159265358979323846 / 180));
														mol.SetTorsion(rotbonds.get(32), rotbonds.get(33), rotbonds.get(34),
																rotbonds.get(35), q * (3.14159265358979323846 / 180));
														mol.SetTorsion(rotbonds.get(36), rotbonds.get(37),
																rotbonds.get(38), rotbonds.get(39),
																r * (3.14159265358979323846 / 180));
														align.SetRefMol(mol_ref);
														align.SetTargetMol(mol);
														align.Align();
														force.Setup(mol);
														double energy = force.Energy();
														double rmsd = align.GetRMSD();
														if (type.equals("pdb")) {
															bufferWriter.write(counter + "," + rmsd
																	+ "," + energy + "\n");
															append = conv.Write(mol);
														} else {
															bufferWriter.write(counter + "," + rmsd + "\n");
															append = conv.Write(mol);
														}
														if(rmsd<=1.0){
															count_rmsd++;
														}
														counter++;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
				bufferfre.write(filename+"\t"+rotnum+"\t"+(counter-1)+"\t"+count_rmsd+"\n");
			}
			
			bufferWriter.close();
			

			long endTime = System.nanoTime() - startTime;
			double seconds = (double) endTime / 1000000000.0;
			DecimalFormat d = new DecimalFormat(".###");
			System.out
					.println("Conformers generated for file " + filename + " in " + d.format(seconds) + " seconds\n\n");
		}
System.out.println("SysGen Completed! Conformers generated for all files ");bufferfre.close();
	}

}