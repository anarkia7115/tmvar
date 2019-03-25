package tmVarlib;
//
// tmVar - Java version
//

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import javax.xml.stream.XMLStreamException;

public class CorrespondGene
{
	public static void main(String [] args) throws IOException, InterruptedException, XMLStreamException, SQLException 
	{
		if(args.length<1)
		{
			System.out.println("\n$ java -Xmx5G -Xms5G -jar CorrespondGene.jar [InputFile] [OutputFile]");
		}
		else
		{
			String InputFile=args [0];
			String OutputFile=args [1];
			
			//String InputFile="../../Project_COSMIC/output/COSMIC_PMIDs_1week.full.xml.BioC.XML";
			//String OutputFile="../../Project_COSMIC/output/COSMIC_PMIDs_1week.full.withCorrespondingGene.xml";
			
			double startTime,endTime,totTime;
			startTime = System.currentTimeMillis();//start time
			
			File f = new File(OutputFile);
			if(f.exists() && !f.isDirectory()) 
			{ 
				System.out.println(OutputFile+" - Done. (The output file exists)");
			}
			else
			{
				BioCConverter BC= new BioCConverter();
				
				/*
				 * Format Check 
				 */
				String Format = "";
				String checkR=BC.BioCFormatCheck(InputFile);
				if(checkR.equals("BioC"))
				{
					Format = "BioC";
				}
				else if(checkR.equals("PubTator"))
				{
					Format = "PubTator";
				}
				else
				{
					System.out.println("Input format is BioC-XML only.");
					System.exit(0);
				}
				
				System.out.print(InputFile+" - ("+Format+" format) : Processing ... \r");
				 
				BC.BioCReaderWithAnnotation(InputFile);
				for (int i = 0; i < BC.PMIDs.size(); i++)
				{
					String Pmid = BC.PMIDs.get(i);
					String FullText=" "; // not cover ref
					HashMap<String,HashMap<String,String>> GeneAnno=new HashMap<String,HashMap<String,String>>();
					HashMap<String,HashMap<String,String>> VariationAnno=new HashMap<String,HashMap<String,String>>();
					HashMap<Integer,Integer> TablePassage=new HashMap<Integer,Integer>();
					HashMap<String,Integer> Gene2Num=new HashMap<String,Integer>(); 
					/** Paragraphs : j */
					for (int j = 0; j < BC.PassageNames.get(i).size(); j++)
					{
						String PassageName= BC.PassageNames.get(i).get(j); // Passage name
						int PassageOffset = BC.PassageOffsets.get(i).get(j); // Passage offset
						String PassageContext = BC.PassageContexts.get(i).get(j); // Passage context
						ArrayList<String> Annotation = BC.Annotations.get(i).get(j); // Annotation
						if(PassageName.toLowerCase().equals("table"))
						{
							TablePassage.put(j,PassageOffset);
						}
						if(!PassageName.toLowerCase().equals("ref"))
						{
							FullText=FullText+PassageContext+" ";
							for(int a=0;a<Annotation.size();a++)
							{
								String anno[]=Annotation.get(a).split("\t");
								String start=anno[0];
								String last=anno[1];
								String type=anno[3];
								String id=anno[4];
								
								String mention = PassageContext.substring(Integer.parseInt(start), Integer.parseInt(last));
								if(type.equals("Gene"))
								{
									if(!GeneAnno.containsKey(id))
									{
										GeneAnno.put(id, new HashMap<String,String>());
									}
									GeneAnno.get(id).put(j+"\t"+start+"\t"+last+"\t"+mention, "");
									if(Gene2Num.containsKey(id))
									{
										Gene2Num.put(id, Gene2Num.get(id)+1);
									}
									else
									{
										Gene2Num.put(id, 1);
									}
								}
								else if(type.matches("(ProteinMutation|DNAMutation)"))
								{
									if(!id.matches(".*RS#:.*"))
									{
										if(!VariationAnno.containsKey(id))
										{
											VariationAnno.put(id, new HashMap<String,String>());
										}
										VariationAnno.get(id).put(j+"\t"+start+"\t"+last+"\t"+mention, "");
									}
								}
							}
						}
						else if(PassageName.toLowerCase().equals("ref")) // only the mentions also in the full text will be included
						{
							for(int a=0;a<Annotation.size();a++)
							{
								String anno[]=Annotation.get(a).split("\t");
								String start=anno[0];
								String last=anno[1];
								String type=anno[3];
								String id=anno[4];
								
								//System.out.println(Pmid+"\t"+PassageContext.length()+"\t"+start+"\t"+last+"\t"+Annotation.get(a));
								
								String mention = PassageContext.substring(Integer.parseInt(start), Integer.parseInt(last));
								if(type.equals("Gene"))
								{
									if(!GeneAnno.containsKey(id))
									{
										GeneAnno.put(id, new HashMap<String,String>());
									}
									GeneAnno.get(id).put(j+"\t"+start+"\t"+last+"\t"+mention, "");
								}
								else if(type.matches("(ProteinMutation|DNAMutation)"))
								{
									if(FullText.toLowerCase().matches(".*"+mention.toLowerCase()+".*"))
									{
										if(!id.matches(".*RS#:.*"))
										{
											if(!VariationAnno.containsKey(id))
											{
												VariationAnno.put(id, new HashMap<String,String>());
											}
											VariationAnno.get(id).put(j+"\t"+start+"\t"+last+"\t"+mention, "");
										}
									}
								}
							}
						}
					}
					
					// find major Gene
					String MajorGene="";
					int MajorGene_Num=0;
					for(String Gene : Gene2Num.keySet())
					{
						if(Gene2Num.get(Gene)>MajorGene_Num)
						{
							MajorGene=Gene;
							MajorGene_Num=Gene2Num.get(Gene);
						}
					}

					// find closed gene
					for(String Varid : VariationAnno.keySet())
					{
						int ClosedDistance=10000;
						String ClosedGene="";
						for(String Geneid : GeneAnno.keySet())
						{
							for(String Varlocation : VariationAnno.get(Varid).keySet())
							{
								String Var_anno[]=Varlocation.split("\t");
								int Var_paragraph_id=Integer.parseInt(Var_anno[0]);
								int Var_start=Integer.parseInt(Var_anno[1]);
								int Var_last=Integer.parseInt(Var_anno[2]);
								String Var_mention=Var_anno[3];
								
								int SurroundText_start=0;
								int SurroundText_last=Var_last+100;
								if(Var_start>100){SurroundText_start=Var_start-100;}
								if((Var_last+100)>BC.PassageContexts.get(i).get(Var_paragraph_id).length()){SurroundText_last=BC.PassageContexts.get(i).get(Var_paragraph_id).length();}
								String SurroundText=BC.PassageContexts.get(i).get(Var_paragraph_id).substring(SurroundText_start, SurroundText_last);
								
								Pattern ptmp = Pattern.compile(".*\\(([^\\)]*)"+Var_mention+"([^\\(]*)\\).*"); //if ....(...Var...).... ; re-assign offset for variations
								Matcher mtmp = ptmp.matcher(SurroundText);
								if(mtmp.find())
								{
									String pre=mtmp.group(1);
									for(String Genelocation : GeneAnno.get(Geneid).keySet())
									{
										String Gene_anno[]=Genelocation.split("\t");
										String Gene_mention=Gene_anno[3];
										Gene_mention=Gene_mention.replaceAll("([\\W\\-\\_])", ".");
										if(!pre.matches(".*"+Gene_mention+".*"))
										{
											Var_start=Var_start-pre.length();
											Var_last=Var_start;
											break;
										}
									}
								}
								
								for(String Genelocation : GeneAnno.get(Geneid).keySet())
								{
									String Gene_anno[]=Genelocation.split("\t");
									int Gene_paragraph_id=Integer.parseInt(Gene_anno[0]);
									int Gene_start=Integer.parseInt(Gene_anno[1]);
									int Gene_last=Integer.parseInt(Gene_anno[2]);
									String Gene_mention=Gene_anno[3];
									if(Var_paragraph_id == Gene_paragraph_id) // paragraph the same
									{
										if(Var_start>Gene_last) // ... Gene ... Var ...
										{
											if(ClosedDistance>(Var_start-Gene_last))
											{
												ClosedDistance=Var_start-Gene_last;
												ClosedGene=Geneid;
											}
										}
										else if(Gene_start>Var_last && (!TablePassage.containsKey(Var_paragraph_id))) // ... Var ... Gene ...
										{
											String betweenString=BC.PassageContexts.get(i).get(Var_paragraph_id).substring(Var_last, Gene_start);
											Gene_mention=Gene_mention.replaceAll("[\\(\\)]", "\\.");
											if(!betweenString.matches(".*[a-z]\\. .*"+Gene_mention+".*")) // Gene is in the next sentence
											{
												if(ClosedDistance>(Gene_start-Var_last))
												{
													ClosedDistance=Gene_start-Var_last;
													ClosedGene=Geneid;
												}
											}
										}
										else if(!TablePassage.containsKey(Var_paragraph_id))// overlap
										{
											if(ClosedDistance>0)
											{
												ClosedDistance=0;
												ClosedGene=Geneid;
											}
										}
									}
								}
							}
						}
						if(ClosedGene.equals(""))
						{
							ClosedGene=MajorGene;
						}
						
						for (int j = 0; j < BC.PassageNames.get(i).size(); j++)
						{
							ArrayList<String> Annotation = BC.Annotations.get(i).get(j); // Annotation
							for(int a=0;a<Annotation.size();a++)
							{
								String anno[]=Annotation.get(a).split("\t");
								String id=anno[4];
								if(Varid.equals(id))
								{
									if(!ClosedGene.equals("")) // find closed gene
									{
										String SequenceCheck="";
										//Check sequence allele
										{
											String type="DNA";
											if(anno[3].equals("ProteinMutation"))
											{
												type="Protein";
											}
											Connection c = null;
											Statement stmt = null;
											try {
												Class.forName("org.sqlite.JDBC");
											} 
											catch ( Exception e ) 
											{
												System.err.println( e.getClass().getName() + ": " + e.getMessage() );
												System.exit(0);
											}
											c = DriverManager.getConnection("jdbc:sqlite:Database/refseq.db");
											stmt = c.createStatement();
											String SQL="SELECT seqid,Seq FROM refseq WHERE geneid='"+ClosedGene+"' and type='"+type+"' order by seqid desc";
											ResultSet data = stmt.executeQuery(SQL);
											while ( data.next() ) 
											{
												String component[]=Varid.split("\\|",-1);
												if(component[1].equals("SUB"))
												{
													String W=component[2];
													String P=component[3];
													String M=component[4];
													P=P.replaceAll("^(IVS|EX)[0-9IV]+[\\+\\-]*", "");
													if(W.matches("[A-Z]") && M.matches("[A-Z]") && P.matches("[0-9]+") && data.getString("Seq").length()>Integer.parseInt(P))
													{
														if(data.getString("Seq").substring(Integer.parseInt(P),Integer.parseInt(P)+1).equals(W))
														{
															SequenceCheck=data.getString("seqid")+"|Matched-W";
															break;
														}
														else if(data.getString("Seq").substring(Integer.parseInt(P),Integer.parseInt(P)+1).equals(M))
														{
															SequenceCheck=data.getString("seqid")+"|Matched-M";
															break;
														}
														if(data.getString("Seq").substring(Integer.parseInt(P)-1,Integer.parseInt(P)).equals(W))
														{
															SequenceCheck=data.getString("seqid")+"|PlusMatched-W";
															break;
														}
														else if(data.getString("Seq").substring(Integer.parseInt(P)-1,Integer.parseInt(P)).equals(M))
														{
															SequenceCheck=data.getString("seqid")+"|PlusMatched-M";
															break;
														}
													}
												}
												if(component[1].matches("(DEL|INS)"))
												{
													String P=component[2];
													String M=component[3];
													P=P.replaceAll("^(IVS|EX)[0-9IV]+[\\+\\-]*", "");
													if(M.matches("[A-Z]") && P.matches("[0-9]+") && data.getString("Seq").length()>Integer.parseInt(P))
													{
														if(data.getString("Seq").substring(Integer.parseInt(P),Integer.parseInt(P)+1).equals(M))
														{
															SequenceCheck=data.getString("seqid")+"|Matched-M";
															break;
														}
														else if(data.getString("Seq").substring(Integer.parseInt(P)-1,Integer.parseInt(P)).equals(M))
														{
															SequenceCheck=data.getString("seqid")+"|PlusMatched-M";
															break;
														}
													}
												}
												else
												{
													//Give Up
												}
											}
											stmt.close();
											c.close();
										}
										if(SequenceCheck.equals(""))
										{
											BC.Annotations.get(i).get(j).set(a, Annotation.get(a)+";CorrespondingGene:"+ClosedGene);
										}
										else
										{
											BC.Annotations.get(i).get(j).set(a, Annotation.get(a)+";CorrespondingGene:"+ClosedGene+";"+SequenceCheck);
										}
									}
								}
							}
						}
					}
				}
				BC.BioCOutput(InputFile,OutputFile);
				
				/*
				 * Time stamp - last
				 */
				endTime = System.currentTimeMillis();//ending time
				totTime = endTime - startTime;
				System.out.println(InputFile+" - ("+Format+" format) : Processing Time:"+totTime/1000+"sec");
			}
		}
	}
}
