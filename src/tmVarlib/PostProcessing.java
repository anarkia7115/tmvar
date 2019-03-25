package tmVarlib;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.sql.*;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.text.BreakIterator;

public class PostProcessing 
{

	public void toME(String Filename,String FilenameOutput,String FilenameLoca,String FilenameME)throws IOException 
	{
		HashMap<String,String> nametothree = new HashMap<String,String>();
		nametothree.put("THYMINE", "T");nametothree.put("ALANINE", "ALA");nametothree.put("ARGININE", "ARG");nametothree.put("ASPARAGINE", "ASN");nametothree.put("ASPARTICACID", "ASP");nametothree.put("ASPARTATE", "ASP");nametothree.put("CYSTEINE", "CYS");nametothree.put("GLUTAMINE", "GLN");nametothree.put("GLUTAMICACID", "GLU");nametothree.put("GLUTAMATE", "GLU");nametothree.put("GLYCINE", "GLY");nametothree.put("HISTIDINE", "HIS");nametothree.put("ISOLEUCINE", "ILE");nametothree.put("LEUCINE", "LEU");nametothree.put("LYSINE", "LYS");nametothree.put("METHIONINE", "MET");nametothree.put("PHENYLALANINE", "PHE");nametothree.put("PROLINE", "PRO");nametothree.put("SERINE", "SER");nametothree.put("THREONINE", "THR");nametothree.put("TRYPTOPHAN", "TRP");nametothree.put("TYROSINE", "TYR");nametothree.put("VALINE", "VAL");nametothree.put("STOP", "XAA");
		HashMap<String,String> threetone = new HashMap<String,String>();
		threetone.put("ALA", "A");threetone.put("ARG", "R");threetone.put("ASN", "N");threetone.put("ASP", "D");threetone.put("CYS", "C");threetone.put("GLN", "Q");threetone.put("GLU", "E");threetone.put("GLY", "G");threetone.put("HIS", "H");threetone.put("ILE", "I");threetone.put("LEU", "L");threetone.put("LYS", "K");threetone.put("MET", "M");threetone.put("PHE", "F");threetone.put("PRO", "P");threetone.put("SER", "S");threetone.put("THR", "T");threetone.put("TRP", "W");threetone.put("TYR", "Y");threetone.put("VAL", "V");threetone.put("ASX", "B");threetone.put("GLX", "Z");threetone.put("XAA", "X");threetone.put("TER", "X");
		HashMap<String,String> sentence_hash = new HashMap<String,String>();
		HashMap<String,String> article_hash = new HashMap<String,String>();
		HashMap<String,String> TypeOrder_hash = new HashMap<String,String>();
		ArrayList<String> pmidARR = new ArrayList<String>(); 
		
		try {
			/*
			 * load input sentences
			 */
			int ParagraphType_count=0;
			BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(Filename), "UTF-8"));
			String line;
			while ((line = inputfile.readLine()) != null)  
			{
				if(line.contains("|")) //Title|Abstract
	        	{
					String Pmid="";
					String ParagraphType="";
					String ParagraphContent="";
					Pattern pat = Pattern.compile("^([^\\|\\t]+)\\|([^\\|\\t]+)\\|(.*)$");
					Matcher mat = pat.matcher(line);
					if(mat.find()) //Title|Abstract
		        	{
						Pmid = mat.group(1);
						ParagraphType=mat.group(2);
						ParagraphContent=mat.group(3);
						if(ParagraphContent.equals(""))
						{
							ParagraphContent="- No text -";
						}
					}
					sentence_hash.put(Pmid+"\t"+ParagraphType+"\t"+ParagraphType_count,ParagraphContent);
					if(article_hash.get(Pmid) != null)
					{
						article_hash.put(Pmid,article_hash.get(Pmid)+" "+ParagraphContent);
					}
					else
					{
						article_hash.put(Pmid,ParagraphContent);
					}
					if(TypeOrder_hash.get(Pmid) != null)
					{
						TypeOrder_hash.put(Pmid,TypeOrder_hash.get(Pmid)+"\t"+ParagraphType+"|"+ParagraphType_count);
					}
					else
					{
						TypeOrder_hash.put(Pmid,ParagraphType+"|"+ParagraphType_count);
					}
					if(!pmidARR.contains(Pmid))
					{
						pmidARR.add(Pmid);
					}
					ParagraphType_count++;
	        	}
			}
			inputfile.close();
			
			/*
			 * load CRF output
			 */
			ArrayList<String> outputArr = new ArrayList<String>(); 
			inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenameOutput), "UTF-8"));
			while ((line = inputfile.readLine()) != null)  
			{
				outputArr.add(line);
			}
			inputfile.close();
			
			/*
			 * load location
			 */
			ArrayList<String> locationArr = new ArrayList<String>(); 
			inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenameLoca), "UTF-8"));
			while ((line = inputfile.readLine()) != null)  
			{
				locationArr.add(line);
			}
			inputfile.close();
			
			BufferedWriter FileME = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(FilenameME), "UTF-8"));
			
			String pmid="";
			for(int i=0;i<outputArr.size();i++)
			{
				String outputs[]=outputArr.get(i).split("\\t");
				/*
				 * Extract the states and tokens from CRF++ output file
				 */
				int start=100000;
				int last=0;
				String mention="";
				String location[]=locationArr.get(i).split("\\t");
				HashMap<String,String> component_hash = new HashMap<String,String>();
				if(component_hash.get("A") == null){component_hash.put("A", "");}
				if(component_hash.get("T") == null){component_hash.put("T", "");}
				if(component_hash.get("P") == null){component_hash.put("P", "");}
				if(component_hash.get("W") == null){component_hash.put("W", "");}
				if(component_hash.get("M") == null){component_hash.put("M", "");}
				if(component_hash.get("F") == null){component_hash.put("F", "");}
				if(component_hash.get("S") == null){component_hash.put("S", "");}
				if(component_hash.get("D") == null){component_hash.put("D", "");}
				if(component_hash.get("I") == null){component_hash.put("I", "");}
				if(component_hash.get("R") == null){component_hash.put("R", "");}
				
				// print title/abstracts
				if(location.length > 1 && !pmid.equals(location[0]))
				{
					if(!pmid.equals("")){FileME.write("\n");}
					pmid=location[0];
					String TypeOrder[]=TypeOrder_hash.get(pmid).split("\\t");
					for(int j=0;j<TypeOrder.length;j++)
					{
						String[] TypeOrder_component=TypeOrder[j].split("\\|");
						FileME.write(pmid+"|"+TypeOrder_component[0]+"|"+sentence_hash.get(pmid+"\t"+TypeOrder_component[0]+"\t"+TypeOrder_component[1])+"\n");
					}
				}
				
				//find mention
				Pattern pat = Pattern.compile("([ATPWMFSDIR])$");
				Matcher mat = pat.matcher(outputs[outputs.length-1]);
				if((!outputs[0].equals(";")) && mat.find())
				{
					String prestate="";
					mat = pat.matcher(outputs[outputs.length-1]);
					
					//until the end token of the mention
					while((!outputs[0].equals(";")) && mat.find())
					{
						String state=mat.group(1);	
						String locationWhile[]=locationArr.get(i).split("\\t");
						String tkn=locationWhile[1];
						pmid=locationWhile[0];
						int start_tmp=Integer.parseInt(locationWhile[2]);
						int last_tmp=Integer.parseInt(locationWhile[3]);
						mention = mention + tkn;
						if(!component_hash.get(state).equals("") && !state.equals(prestate))
						{
							component_hash.put(state, component_hash.get(state)+","+tkn);
						}
						else
						{
							component_hash.put(state, component_hash.get(state)+tkn);
						}
						if(start_tmp<start){start=start_tmp;}
						if(last_tmp>last){last=last_tmp;}
						prestate=state;
						i++;
						outputs=outputArr.get(i).split("\\t");
						mat = pat.matcher(outputs[outputs.length-1]);
					}
					
					/*
					 * Recognize the components(identifiers)
					 */
					String identifier;
					String type="";
					if(!component_hash.get("D").equals(""))
					{
						identifier=component_hash.get("A")+"|"+component_hash.get("T")+"|"+component_hash.get("P")+"|"+component_hash.get("M")+"|"+component_hash.get("D");
					}
					else if(!component_hash.get("T").equals(""))
					{
						identifier=component_hash.get("A")+"|"+component_hash.get("T")+"|"+component_hash.get("P")+"|"+component_hash.get("M");
					}
					else if(!component_hash.get("S").equals(""))
					{
						identifier=component_hash.get("A")+"|"+component_hash.get("W")+"|"+component_hash.get("P")+"|"+component_hash.get("M")+"|"+component_hash.get("F")+"|"+component_hash.get("S");
						type="ProteinMutation";
					}
					else if(!component_hash.get("F").equals(""))
					{
						identifier=component_hash.get("A")+"|"+component_hash.get("W")+"|"+component_hash.get("P")+"|"+component_hash.get("M")+"|"+component_hash.get("F");
						type="ProteinMutation";
					}
					else if(!component_hash.get("R").equals(""))
					{
						identifier=mention;
						type="SNP";
					}
					else
					{
						identifier=component_hash.get("A")+"|"+component_hash.get("W")+"|"+component_hash.get("P")+"|"+component_hash.get("M");
					}
					
					/*
					 * Recognize the type of mentions: ProteinMutation | DNAMutation | SNP 
					 */
					if(type.equals(""))
					{
						if(!component_hash.get("A").equals("") &&
						   (
						   component_hash.get("A").toLowerCase().equals("c") ||
						   component_hash.get("A").toLowerCase().equals("r") ||
						   component_hash.get("A").toLowerCase().equals("m") ||
						   component_hash.get("A").toLowerCase().equals("g")
						   )
						  )
						{
							type="DNAMutation";
						}
						else if(!component_hash.get("A").equals("") && component_hash.get("A").toLowerCase().equals("p"))
						{
							type="ProteinMutation";
						}
						else if(!component_hash.get("T").equals("") && component_hash.get("T").toLowerCase().equals("delta"))
						{
							type="DNAMutation";
						}
						else if(
								!component_hash.get("P").equals("") && 
								(
								component_hash.get("P").toLowerCase().equals("ex") ||
								component_hash.get("P").toLowerCase().equals("intron") ||
								component_hash.get("P").toLowerCase().equals("ivs")
								)
							   )
						{
							type="DNAMutation";
						}
						else if(!component_hash.get("M").equals("") && !component_hash.get("W").equals("")) //others
						{
							pat = Pattern.compile("^[ATCGatcgu]+$");
							Matcher mat_M = pat.matcher(component_hash.get("M"));
							Matcher mat_W = pat.matcher(component_hash.get("W"));
							
							if(!mat_M.find() || !mat_W.find())
							{
								type="ProteinMutation";
							}
							else //default
							{
								type="DNAMutation";
							}
						}
					}
					
					/*
					 * filtering and Print out
					 */
					if( (component_hash.get("W").length() == 3 || component_hash.get("M").length() == 3) 
							&& component_hash.get("W").length() != component_hash.get("M").length()
							&& !component_hash.get("W").equals("") && !component_hash.get("M").equals("") && component_hash.get("W").indexOf(",")!=-1 && component_hash.get("M").indexOf(",")!=-1
							&& ((component_hash.get("W").indexOf("A")!=-1 && component_hash.get("W").indexOf("T")!=-1 && component_hash.get("W").indexOf("C")!=-1 && component_hash.get("W").indexOf("G")!=-1) || (component_hash.get("M").indexOf("A")!=-1 && component_hash.get("M").indexOf("T")!=-1 && component_hash.get("M").indexOf("C")!=-1 && component_hash.get("M").indexOf("G")!=-1))
							&& component_hash.get("T").equals("")
						)
							{/*System.out.println("filtering 1:"+article_hash.get(pmid).substring(start-1,last));*/}
					else if((component_hash.get("M").matches("[ISQMNPKDFHLRWVEYX]") || component_hash.get("W").matches("[ISQMNPKDFHLRWVEYX]")) && component_hash.get("P").matches("[6-9][0-9][0-9][0-9]+")){/*System.out.println("filtering 2:"+article_hash.get(pmid).substring(start-1,last));*/} //length > 3000, if protein mutation
					else if(component_hash.get("M").equals("") && component_hash.get("T").equals("") && component_hash.get("F").equals("") && component_hash.get("P").equals("") && !type.equals("SNP")){/*System.out.println("filtering 2.2:"+article_hash.get(pmid).substring(start-1,last));*/} //Arg235
					else if(component_hash.get("W").equals("") && component_hash.get("T").equals("") && component_hash.get("F").equals("") && component_hash.get("D").equals("") && (!component_hash.get("M").equals("")) && !type.equals("SNP")){/*System.out.println("filtering 2.2:"+article_hash.get(pmid).substring(start-1,last));*/} //314C
					else if(component_hash.get("M").equals("") && component_hash.get("T").equals("") && component_hash.get("F").equals("")&& component_hash.get("D").equals("") && (!component_hash.get("W").equals("")) && !type.equals("SNP")){/*System.out.println("filtering 2.2:"+article_hash.get(pmid).substring(start-1,last));*/} //C314
					else if(component_hash.get("T").toLowerCase().equals("delta") && (!component_hash.get("M").matches("[ATCG0-9]*"))  && component_hash.get("P").equals("")){/*System.out.println("filtering 3:"+article_hash.get(pmid).substring(start-1,last));*/} //DeltaEVIL
					else if(component_hash.get("T").toLowerCase().equals("delta") && component_hash.get("M").equals("")  && component_hash.get("P").equals("")){/*System.out.println("filtering 3:"+article_hash.get(pmid).substring(start-1,last));*/} //Delta
					else if(component_hash.get("P").matches("^-") && type.equals("ProteinMutation")){/*System.out.println("filtering 4:"+article_hash.get(pmid).substring(start-1,last));*/} //negative protein mutation
					else if(component_hash.get("W").matches("^[BJOUZ]") || component_hash.get("M").matches("^[BJOUZ]")){/*System.out.println("filtering 5:"+article_hash.get(pmid).substring(start-1,last));*/} //not a mutation
					else if(component_hash.get("W").matches("^[A-Za-z][a-z]") || component_hash.get("M").matches("^[A-Za-z][a-z]") ){/*System.out.println("filtering 7:"+article_hash.get(pmid).substring(start-1,last));*/} //not a mutation
					else if( component_hash.get("W").matches("[A-Za-z]") && component_hash.get("M").matches("[A-Za-z][a-z][a-z]+") && (!nametothree.containsKey(component_hash.get("M").toUpperCase())) && (!threetone.containsKey(component_hash.get("M").toUpperCase())) ){} //T-->infinity
					else if( component_hash.get("M").matches("[A-Za-z]") && component_hash.get("W").matches("[A-Za-z][a-z][a-z]+") && (!nametothree.containsKey(component_hash.get("W").toUpperCase())) && (!threetone.containsKey(component_hash.get("W").toUpperCase())) ){} //T-->infinity
					else if(article_hash.get(pmid).length()>start-1 && start>15 && article_hash.get(pmid).substring(start-15,start-1).matches(".*(Figure|Table|Fig|Figs|Tab|Tabs|figure|table|fig|figs|tab|tabs)[ \\.].*")){/*System.out.println("filtering 6.1:"+start+"\t"+last+"\t"+article_hash.get(pmid).substring(start-1,last));*/} //Figure S13A
					else if(article_hash.get(pmid).length()>start-1 && start>2 && article_hash.get(pmid).substring(start-2,start-1).matches("[0-9A-Za-z]") && (!mention.matches("-[0-9]+.*"))){/*System.out.println("filtering 6.4:"+article_hash.get(pmid).substring(start-1,last));*/} //not empty
					else if(article_hash.get(pmid).length()>last+1 && article_hash.get(pmid).substring(last,last+1).matches("[A-Za-z0-9\\>]") && !type.equals("SNP")){/*System.out.println("filtering 7:"+article_hash.get(pmid).substring(start-1,last));*/} //V79 Chinese
					else if(mention.matches(".+\\).+\\(.+")){/*System.out.println("filtering - Mention:"+mention);*/} //Arg)3(Ser
					else if(mention.matches(".+\\>.+\\>.+")){/*System.out.println("filtering - Mention:"+mention);*/} //G>T>C
					else if(mention.matches("[A-Z][a-z][a-z][0-9]+,[A-Z][a-z][a-z][0-9]+")){/*System.out.println("filtering - Mention:"+mention);*/} //Glu207, Thr209
					else if(mention.matches("[A-Z][a-z][a-z][0-9]+,[A-Z][a-z][a-z][0-9]+,[A-Z][a-z][a-z][0-9]+")){/*System.out.println("filtering - Mention:"+mention);*/} //Glu207, Thr209, Gly211
					else if(mention.matches(".*[gG]\\/L.*")){/*System.out.println("filtering - Mention:"+mention);*/} //450 G/L
					else if(tmVar.PAM_lowerScorePair.contains(component_hash.get("M")+"\t"+component_hash.get("W"))){/*System.out.println("filtering - Mention:"+mention);*/} //unlikely to occur
					else if(article_hash.get(pmid).length()>last)
					{
						FileME.write(pmid+"\t"+(start-1)+"\t"+last+"	"+article_hash.get(pmid).substring(start-1,last)+"\t"+type+"\t"+identifier+"\n");
					}
					
					if(!outputs[outputs.length-1].equals("O"))
					{
						i--;
					}
				}
			}
			FileME.write("\n");
			FileME.close();
		}
		catch(IOException e1){ System.out.println("[toME]: Input file is not exist.");}
	}

	public void toPostME(String FilenameME,String FilenamePostME)throws IOException
	{
		Pattern Pattern_Component_1 = Pattern.compile("^([RrSs][Ss][ ]*[0-9]+)[ ]*(and|/|,|or)[ ]*([RrSs][Ss][ ]*[0-9]+)$");
		Pattern Pattern_Component_2 = Pattern.compile("^(.*[^0-9])[ ]*([RrSs][Ss][ ]*[0-9]+)$");
		Pattern Pattern_Component_3 = Pattern.compile("^([RrSs][Ss][ ]*[0-9]+)[ ]*([^0-9].*)$");
		
		try {
			ArrayList<String> MF_Pattern = new ArrayList<String>();
			ArrayList<String> MF_Type = new ArrayList<String>();
			
			/*
			 * Append-pattern
			 */
			BufferedReader RegEx_NL = new BufferedReader(new InputStreamReader(new FileInputStream("lib/RegEx/MF.RegEx.2.txt"), "UTF-8"));
			String line="";
			while ((line = RegEx_NL.readLine()) != null)  
			{
				String RegEx[]=line.split("\t");
				MF_Pattern.add(RegEx[0]);
				MF_Type.add(RegEx[1]);
			}
			RegEx_NL.close();
					
			BufferedWriter FilePostME = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(FilenamePostME), "UTF-8")); // .location
			BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenameME), "UTF-8"));
			ArrayList<String> Annotation = new ArrayList<String>();
			String article="";
			String Pmid="";
			while ((line = inputfile.readLine()) != null)  
			{
				Pattern pat = Pattern.compile("^([^\\|\\t]+)\\|([^\\|\\t]+)\\|(.*)$");
				Matcher mat = pat.matcher(line);
				if(mat.find()) //Title|Abstract
	        	{
					Pmid = mat.group(1);
					String ParagraphContent=mat.group(3);
					article=article+ParagraphContent+"\t";
					FilePostME.write(line+"\n");
				}
				else if (line.contains("\t")) //Annotation
		    	{
					Annotation.add(line);
		    	}
				else if(line.length()==0) //Processing
				{
					ArrayList<Integer> RemoveAnno = new ArrayList<Integer>();
					
					/*
					 * Split RSnumber
					 */
					for(int i=0;i<Annotation.size();i++)
					{
						String anno[]=Annotation.get(i).split("\t");
						
			    		int start = Integer.parseInt(anno[1]);
	        			int last = Integer.parseInt(anno[2]);
	        			String mention = anno[3];
	        			Matcher m1 = Pattern_Component_1.matcher(mention);
	        			Matcher m2 = Pattern_Component_2.matcher(mention);
	        			Matcher m3 = Pattern_Component_3.matcher(mention);
	        			if(m1.find())
	        			{
	        				String sub_mention1=m1.group(1);
	        				String sub_mention2=m1.group(3);
	        				int start1=start;
	        				int last1=start+sub_mention1.length();
	        				int start2=last-sub_mention2.length();
	        				int last2=last;
	        				RemoveAnno.add(i);
	        				Annotation.add(Pmid+"\t"+start1+"\t"+last1+"\t"+sub_mention1+"\tSNP\t"+sub_mention1);
	        				Annotation.add(Pmid+"\t"+start2+"\t"+last2+"\t"+sub_mention2+"\tSNP\t"+sub_mention2);
	        			}
	        			else if(m2.find())
	        			{
	        				String sub_mention1=m2.group(1);
	        				String sub_mention2=m2.group(2);
	        				int start1=start;
	        				int last1=start+sub_mention1.length();
	        				int start2=last-sub_mention2.length();
	        				int last2=last;
	        				RemoveAnno.add(i);
	        				Annotation.add(Pmid+"\t"+start1+"\t"+last1+"\t"+sub_mention1+"\tDNAMutation\t"+sub_mention1);
	        				Annotation.add(Pmid+"\t"+start2+"\t"+last2+"\t"+sub_mention2+"\tSNP\t"+sub_mention2);
	        			}
	        			else if(m3.find())
	        			{
	        				String sub_mention1=m3.group(1);
	        				String sub_mention2=m3.group(2);
	        				int start1=start;
	        				int last1=start+sub_mention1.length();
	        				int start2=last-sub_mention2.length();
	        				int last2=last;
	        				RemoveAnno.add(i);
	        				Annotation.add(Pmid+"\t"+start1+"\t"+last1+"\t"+sub_mention1+"\tSNP\t"+sub_mention1);
	        				Annotation.add(Pmid+"\t"+start2+"\t"+last2+"\t"+sub_mention2+"\tDNAMutation\t"+sub_mention2);
	        			}
	        		}
					for(int i=RemoveAnno.size()-1;i>=0;i--)
					{
						int RemA=RemoveAnno.get(i);
						Annotation.remove(RemA);
					}
					
					/*
					 * Boundary
					 */
					for(int i=0;i<Annotation.size();i++)
					{
						String anno[]=Annotation.get(i).split("\t");
						int start = Integer.parseInt(anno[1]);
	        			int last = Integer.parseInt(anno[2]);
	        			String mention = anno[3];
	        			String type = anno[4];
	        			String identifier = anno[5];
	        			int check=0;
	        			
	        			if(mention.matches("^[0-9]") && article.substring(start-1,start).matches("[+-]"))	//17000021	251	258	1858C>T --> +1858C>T
	        			{
	        				check=1;
	        				start=start-1;
	        				mention=article.substring(start-1,start)+mention;
	        			}
	        			if(mention.matches("^[^0-9A-Za-z][RrSs][Ss][0-9]+$"))	//_rs7207916
						{
	        				check=1;
	        				start=start+1;
	        				mention=mention.substring(1,mention.length());
	        			}
	        			if(mention.startsWith("(") && !mention.contains(")")) // delete (
						{
	        				check=1;
	        				start=start+1;
	        				mention=mention.substring(1,mention.length());
						}
	        			
	        			//mention.substring(mention.length()-1, mention.length()) : the last character of mention
	        			
	        			else if(mention.substring(mention.length()-1, mention.length()).equals(")") && !mention.contains("(")) // delete )
	        			{
	        				check=1;
	        				last=last-1;
	        				mention=mention.substring(0,mention.length()-1);
	        			}
	        			else if(start > 0 && article.substring(start-1,start).equals("(") && mention.contains(")") && !mention.contains("(")) // add (
	        			{
	        				check=1;
	        				start=start-1;
	        				mention="("+mention;
	        			}
	        			else if(start > 0 && article.substring(start-1,start).equals("[") && mention.contains("]") && !mention.contains("[")) // add [
	        			{
	        				check=1;
	        				start=start-1;
	        				mention="["+mention;
	        			}
	        			else if(article.substring(last,last+1).equals(")") && mention.contains("(") && !mention.contains(")")) // add )
						{
	        				check=1;
	        				last=last+1;
	        				mention=mention+")";
						}
	        			else if(article.substring(last,last+1).equals("]") && mention.contains("[") && !mention.contains("]")) // add ]
						{
	        				check=1;
	        				last=last+1;
	        				mention=mention+"]";
						}
	        			else if(mention.startsWith("(") && mention.substring(mention.length()-1, mention.length()).equals(")")) // delete (  )
	        			{
	        				check=1;
	        				start=start+1;
	        				last=last-1;
	        				mention=mention.substring(1,mention.length()-1);
	        			}
	        			if(check == 1)
	        			{	
	        				Annotation.remove(i);
	        				Annotation.add(i,Pmid+"\t"+start+"\t"+last+"\t"+mention+"\t"+type+"\t"+identifier); //problem not solve!!!
	        			}
					}
					
	    			/*
	    			 *  Mention Recognition by Pattern
	    			 */
					
					/*
					 * Self-pattern: should close in Full text
					 */
					/*
        			HashMap<String,String> SelfPattern2type_hash = new HashMap<String,String>();
					for(int i=0;i<Annotation.size();i++)
					{
						String anno[]=Annotation.get(i).split("\t");
						
			    		String mention = anno[3];
	        			String type = anno[4];
	        			if(!type.equals("SNP"))
	        			{
	        				mention = mention.replaceAll("([^A-Za-z0-9])","\\$1");
	        				mention = mention.replaceAll("[0-9]+","[0-9]+");
	        				mention = mention.replaceAll("(IVS|EX)","@@@@");
	        				mention = mention.replaceAll("(rs|ss)","@@@");
	        				mention = mention.replaceAll("[A-Z]","[A-Z]");
	        				mention = mention.replaceAll("@@@@","(IVS|EX)");
	        				mention = mention.replaceAll("@@@","(rs|ss)");
	        				SelfPattern2type_hash.put(mention, type);
	        			}
					}
	    			*/
	    			
					/*
	    			 * Pattern match
	    			 */
					HashMap <String,String> AddList = new HashMap <String,String>();
					ArrayList <String> removeList = new ArrayList <String>();
					String article_tmp=article;
					for(int pat_i=0;pat_i<MF_Pattern.size();pat_i++)
					{
						String pattern_RegEx = MF_Pattern.get(pat_i); 
						String pattern_Type = MF_Type.get(pat_i);
						
						Pattern PATTERN_RegEx = Pattern.compile("^(.*[^A-Za-z0-9])("+pattern_RegEx+")[^A-Za-z0-9]");
						Matcher mp = PATTERN_RegEx.matcher(article_tmp);
						boolean ExistLarger = false;
						while (mp.find()) 
						{
							String pre = mp.group(1);
							String mention = mp.group(2);
							
							int start=pre.length();
							int last=start+mention.length();
							//Check if overlap with previous annotation
							for(int a=0;a<Annotation.size();a++)
							{
								String Exist[] = Annotation.get(a).split("\t");
								int Exist_start = Integer.parseInt(Exist[1]);
								int Exist_last = Integer.parseInt(Exist[2]);
								if( (start<Exist_start && last>=Exist_last) || (start<=Exist_start && last>Exist_last) )
								{
									removeList.add(Annotation.get(a));
								}
								else if (((start>Exist_start && start<Exist_last) || (last>Exist_start && last<Exist_last)) && ((last-start)>(Exist_last-Exist_start))) //overlap, but not subset
								{
									removeList.add(Annotation.get(a));
								}
								else if( (Exist_start<start && Exist_last>=last) || (Exist_start<=start && Exist_last>last) )
								{
									ExistLarger = true;
								}
							}
							if(ExistLarger == false)
							{	
								//Check if overlap with previous annotation
								for(String added : AddList.keySet())
								{
									String Exist[] = added.split("\t");
									int Exist_start = Integer.parseInt(Exist[1]);
									int Exist_last = Integer.parseInt(Exist[2]);
									if( (start<Exist_start && last>=Exist_last) || (start<=Exist_start && last>Exist_last) )
									{
										AddList.put(added,"remove");
									}
									else if ( ((start>Exist_start && start<Exist_last) || (last>Exist_start && last<Exist_last)) && ((last-start)>(Exist_last-Exist_start)) ) //overlap, but not subset
									{
										AddList.put(added,"remove");
									}
									else if( (Exist_start<start && Exist_last>=last) || (Exist_start<=start && Exist_last>last) )
									{
										ExistLarger = true;
									}
								}
								if(ExistLarger == false && !AddList.containsKey(Pmid+"\t"+start+"\t"+last+"\t"+mention))
								{
									AddList.put(Pmid+"\t"+start+"\t"+last+"\t"+mention,pattern_Type);
								}
							}
							String tmp="";
	        				for(int j=0;j<mention.length();j++){tmp=tmp+"@";}
	        				article_tmp=article_tmp.substring(0,start)+tmp+article_tmp.substring(last);
							mp = PATTERN_RegEx.matcher(article_tmp);
						}
					}
					for(int r=0;r<removeList.size();r++)
					{
						Annotation.remove(removeList.get(r));
					}
					for(String added : AddList.keySet())
					{
						if(!AddList.get(added).equals("remove"))
						{
							boolean found= false;
							for(int a=0;a<Annotation.size();a++)
							{
								String Anno[] = Annotation.get(a).split("\t");
								if(added.equals(Anno[0]+"\t"+Anno[1]+"\t"+Anno[2]+"\t"+Anno[3]))
								{
									found = true;
								}
							}
							if(found == false)
							{
								Annotation.add(added+"\t"+AddList.get(added));
							}
						}
					}
    				
    				/*
	    			 * Pattern match : Self pattern
	    			 */
    				/*
    				for(String pattern_RegEx : SelfPattern2type_hash.keySet())
					{
    					String pattern_Type = SelfPattern2type_hash.get(pattern_RegEx);
				   		Pattern PATTERN_RegEx = Pattern.compile("^(.*[^A-Za-z0-9])("+pattern_RegEx+")[^A-Za-z0-9]");
						Matcher mp = PATTERN_RegEx.matcher(article);
						while (mp.find()) 
						{
							String pre = mp.group(1);
							String mention = mp.group(2);
							int start=pre.length();
							int last=start+mention.length();
							Annotation.add(Pmid+"\t"+start+"\t"+last+"\t"+mention+"\t"+pattern_Type);
							String tmp="";
	        				for(int j=0;j<mention.length();j++){tmp=tmp+"@";}
							article=article.substring(0,start)+tmp+article.substring(last);
							System.out.println(mention);
						}
					}
					*/
	    			
					for(int i=0;i<Annotation.size();i++)
					{
						FilePostME.write(Annotation.get(i)+"\n");
					}
					
					FilePostME.write("\n");
					article="";
					Annotation.clear();
				}
			}
			FilePostME.close();
		}
		catch(IOException e1){ System.out.println("[toPostME]: Input file is not exist.");}
	}
	public void toPostMEData(String Filename,String FilenamePostME,String FilenamePostMeML,String FilenamePostMeData,String TrainTest) throws IOException
	{
		try
		{
			//Parse identifier (components)
			Pattern Pattern_Component_1 = Pattern.compile("^([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)\\|(fs[^|]*)\\|([^|]*)$");
			Pattern Pattern_Component_1_1 = Pattern.compile("^([^|]*)\\|([^|]*(ins|del|Del|dup|-)[^|]*)\\|([^|]*)\\|([^|]*)\\|(fs[^|]*)$"); //append for p.G352fsdelG	p|del|352|G,G|fs
			Pattern Pattern_Component_2 = Pattern.compile("^([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)\\|(fs[^|]*)$");
			Pattern Pattern_Component_3 = Pattern.compile("^([^|]*)\\|([^|]*(ins|del|Del|dup|-)[^|]*)\\|([^|]*)\\|([^|]*)$");
			Pattern Pattern_Component_4 = Pattern.compile("^([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)$");
			Pattern Pattern_Component_5 = Pattern.compile("^([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)$");
			Pattern Pattern_Component_6 = Pattern.compile("^((\\[rs\\]|[RrSs][Ss]|reference SNP no[.] )[0-9][0-9][0-9]+)$");
			
			HashMap<String,String> mention_hash = new HashMap<String,String>();
			HashMap<String,String> mention2type_hash = new HashMap<String,String>();
			BufferedReader inputfile;
			if(TrainTest.equals("Train"))
			{
				inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(Filename), "UTF-8"));
			}
			else
			{
				inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenamePostME), "UTF-8"));
			}
			
			String line;
			while ((line = inputfile.readLine()) != null)  
			{
				Pattern pat = Pattern.compile("^([^\\|\\t]+)\\|([^\\|\\t]+)\\|(.*)$");
				Matcher mat = pat.matcher(line);
				if(mat.find()) //Title|Abstract
	        	{
					
	        	}
				else if (line.contains("\t")) //Annotation
				{
					String anno[]=line.split("\t");
					if(anno.length>=6 && TrainTest.equals("Train"))
	        		{
	        			String mention=anno[3];
	        			String type=anno[4];
	        			String identifier=anno[5];
	        			mention2type_hash.put(mention, type);
	        			mention_hash.put(mention, identifier);
	        		}
					else if(anno.length>=5)
	        		{
	        			String mention=anno[3];
	        			String type=anno[4];
	        			mention2type_hash.put(mention, type);
	        			mention_hash.put(mention, type);
	        		}
				}
			}
			inputfile.close();
			
			BufferedWriter mentionlistbw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(FilenamePostMeML), "UTF-8")); // .ml
			BufferedWriter mentiondata = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(FilenamePostMeData), "UTF-8")); // .location
			for(String mention : mention_hash.keySet() )
			{
				HashMap<Integer, String> character_hash = new HashMap<Integer, String>();
				
				int start=0;
    			int last=mention.length();
    			for(int s=start;s<last;s++)
    			{
    				character_hash.put(s,"I");
    			}
    			
    			if(TrainTest.equals("Train"))
    			{
    				Matcher m1 = Pattern_Component_1.matcher(mention_hash.get(mention));
    				Matcher m1_1 = Pattern_Component_1_1.matcher(mention_hash.get(mention));
        			Matcher m2 = Pattern_Component_2.matcher(mention_hash.get(mention));
        			Matcher m3 = Pattern_Component_3.matcher(mention_hash.get(mention));
        			Matcher m4 = Pattern_Component_4.matcher(mention_hash.get(mention));
        			Matcher m5 = Pattern_Component_5.matcher(mention_hash.get(mention));
        			Matcher m6 = Pattern_Component_6.matcher(mention_hash.get(mention));
        			
        			if(m1.find())
	    			{
        				String type[]=m1.group(1).split(",");
	    				String W[]=m1.group(2).split(",");
	    				String P[]=m1.group(3).split(",");
	    				String M[]=m1.group(4).split(",");
	    				String F[]=m1.group(5).split(",");
	    				String S[]=m1.group(6).split(",");
	    				String mention_tmp=mention;
	    				for(int i=0;i<type.length;i++)
	    				{
	    					String patt="^(.*?)("+type[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    		    				character_hash.put(j,"A");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1\ttype\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<W.length;i++)
	    				{
	    					String patt="^(.*?)("+W[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"W");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1\tW\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<P.length;i++)
	    				{
	    					P[i]=P[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*?)("+P[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						for(int j=mtmp.group(1).length();j<(mtmp.group(1).length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"P");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1\tP\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<M.length;i++)
	    				{
	    					M[i]=M[i].replace("*", "\\*");
	    					String patt="^(.*?)("+M[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						for(int j=mtmp.group(1).length();j<(mtmp.group(1).length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"M");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1\tM\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<F.length;i++)
	    				{
	    					F[i]=F[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*?)("+F[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"F");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1\tF\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<S.length;i++)
	    				{
	    					String patt="^(.*?)("+S[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"S");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1\tS\t"+mention);
	    					}
	    				}
	    			}
        			else if(m1_1.find())
	    			{
        				String type[]=m1_1.group(1).split(",");
	    				String T[]=m1_1.group(2).split(",");
	    				String P[]=m1_1.group(4).split(",");
	    				String M[]=m1_1.group(5).split(",");
	    				String F[]=m1_1.group(6).split(",");
	    				String mention_tmp=mention;
	    				for(int i=0;i<type.length;i++)
	    				{
	    					String patt="^(.*?)("+type[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    		    				character_hash.put(j,"A");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1_1\ttype\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<T.length;i++)
	    				{
	    					String patt="^(.*?)("+T[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"T");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1_1\tT\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<P.length;i++)
	    				{
	    					P[i]=P[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*)("+P[i]+")(.*)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						for(int j=mtmp.group(1).length();j<(mtmp.group(1).length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"P");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1_1\tP\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<M.length;i++)
	    				{
	    					M[i]=M[i].replace("*", "\\*");
	    					String patt="^(.*?)("+M[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						for(int j=mtmp.group(1).length();j<(mtmp.group(1).length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"M");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1_1\tM\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<F.length;i++)
	    				{
	    					F[i]=F[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*?)("+F[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"F");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m1_1\tF\t"+mention);
	    					}
	    				}
	    			}
        			else if(m2.find())
	    			{
	    				String type[]=m2.group(1).split(",");
	    				String W[]=m2.group(2).split(",");
	    				String P[]=m2.group(3).split(",");
	    				String M[]=m2.group(4).split(",");
	    				String F[]=m2.group(5).split(",");
	    				String mention_tmp=mention;
	    				for(int i=0;i<type.length;i++)
	    				{
	    					String patt="^(.*?)("+type[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"A");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m2\tType\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<W.length;i++)
	    				{
	    					String patt="^(.*?)("+W[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"W");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m2\tW\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<P.length;i++)
	    				{
	    					P[i]=P[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*?)("+P[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"P");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m2\tP\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<M.length;i++)
	    				{
	    					M[i]=M[i].replace("*", "\\*");
	    					String patt="^(.*?)("+M[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"M");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m2\tM\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<F.length;i++)
	    				{
	    					F[i]=F[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*?)("+F[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"F");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m2\tF\t"+mention);
	    					}
	    				}
	    			}
	    			else if(m3.find())
	    			{
	    				String type[]=m3.group(1).split(",");
	    				String T[]=m3.group(2).split(",");
	    				String P[]=m3.group(4).split(",");
	    				String M[]=m3.group(5).split(",");
	    				String mention_tmp=mention;
	    				for(int i=0;i<type.length;i++)
	    				{
	    					String patt="^(.*?)("+type[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"A");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m3\tType\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<P.length;i++)
	    				{
	    					P[i]=P[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*?)("+P[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"P");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m3\tP\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<T.length;i++)
	    				{
	    					String patt="^(.*?)("+T[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"T");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m3\tT\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<M.length;i++)
	    				{
	    					M[i]=M[i].replace("*", "\\*");
	    					String patt="^(.*?)("+M[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"M");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m3\tM\t"+mention);
	    					}
	    				}
	    			}
	    			else if(m4.find())
	    			{
	    				String type[]=m4.group(1).split(",");
	    				String W[]=m4.group(2).split(",");
	    				String P[]=m4.group(3).split(",");
	    				String M[]=m4.group(4).split(",");
	    				String mention_tmp=mention;
	    				for(int i=0;i<type.length;i++)
	    				{
	    					String patt="^(.*?)("+type[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"A");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m4\tType\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<W.length;i++)
	    				{
	    					String patt="^(.*?)("+W[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"W");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m4\tW\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<P.length;i++)
	    				{
	    					P[i]=P[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*?)("+P[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"P");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m4\tP\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<M.length;i++)
	    				{
	    					M[i]=M[i].replace("*", "\\*");
	    					String patt="^(.*?)("+M[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"M");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m4\tM\t"+mention);
	    					}
	    				}
	    			}
	    			else if(m5.find())
	    			{
	    				String type[]=m5.group(1).split(",");
	    				String T[]=m5.group(2).split(",");
	    				String P[]=m5.group(3).split(",");
	    				String M[]=m5.group(4).split(",");
	    				String D[]=m5.group(5).split(",");
	    				String mention_tmp=mention;
	    				for(int i=0;i<type.length;i++)
	    				{
	    					String patt="^(.*?)("+type[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"A");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m5\tType\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<T.length;i++)
	    				{
	    					String patt="^(.*?)("+T[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"T");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m5\tT\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<P.length;i++)
	    				{
	    					P[i]=P[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*?)("+P[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"P");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m5\tP\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<M.length;i++)
	    				{
	    					M[i]=M[i].replace("*", "\\*");
	    					String patt="^(.*?)("+M[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"M");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m5\tM\t"+mention);
	    					}
	    				}
	    				for(int i=0;i<D.length;i++)
	    				{
	    					D[i]=D[i].replaceAll("([^A-Za-z0-9@])", "\\\\$1");
	    					String patt="^(.*?)("+D[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"D");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m5\tD\t"+mention);
	    					}
	    				}
	    			}
	    			else if(m6.find())
	    			{
	    				String RS[]=m6.group(1).split(",");
	    				String mention_tmp=mention;
	    				for(int i=0;i<RS.length;i++)
	    				{
	    					RS[i]=RS[i].replaceAll("([\\[\\]])", "\\\\$1");
	    					String patt="^(.*?)("+RS[i]+")(.*?)$";
	    					Pattern ptmp = Pattern.compile(patt);
	    					Matcher mtmp = ptmp.matcher(mention_tmp);
	    					if(mtmp.find())
	    					{
	    						String mtmp1=mtmp.group(1);
	    						for(int j=mtmp1.length();j<(mtmp1.length()+mtmp.group(2).length());j++)
	    						{
	    							character_hash.put(j,"R");
	    						}
	    						String mtmp2_tmp="";
	    						for(int j=0;j<mtmp.group(2).length();j++){mtmp2_tmp=mtmp2_tmp+"@";}
	    						mention_tmp=mtmp.group(1)+mtmp2_tmp+mtmp.group(3);
	    					}
	    					else
	    					{
	    						System.out.println("Error! Cannot find component: m6\tType\t"+mention);
	    					}
	    				}
	    			}
	    			else
	        		{
	        			System.out.println("Error! Annotation component cannot match RegEx. " + mention);
	        		}
    			}
    			
    			mentionlistbw.write(mention+"\t"+mention2type_hash.get(mention)+"\n");
				if(TrainTest.equals("Train"))
    			{
					mentiondata.write("I I I I I I I I I I I I I I I I I I I\n");
				}
				else
				{
					mentiondata.write("I I I I I I I I I I I I I I I I I I\n");
				}
				String mention_tmp=mention;
				String mention_org=mention;
				mention_tmp = mention_tmp.replaceAll("([0-9])([A-Za-z])", "$1 $2");
				mention_tmp = mention_tmp.replaceAll("([A-Za-z])([0-9])", "$1 $2");
				mention_tmp = mention_tmp.replaceAll("([A-Z])([a-z])", "$1 $2");
				mention_tmp = mention_tmp.replaceAll("([a-z])([A-Z])", "$1 $2");
				mention_tmp = mention_tmp.replaceAll("(.+)fs", "$1 fs");
				mention_tmp = mention_tmp.replaceAll("fs(.+)", "fs $1");
				mention_tmp = mention_tmp.replaceAll("[ ]+", " ");
				String regex="\\s+|(?=\\p{Punct})|(?<=\\p{Punct})";
				String Tokens[]=mention_tmp.split(regex);
				
				start=0;
				last=0;
				
    			for(int i=0;i<Tokens.length;i++)
				{
					if(Tokens[i].length()>0)
					{
						String tkni=Tokens[i].replaceAll("([\\{\\}\\[\\]\\+\\-\\(\\)\\*\\?\\/\\\\])", "\\\\$1");
						Pattern Pcn = Pattern.compile("^([ ]*)("+tkni+")(.*)$");
						Matcher mPcn = Pcn.matcher(mention_org);
		    			if(mPcn.find())
						{
		    				last=last+mPcn.group(1).length();
		    				mention_org=mPcn.group(3);
		    			}
						start=last+1;
						last=start+Tokens[i].length()-1;
					
						//Number of Numbers [0-9]
						String Num_num="";
						String tmp=Tokens[i];
						tmp=tmp.replaceAll("[^0-9]","");
						if(tmp.length()>3){Num_num="N:4+";}else{Num_num="N:"+ tmp.length();}
						
						//Number of Uppercase [A-Z]
						String Num_Uc="";
						tmp=Tokens[i];
						tmp=tmp.replaceAll("[^A-Z]","");
						if(tmp.length()>3){Num_Uc="U:4+";}else{Num_Uc="U:"+ tmp.length();}
						
						//Number of Lowercase [a-z]
						String Num_Lc="";
						tmp=Tokens[i];
						tmp=tmp.replaceAll("[^a-z]","");
						if(tmp.length()>3){Num_Lc="L:4+";}else{Num_Lc="L:"+ tmp.length();}
						
						//Number of ALL char
						String Num_All="";
						if(Tokens[i].length()>3){Num_All="A:4+";}else{Num_All="A:"+ Tokens[i].length();}
						
						//specific character (;:,.->+_)
						String SpecificC="";
						tmp=Tokens[i];
						
						if(Tokens[i].equals(";") || Tokens[i].equals(":") || Tokens[i].equals(",") || Tokens[i].equals(".") || Tokens[i].equals("-") || Tokens[i].equals(">") || Tokens[i].equals("+") || Tokens[i].equals("_"))
						{
							SpecificC="-SpecificC1-";
						}
						else if(Tokens[i].equals("(") || Tokens[i].equals(")"))
						{
							SpecificC="-SpecificC2-";
						}
						else if(Tokens[i].equals("{") || Tokens[i].equals("}"))
						{
							SpecificC="-SpecificC3-";
						}
						else if(Tokens[i].equals("[") || Tokens[i].equals("]"))
						{
							SpecificC="-SpecificC4-";
						}
						else if(Tokens[i].equals("\\") || Tokens[i].equals("/"))
						{
							SpecificC="-SpecificC5-";
						}
						else
						{
							SpecificC="__nil__";
						}
						
						//mutation level
						String Mlevel="";
						if(Tokens[i].equals("p")){Mlevel="-ProteinLevel-";}
						else if(Tokens[i].matches("^[cgmr]$")){Mlevel="-DNALevel-";}
						else{Mlevel="__nil__";}
						
						//mutation type
						String Mtype="";
						String tkn=Tokens[i].toLowerCase();
						String last2_tkn="";
						String last_tkn="";
						String next_tkn="";
						String next2_tkn="";
						if(i>1){last2_tkn=Tokens[i-2];}
						if(i>0){last_tkn=Tokens[i-1];}
						if(Tokens.length>1 && i<Tokens.length-1){next_tkn=Tokens[i+1];}
						if(Tokens.length>2 && i<Tokens.length-2){next2_tkn=Tokens[i+2];}
						
						if(tkn.matches("^(deletion|insertion|duplication|repeat|inversion|delta)$")) {Mtype="-Mtype- -MtypeFull-";}
						else if(tkn.matches("^(eletion|elta|nsertion|uplication|epeat|nversion)$")) {Mtype="-Mtype- -MtypeFull_suffix-";}
						else if(tkn.matches("^(del|ins|delins|indel|dup|inv)$")) {Mtype="-Mtype- -MtypeTri-";}
						else if(tkn.equals("/") && last_tkn.matches("(ins|del)") && next_tkn.toLowerCase().matches("(ins|del)")) {Mtype="-Mtype- -MtypeTri-";}
						else if((last2_tkn.toLowerCase().equals("/") && last_tkn.toLowerCase().equals("\\")) || (next_tkn.toLowerCase().equals("/") && next2_tkn.toLowerCase().equals("\\"))) {Mtype="-Mtype- -MtypeTri-";}
						else {Mtype="__nil__ __nil__";}
						
						//DNA symbols
						String DNASym="";
						if(tkn.matches("^(adenine|guanine|thymine|cytosine)$")){DNASym="-DNASym- -DNASymFull-";}
						else if(tkn.matches("^(denine|uanine|hymine|ytosine)$")){DNASym="-DNASym- -DNASymFull_suffix-";}
						else if(tkn.matches("^[atcgu]+$")){DNASym="-DNASym- -DNASymChar-";}
						else {DNASym="__nil__ __nil__";}
						
						//Protein symbols
						String ProteinSym="";
						if(tkn.matches("^(glutamine|glutamic|leucine|valine|isoleucine|lysine|alanine|glycine|aspartate|methionine|threonine|histidine|aspartic|asparticacid|arginine|asparagine|tryptophan|proline|phenylalanine|cysteine|serine|glutamate|tyrosine|stop|frameshift)$")){ProteinSym="-ProteinSym- -ProteinSymFull-";}
						else if(tkn.matches("^(lutamine|lutamic|eucine|aline|soleucine|ysine|lanine|lycine|spartate|ethionine|hreonine|istidine|spartic|sparticacid|rginine|sparagine|ryptophan|roline|henylalanine|ysteine|erine|lutamate|yrosine|top|rameshift)$")){ProteinSym="-ProteinSym- -ProteinSymFull_suffix-";}
						else if(tkn.matches("^(cys|ile|ser|gln|met|asn|pro|lys|asp|thr|phe|ala|gly|his|leu|arg|trp|val|glu|tyr)$")){ProteinSym="-ProteinSym- -ProteinSymTri-";}
						else if(tkn.matches("^(ys|le|er|ln|et|sn|ro|ys|sp|hr|phe|la|ly|is|eu|rg|rp|al|lu|yr)$")){ProteinSym="-ProteinSym- -ProteinSymTri_suffix-";}
						else if(Tokens[i].matches("^[CISQMNPKDTFAGHLRWVEYX]$") && !next_tkn.toLowerCase().matches("^[ylsrhpiera]")) {ProteinSym="-ProteinSym- -ProteinSymChar-";}
						else if(Tokens[i].matches("^[CISGMPLTHAVF]$") && next_tkn.toLowerCase().matches("^[ylsrhpiera]")) {ProteinSym="-ProteinSym- -ProteinSymChar-";}
						else {ProteinSym="__nil__ __nil__";}
						
						//IVS/EX
						String IVSEX="";
						if(tkn.matches("^(ivs|ex)$")){IVSEX="-IVSEX-";}
						else if(Tokens[i].equals("E") && last_tkn.equals("x")){IVSEX="-IVSEX-";}
						else if(last_tkn.equals("E") && Tokens[i].equals("x")){IVSEX="-IVSEX-";}
						else {IVSEX="__nil__";}
						
						//FSX feature
						String FSXfeature="";
						if(tkn.matches("^(fs|fsx|x|\\*)$")){FSXfeature="-FSX-";}
						else if(last_tkn.toLowerCase().equals("s") && tkn.equals("x")){FSXfeature="-FSX-";}
						else {FSXfeature="__nil__";}
						
						//position type
						String PositionType="";
						if(tkn.matches("^(nucleotide|codon|amino|acid|position|bp|b)$")){PositionType="-PositionType-";}
						else {PositionType="__nil__";}
						
						//sequence location
						String SeqLocat="";
						if(tkn.matches("^(intron|exon|promoter|utr)$")){SeqLocat="-SeqLocat-";}
						else {SeqLocat="__nil__";}
						
						//RS
						String RScode="";
						if(tkn.equals("rs")){RScode="-RScode-";}
						else {RScode="__nil__";}
						
						if(TrainTest.equals("Train"))
		    			{
							mentiondata.write(Tokens[i]+" "+Num_num+" "+Num_Uc+" "+Num_Lc+" "+Num_All+" "+SpecificC+" "+Mlevel+" "+Mtype+" "+DNASym+" "+ProteinSym+" "+IVSEX+" "+FSXfeature+" "+PositionType+" "+SeqLocat+" "+RScode+" "+character_hash.get(start-1)+"\n");
		    			}
						else
						{
							mentiondata.write(Tokens[i]+" "+Num_num+" "+Num_Uc+" "+Num_Lc+" "+Num_All+" "+SpecificC+" "+Mlevel+" "+Mtype+" "+DNASym+" "+ProteinSym+" "+IVSEX+" "+FSXfeature+" "+PositionType+" "+SeqLocat+" "+RScode+"\n");
						}
					}
				}
				mentiondata.write("\n");
			}
			mentionlistbw.close();
			mentiondata.close();
		}
		catch(IOException e1){ System.out.println("[toPostMEData]: "+e1+" Input file is not exist.");}
	}
	
	public void toPostMEModel(String FilenamePostMEdata) throws IOException
	{
		Process process = null;
	    String line = null;
	    InputStream is = null;
	    InputStreamReader isr = null;
	    BufferedReader br = null;
	   
	    Runtime runtime = Runtime.getRuntime();
	    String OS=System.getProperty("os.name").toLowerCase();
		String cmd="";
	    if(OS.contains("windows"))
	    {
	    	cmd ="CRF/crf_learn -f 3 -c 4.0 CRF/template_UB.mention "+FilenamePostMEdata+" CRF/ComponentExtraction.Model.new";
	    }
	    else if(OS.contains("nux")||OS.contains("nix"))
	    {
	    	cmd ="./CRF/crf_learn -f 3 -c 4.0 CRF/template_UB.mention "+FilenamePostMEdata+" CRF/ComponentExtraction.Model.new";
	    }
	    
	    try {
	    	process = runtime.exec(cmd);
		    is = process.getInputStream();
		    isr = new InputStreamReader(is);
		    br = new BufferedReader(isr);
		    while ( (line = br.readLine()) != null) 
		    {
		    	System.out.println(line);
		        System.out.flush();
		    }
		    is.close();
		    isr.close();
		    br.close();
	    }
	    catch (IOException e) {
	    	System.out.println(e);
	    	runtime.exit(0);
	    }
	}
	public void toPostMEoutput(String FilenamePostMEdata,String FilenamePostMEoutput) throws IOException
	{
		/* 
		 * Recognizing components
		 */
		Runtime runtime = Runtime.getRuntime();
	    String OS=System.getProperty("os.name").toLowerCase();
		String cmd="";
	    if(OS.contains("windows"))
	    {
	    	cmd ="CRF/crf_test -m CRF/ComponentExtraction.Model -o "+FilenamePostMEoutput+" "+FilenamePostMEdata;
	    }
	    else if(OS.contains("nux")||OS.contains("nix"))
	    {
	    	cmd ="./CRF/crf_test -m CRF/ComponentExtraction.Model -o "+FilenamePostMEoutput+" "+FilenamePostMEdata;
	    }
	    
	    try {
	    	File f = new File(FilenamePostMEoutput);
	        BufferedWriter fr = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(f), "UTF-8"));
	    	Process process = runtime.exec(cmd);
	    	InputStream is = process.getInputStream();
	    	InputStreamReader isr = new InputStreamReader(is);
	    	BufferedReader br = new BufferedReader(isr);
	    	String line="";
		    while ( (line = br.readLine()) != null) 
		    {
		    	fr.write(line);
		    	fr.newLine();
		        fr.flush();
		    }
		    is.close();
		    isr.close();
		    br.close();
		    fr.close();
	    }
	    catch (IOException e) {
	    	System.out.println(e);
	    	runtime.exit(0);
	    }
	}
	public void output2PubTator(String FilenamePostMEml,String FilenamePostMEoutput,String FilenamePostME,String FilenamePubTator) throws IOException
	{
		try {
			
			/*
			HashMap<String,String> setup_hash = new HashMap<String,String>();
			BufferedReader setupfile = new BufferedReader(new InputStreamReader(new FileInputStream(SetupFile), "UTF-8"));
			String line;
			while ((line = setupfile.readLine()) != null)  
			{
				if(line.contains(" = "))
	        	{
					String str[]=line.split(" = ");//ProteinMutation = True
					setup_hash.put(str[0], str[1]);
	        	}
			}
			setupfile.close();
			*/
			
			ArrayList<String> mentionlist = new ArrayList<String>(); 
			ArrayList<String> identifierlist = new ArrayList<String>(); 
			ArrayList<String> typelist = new ArrayList<String>(); 
			BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenamePostMEml), "UTF-8"));
			int count=0;
			HashMap<Integer,Integer> boundary_hash = new HashMap<Integer,Integer>();
			HashMap<Integer,String> WMstate_hash = new HashMap<Integer,String>();
			String line;
			while ((line = inputfile.readLine()) != null)  
			{
				String columns[]=line.split("\\t");
				String line_nospace=columns[0].replaceAll(" ","");
				Pattern pat1 = Pattern.compile("(.+?)(for|inplaceof|insteadof|mutantof|mutantsof|ratherthan|ratherthan|replacementof|replacementsof|replaces|replacing|residueatposition|residuefor|residueinplaceof|residueinsteadof|substitutioat|substitutionfor|substitutionof|substitutionsat|substitutionsfor|substitutionsof|substitutedfor|toreplace)(.+)$");
				Matcher mat1 = pat1.matcher(line_nospace.toLowerCase());
				Pattern pat2 = Pattern.compile("^(.+?)(>|to|into|of|by|with|at)(.+)$");
				Matcher mat2 = pat2.matcher(line_nospace.toLowerCase());
				
				if(mat1.find())
				{
					boundary_hash.put(count,mat1.group(1).length());
					WMstate_hash.put(count,"Backward");
				}
				else if(mat2.find())
				{
					boundary_hash.put(count,mat2.group(1).length());
					WMstate_hash.put(count,"Forward");
				}
				mentionlist.add(columns[0]);
				if(columns.length==2)
				{
					typelist.add(columns[1]);
				}
				else
				{
					typelist.add("DNAMutation");
				}
				count++;
			}
			inputfile.close();
			
			HashMap<String,String> component_hash = new HashMap<String,String>();
			component_hash.put("A", ""); //type
			component_hash.put("T", ""); //Method
			component_hash.put("P", ""); //Position
			component_hash.put("W", ""); //Wide type
			component_hash.put("M", ""); //Mutant
			component_hash.put("F", ""); //frame shift
			component_hash.put("S", ""); //frame shift position
			component_hash.put("D", ""); //
			component_hash.put("I", ""); //Inside
			component_hash.put("R", ""); //RS number
			
			BufferedReader PostMEfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenamePostMEoutput), "UTF-8"));
			String prestate="";
			count=0;
			int start_count=0;
			HashMap<Integer,String> filteringNum_hash = new HashMap<Integer,String>();
			while ((line = PostMEfile.readLine()) != null)  
			{
				String outputs[]=line.split("\\t");
					
				/*
				 * recognize status and mention
				 */
				if(outputs.length<=1)
				{
					/*
					 *  Translate : nametothree | threetone | etc
					 */
					component_hash.put("W",component_hash.get("W").toUpperCase());
					component_hash.put("M",component_hash.get("M").toUpperCase());
					component_hash.put("P",component_hash.get("P").toUpperCase());
					HashMap<String,String> nametothree = new HashMap<String,String>();
					nametothree.put("THYMINE", "T");nametothree.put("ALANINE", "ALA");nametothree.put("ARGININE", "ARG");nametothree.put("ASPARAGINE", "ASN");nametothree.put("ASPARTICACID", "ASP");nametothree.put("ASPARTATE", "ASP");nametothree.put("CYSTEINE", "CYS");nametothree.put("GLUTAMINE", "GLN");nametothree.put("GLUTAMICACID", "GLU");nametothree.put("GLUTAMATE", "GLU");nametothree.put("GLYCINE", "GLY");nametothree.put("HISTIDINE", "HIS");nametothree.put("ISOLEUCINE", "ILE");nametothree.put("LEUCINE", "LEU");nametothree.put("LYSINE", "LYS");nametothree.put("METHIONINE", "MET");nametothree.put("PHENYLALANINE", "PHE");nametothree.put("PROLINE", "PRO");nametothree.put("SERINE", "SER");nametothree.put("THREONINE", "THR");nametothree.put("TRYPTOPHAN", "TRP");nametothree.put("TYROSINE", "TYR");nametothree.put("VALINE", "VAL");nametothree.put("STOP", "XAA");
					HashMap<String,String> threetone = new HashMap<String,String>();
					threetone.put("ALA", "A");threetone.put("ARG", "R");threetone.put("ASN", "N");threetone.put("ASP", "D");threetone.put("CYS", "C");threetone.put("GLN", "Q");threetone.put("GLU", "E");threetone.put("GLY", "G");threetone.put("HIS", "H");threetone.put("ILE", "I");threetone.put("LEU", "L");threetone.put("LYS", "K");threetone.put("MET", "M");threetone.put("PHE", "F");threetone.put("PRO", "P");threetone.put("SER", "S");threetone.put("THR", "T");threetone.put("TRP", "W");threetone.put("TYR", "Y");threetone.put("VAL", "V");threetone.put("ASX", "B");threetone.put("GLX", "Z");threetone.put("XAA", "X");threetone.put("TER", "X");
					boolean translate=false;
					boolean NotAaminoacid=true;
					
					//M
					String components[]=component_hash.get("M").split(",");
					component_hash.put("M","");
					String component="";
					for(int i=0;i<components.length;i++)
					{
						if(nametothree.containsKey(components[i]))
						{
							component=nametothree.get(components[i]);
							translate=true;
						}
						else
						{
							component=components[i];
						}	
						if(component_hash.get("M").equals(""))
						{
							component_hash.put("M",component);
						}
						else
						{
							component_hash.put("M",component_hash.get("M")+","+component);
						}
					}
					String components2[]=component_hash.get("M").split(",");
					component_hash.put("M","");
					component="";
					for(int i=0;i<components2.length;i++)
					{
						if(threetone.containsKey(components2[i]))
						{
							component=threetone.get(components2[i]);
							translate=true;
						}
						else if(components2[i].length()>1)
						{
							NotAaminoacid=false;
							component=components2[i];
						}
						else
						{
							component=components2[i];
						}	
						if(component_hash.get("M").equals(""))
						{
							component_hash.put("M",component);
						}
						else
						{
							component_hash.put("M",component_hash.get("M")+","+component);
						}
					}
					
					//W
					String components3[]=component_hash.get("W").split(",");
					component_hash.put("W","");
					component="";
					for(int i=0;i<components3.length;i++)
					{
						if(nametothree.containsKey(components3[i]))
						{
							component=nametothree.get(components3[i]);
							translate=true;
						}
						else
						{
							component=components3[i];
						}	
						if(component_hash.get("W").equals(""))
						{
							component_hash.put("W",component);
						}
						else
						{
							component_hash.put("W",component_hash.get("W")+","+component);
						}
					}
					String components4[]=component_hash.get("W").split(",");
					component_hash.put("W","");
					component="";
					for(int i=0;i<components4.length;i++)
					{
						if(threetone.containsKey(components4[i]))
						{
							component=threetone.get(components4[i]);
							translate=true;
						}
						else if(components4[i].length()>1)
						{
							NotAaminoacid=false;
							component=components4[i];
						}
						else
						{
							component=components4[i];
						}	
						if(component_hash.get("W").equals(""))
						{
							component_hash.put("W",component);
						}
						else
						{
							component_hash.put("W",component_hash.get("W")+","+component);
						}
					}
					
					//W/M - 2
					if(component_hash.get("W").matches(",") && (component_hash.get("M").equals("") && component_hash.get("T").equals("")))
					{
						String spl[]=component_hash.get("W").split("-");
						component_hash.put("W",spl[0]);
						component_hash.put("M",spl[1]);
					}
					if(component_hash.get("M").matches(",") && (component_hash.get("W").equals("") && component_hash.get("T").equals("")))
					{
						String spl[]=component_hash.get("M").split("-");
						component_hash.put("W",spl[0]);
						component_hash.put("M",spl[1]);
					}
					
					//A
					Pattern pap = Pattern.compile("[+-][0-9]");
					Matcher mp = pap.matcher(component_hash.get("P"));
					component_hash.put("A",component_hash.get("A").toLowerCase());
					if(component_hash.get("A").equals("") && component_hash.get("P").matches("[\\+]*[0-9]+") && (component_hash.get("P").matches("[\\-\\+]{0,1}[0-9]{1,8}")  && component_hash.get("W").matches("[ATCG]")  && component_hash.get("M").matches("[ATCG]") && Integer.parseInt(component_hash.get("P"))>3000)) 
					{
						component_hash.put("A","g");
					}
					else if(component_hash.get("A").equals("cdna"))
					{
						component_hash.put("A","c");
					}
					else if(component_hash.get("A").equals("") && mp.find() && !typelist.get(count).equals("ProteinMutation"))
					{
						component_hash.put("A","c");
					}
					else if(component_hash.get("W").matches("^[ATCG]*$") && component_hash.get("M").matches("^[ATCG]*$") && translate==false && component_hash.get("A").equals(""))
					{
						component_hash.put("A","c");
					}
					
					//F
					if(component_hash.get("F").equals("*"))
					{
						component_hash.put("F","X");
					}
					
					
					//R
					component_hash.put("R",component_hash.get("R").toLowerCase());
					component_hash.put("R",component_hash.get("R").replaceAll("[\\[\\]]", ""));
					
					//P
					if(component_hash.get("P").matches("^([0-9]+)-([0-9]{1,8})$"))
					{
						String spl[]=component_hash.get("P").split("-");
						if(Integer.parseInt(spl[0])<Integer.parseInt(spl[1]))
						{
							component_hash.put("P",spl[0]+"_"+spl[1]);
						}
					}
					component_hash.put("P",component_hash.get("P").replaceAll("[-\\[]+$", ""));
					component_hash.put("P",component_hash.get("P").replaceAll("^(POSITION|NUCLEOTIDE|:)", ""));
					if(typelist.get(count).equals("ProteinMutation"))
					{
						component_hash.put("P",component_hash.get("P").replaceAll("^CODON", ""));	
					}
					
					//T
					component_hash.put("T",component_hash.get("T").toUpperCase());
					component_hash.put("T",component_hash.get("T").replaceAll("DELTA","DEL"));
					component_hash.put("T",component_hash.get("T").replaceAll("INSERTION","INS"));
					if(component_hash.get("T").equals("DELINS")) {component_hash.put("T","INDEL"); translate=false;}
					else if(component_hash.get("T").matches("INS.*DEL")) {component_hash.put("T","INDEL"); translate=false;}
					else if(component_hash.get("T").matches("DEL.*INS")) {component_hash.put("T","INDEL"); translate=false;}
					if(component_hash.get("M").matches("(DEL|INS|DUP)"))
					{
						component_hash.put("T",component_hash.get("M"));
						component_hash.put("M","");
					}
					else if(component_hash.get("W").matches("(DEL|INS|DUP)"))
					{
						component_hash.put("T",component_hash.get("W"));
						component_hash.put("W","");
					}
					else if(!component_hash.get("D").equals(""))
					{
						component_hash.put("T","DUP");
					}
					
					// Recognize the components(identifiers)
					String identifier="";
					String type="";
					if(component_hash.get("T").equals("DUP") && !component_hash.get("D").equals("")) //dup
					{
						identifier=component_hash.get("A")+"|"+component_hash.get("T")+"|"+component_hash.get("P")+"|"+component_hash.get("M")+"|"+component_hash.get("D");
					}
					else if(component_hash.get("T").equals("DUP")) //dup
					{
						identifier=component_hash.get("A")+"|"+component_hash.get("T")+"|"+component_hash.get("P")+"|"+component_hash.get("M")+"|";
					}
					else if(!component_hash.get("T").equals("")) //DEL|INS|INDEL
					{
						identifier=component_hash.get("A")+"|"+component_hash.get("T")+"|"+component_hash.get("P")+"|"+component_hash.get("M");
					}
					else if(!component_hash.get("F").equals("")) //FS
					{
						identifier=component_hash.get("A")+"|FS|"+component_hash.get("W")+"|"+component_hash.get("P")+"|"+component_hash.get("M")+"|"+component_hash.get("S");
						type="ProteinMutation";
					}
					else if(!component_hash.get("R").equals("")) //RS
					{
						identifier=mentionlist.get(count);
						type="SNP";
					}
					else if(mentionlist.get(count).matches("^I([RSrs][Ss][0-9].+)"))
					{
						String men=mentionlist.get(count);
						men.substring(1, men.length());
						men=men.replaceAll("[\\W-_]", men.toLowerCase());
						type="SNP";
					}
					else
					{
						identifier=component_hash.get("A")+"|SUB|"+component_hash.get("W")+"|"+component_hash.get("P")+"|"+component_hash.get("M");
					}
					
					// filteringNum_hash
					if(component_hash.get("M").equals(component_hash.get("W")) && (!component_hash.get("W").equals("")) && type.equals("DNAMutation")) // remove genotype
					{
						filteringNum_hash.put(count, "");
					}
					else if(NotAaminoacid==false && type.equals("ProteinMutation")) //E 243 ASPARTATE
					{
						filteringNum_hash.put(count, "");
					}
					else if(component_hash.get("W").matches(",") && component_hash.get("M").matches(",") && component_hash.get("P").matches("")) //T,C/T,C
					{
						filteringNum_hash.put(count, "");
					}
					else if(component_hash.get("W").equals("") && component_hash.get("M").equals("") && component_hash.get("T").equals("") && !type.equals("SNP")) //exons 5
					{
						filteringNum_hash.put(count, "");
					}
					else if(mentionlist.get(count).matches("^I[RSrs][Ss]") && type.equals("SNP")) //exons 5
					{
						filteringNum_hash.put(count, "");
					}
					
					
					// Recognize the type of mentions: ProteinMutation | DNAMutation | SNP 
					if(type.equals(""))
					{
						if(translate==true)
						{
							type="ProteinMutation";
						}
						else
						{
							type="DNAMutation";
						}
						
						if(component_hash.get("P").matches("^([Ee]x|EX|[In]ntron|IVS|[Ii]vs)"))
						{
							type="DNAMutation";
							identifier=identifier.replaceAll("^[^|]*\\|","c|");
						}
						else if(component_hash.get("P").matches(".*[\\+\\-].*"))
						{
							type="DNAMutation";
							identifier=identifier.replaceAll("^[^|]*\\|","c|");
						}
						else if(component_hash.get("M").matches("[ISQMNPKDFHLRWVEYX]") || component_hash.get("W").matches("[ISQMNPKDFHLRWVEYX]") )
						{
							type="ProteinMutation";
							identifier=identifier.replaceAll("^[^|]*\\|","p|");
						}
						else if(type.equals(""))
						{
							type="DNAMutation";
						}
						
						if(!component_hash.get("A").equals("") &&
						   (
						   component_hash.get("A").toLowerCase().equals("c") ||
						   component_hash.get("A").toLowerCase().equals("r") ||
						   component_hash.get("A").toLowerCase().equals("m") ||
						   component_hash.get("A").toLowerCase().equals("g")
						   )
						  )
						{
							type="DNAMutation";
							identifier=identifier.replaceAll("^[^|]*\\|",component_hash.get("A")+"|");
						}
						else if(component_hash.get("A").equals("p"))
						{
							type="ProteinMutation";
						}
					}
					
					if(type.equals("ProteinMutation"))
					{
						identifier=identifier.replaceAll("^[^|]*\\|","p|");
					}
					
					//System.out.println(mentionlist.get(count)+"\t"+component_hash.get("T")+"\t"+component_hash.get("M")+"\t"+component_hash.get("W"));
					
					// filtering and Print out
					if( (component_hash.get("W").length() == 3 || component_hash.get("M").length() == 3) 
							&& component_hash.get("W").length() != component_hash.get("M").length()
							&& !component_hash.get("W").equals("") && !component_hash.get("M").equals("") && component_hash.get("W").indexOf(",")!=-1 && component_hash.get("M").indexOf(",")!=-1
							&& ((component_hash.get("W").indexOf("A")!=-1 && component_hash.get("W").indexOf("T")!=-1 && component_hash.get("W").indexOf("C")!=-1 && component_hash.get("W").indexOf("G")!=-1) || (component_hash.get("M").indexOf("A")!=-1 && component_hash.get("M").indexOf("T")!=-1 && component_hash.get("M").indexOf("C")!=-1 && component_hash.get("M").indexOf("G")!=-1))
							&& component_hash.get("T").equals("")
						)
					{/*System.out.println("filtering 1:"+mentionlist.get(count));*/identifierlist.add("_Remove_");}
					else if((component_hash.get("M").matches("[ISQMNPKDFHLRWVEYX]") || component_hash.get("W").matches("[ISQMNPKDFHLRWVEYX]")) && component_hash.get("P").matches("[6-9][0-9][0-9][0-9]+")){/*System.out.println("filtering length:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //M300000X
					else if(component_hash.get("M").equals("") && component_hash.get("T").equals("") && component_hash.get("F").equals("") && component_hash.get("D").equals("") && component_hash.get("P").equals("") && !type.equals("SNP")){/*System.out.println("filtering 2:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //Arg235
					else if(component_hash.get("T").equals("DUP") && component_hash.get("M").matches("") && component_hash.get("P").equals("") && !type.equals("SNP")){/*System.out.println("filtering 2:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //|DUP|||4]q33-->qter	del[4] q33-->qter	del4q33qter
					else if(component_hash.get("T").equals("") && (component_hash.get("M").equals("") || component_hash.get("W").equals("")) && component_hash.get("F").equals("") && component_hash.get("D").equals("") && !type.equals("SNP")){/*System.out.println("filtering 3:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //g. 200 A
					else if(component_hash.get("T").equals("") && component_hash.get("M").matches("[ATCGUatcgu]") && (component_hash.get("M").equals(component_hash.get("W")))){/*System.out.println("filtering 4:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //T --> T
					else if(component_hash.get("P").matches("-[0-9]+") && type.equals("ProteinMutation")){/*System.out.println("filtering 6:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //negative protein mutation
					else if(component_hash.get("P").matches(".*>.*")){/*System.out.println("filtering 6:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //negative protein mutation
					else if(component_hash.get("W").matches("^[BJOUZ]") || component_hash.get("M").matches("^[BJOUZ]")){/*System.out.println("filtering 7:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //not a mutation
					else if(type.equals("SNP") && identifier.matches("RS[0-9][0-9]{0,1}")){/*System.out.println("filtering 2:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //too short rs number
					else if(type.equals("SNP") && identifier.matches("[0-9][0-9]{0,1}[\\W\\-\\_]*delta")){/*System.out.println("filtering 2:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //0 delta
					else if(tmVar.PAM_lowerScorePair.contains(component_hash.get("M")+"\t"+component_hash.get("W"))){/*System.out.println("filtering 7:"+mentionlist.get(count));*/identifierlist.add("_Remove_");} //unlikely to occur
					else
					{
						identifierlist.add(identifier);
						typelist.add(type);
					}
					
					//End
					component_hash.put("A", "");
					component_hash.put("T", "");
					component_hash.put("P", "");
					component_hash.put("W", "");
					component_hash.put("M", "");
					component_hash.put("F", "");
					component_hash.put("S", "");
					component_hash.put("D", "");
					component_hash.put("I", "");
					component_hash.put("R", "");
					count++;
					start_count=0;
				}
				else if(outputs[1].equals("I"))
				{
					//Start
				}
				else if(outputs[outputs.length-1].matches("[ATPWMFSDIR]"))
				{
					if(WMstate_hash.containsKey(count) && WMstate_hash.get(count).equals("Forward"))
					{
						if(start_count<boundary_hash.get(count) && outputs[outputs.length-1].equals("M"))
						{
							outputs[outputs.length-1]="W";
						}
						else if(start_count>boundary_hash.get(count) && outputs[outputs.length-1].equals("W"))
						{
							outputs[outputs.length-1]="M";
						}
					}
					else if(WMstate_hash.containsKey(count) && WMstate_hash.get(count).equals("Backward"))
					{
						if(start_count<boundary_hash.get(count) && outputs[outputs.length-1].equals("W"))
						{
							outputs[outputs.length-1]="M";
						}
						else if(start_count>boundary_hash.get(count) && outputs[outputs.length-1].equals("M"))
						{
							outputs[outputs.length-1]="W";
						}
					}
					String state=outputs[outputs.length-1];	
					String tkn=outputs[0];
					
					if(!component_hash.get(state).equals("") && !state.equals(prestate))
					{
						component_hash.put(state, component_hash.get(state)+","+tkn);
					}
					else
					{
						component_hash.put(state, component_hash.get(state)+tkn);
					}
					prestate=state;
				}
				start_count=start_count+outputs[0].length();
			}
			PostMEfile.close();
			
			HashMap<String,String> mention2id_hash = new HashMap<String,String>();
			for(int i=0;i<count;i++)
			{
				if(!filteringNum_hash.containsKey(i) && (!identifierlist.get(i).equals("_Remove_")))
				{
					mention2id_hash.put(mentionlist.get(i),typelist.get(i)+"\t"+identifierlist.get(i));
				}
			}
			
			//Filtering
			HashMap<String,String> filteringStr_hash = new HashMap<String,String>();
			BufferedReader filterfile = new BufferedReader(new InputStreamReader(new FileInputStream("lib/filtering.txt"), "UTF-8"));
			while ((line = filterfile.readLine()) != null)  
			{
				filteringStr_hash.put(line, "");
			}
			filterfile.close();
			
			BufferedWriter PubTatorfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(FilenamePubTator), "UTF-8")); // .location
			PostMEfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenamePostME), "UTF-8"));
			while ((line = PostMEfile.readLine()) != null)  
			{
				Pattern pat = Pattern.compile("^([^\\|\\t]+)\\|([^\\|\\t]+)\\|(.*)$");
				Matcher mat = pat.matcher(line);
				if(mat.find()) //Title|Abstract
	        	{
					PubTatorfile.write(line+"\n");
	        	}
				else if(line.contains("\t")) //Annotation
	        	{
					String outputs[]=line.split("\\t");
					String pmid=outputs[0];
					String start=outputs[1];
					String last=outputs[2];
					String mention=outputs[3];
					if((!filteringStr_hash.containsKey(mention)) && mention2id_hash.containsKey(mention))
					{
						if(mention.matches(".* at")){}
						else if(mention.matches("C [0-9]+H")){}
						else
						{
							mention2id_hash.put(mention,mention2id_hash.get(mention).replaceAll(" ",""));
							PubTatorfile.write(pmid+"\t"+start+"\t"+last+"\t"+mention+"\t"+mention2id_hash.get(mention)+"\n");
						}
					}
	        	}
				else
				{
					PubTatorfile.write(line+"\n");
				}
			}
			PostMEfile.close();
			PubTatorfile.close();
		}
		catch(IOException e1){ System.out.println("[output2PubTator]: "+e1+" Input file is not exist.");}
	}
	
	public void Normalization(String DisplayRSnumOnly, String input,String outputPubTator,String finalPubTator) throws IOException
	{
		/**
		 * input : gene mentions
		 * outputPubTator : mutation mentions
		 * finalPubTator : normalized result of mutation mentions
		 */
		BreakIterator iterator = BreakIterator.getSentenceInstance(Locale.US);	
		
		/*Database Connection*/
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
		
		/*one2three*/
		HashMap<String,String> one2three = new HashMap<String,String>();
		one2three.put("A", "Ala");
		one2three.put("R", "Arg");
		one2three.put("N", "Asn");
		one2three.put("D", "Asp");
		one2three.put("C", "Cys");
		one2three.put("Q", "Gln");
		one2three.put("E", "Glu");
		one2three.put("G", "Gly");
		one2three.put("H", "His");
		one2three.put("I", "Ile");
		one2three.put("L", "Leu");
		one2three.put("K", "Lys");
		one2three.put("M", "Met");
		one2three.put("F", "Phe");
		one2three.put("P", "Pro");
		one2three.put("S", "Ser");
		one2three.put("T", "Thr");
		one2three.put("W", "Trp");
		one2three.put("Y", "Tyr");
		one2three.put("V", "Val");
		one2three.put("B", "Asx");
		one2three.put("Z", "Glx");
		one2three.put("X", "Xaa");
		one2three.put("X", "Ter");
		
		/*RS_DNA_Protein.txt - Pattern*/
		ArrayList<String> RS_DNA_Protein = new ArrayList<String>();
		BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream("lib/RegEx/RS_DNA_Protein.txt"), "UTF-8"));
		String line;
		while ((line = inputfile.readLine()) != null)  
		{
			if(!line.equals(""))
			{
				RS_DNA_Protein.add(line);
			}
		}
		inputfile.close();
		
		/*Mutation_RS_Geneid.txt*/
		HashMap<String,String> Mutation_RS_Geneid_hash = new HashMap<String,String>();
		inputfile = new BufferedReader(new InputStreamReader(new FileInputStream("lib/RegEx/Mutation_RS_Geneid.txt"), "UTF-8"));
		while ((line = inputfile.readLine()) != null)  
		{
			Pattern pat = Pattern.compile("^(Pattern|Recognized)	([^\t]+)	([0-9]+)	([0-9,]+)$");
			Matcher mat = pat.matcher(line);
			if (mat.find())
        	{
				String geneids[]=mat.group(4).split(",");
				for(int i=0;i<geneids.length;i++)
				{
					//mutation id | geneid --> rs#
					Mutation_RS_Geneid_hash.put(mat.group(2)+"\t"+geneids[i], mat.group(3));
				}
        	}
		}
		inputfile.close();
		
		
		try {
			/*
			 * Gene Mention Extraction
			 */
			inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(input), "UTF-8"));
			HashMap<String,String> tmp_gene = new HashMap<String,String>();
			HashMap<String,String> annotations_gene = new HashMap<String,String>();
			while ((line = inputfile.readLine()) != null)  
			{
				Pattern pat = Pattern.compile("^([^\\|\\t]+)\\|([^\\|\\t]+)\\|(.*)$");
				Matcher mat = pat.matcher(line);
				if(mat.find()) //Title|Abstract
	        	{
					
	        	}
				else if (line.contains("\t")) //Annotation
	        	{
					String anno[]=line.split("\t");
	        		if(anno.length>=6)
	        		{
	        			String mentiontype=anno[4];
	        			if(mentiontype.equals("Gene"))
	        			{
	        				/*
	        				 * Search Database - Gene
	        				 */
	        				String geneids[]=anno[5].split(";");
	        				for(int gi=0;gi<geneids.length;gi++)
	        				{
	        					annotations_gene.put(anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+geneids[gi],"");
		        				
	        					/* tmp_gene(pmid,geneid) -> RS#s */
	        					/* annotations_gene(pmid,start,last,mention,geneid) -> RS#s */
	        					if(tmp_gene.containsKey(anno[0]+"\t"+geneids[gi]))
		        				{
		        					annotations_gene.put(anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+geneids[gi],tmp_gene.get(anno[0]+"\t"+geneids[gi]));
		        				}
		        				else
		        				{
			        				try {
			        					c = DriverManager.getConnection("jdbc:sqlite:Database/gene2rs.db");
			        					stmt = c.createStatement();
			        					ResultSet rs = stmt.executeQuery("SELECT rs FROM gene2rs WHERE gene='"+geneids[gi]+"';");
			        					if ( rs.next() ) 
			        					{
			        				         String  rsNumber = rs.getString("rs");
			        				         annotations_gene.put(anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+geneids[gi],rsNumber);
			        				         tmp_gene.put(anno[0]+"\t"+geneids[gi],rsNumber);
			        				    }
		        						stmt.close();
		        						c.close();
		        					} 
		        					catch ( SQLException e ) 
		        					{
		        						System.err.println( e.getClass().getName() + ": " + e.getMessage() );
		        					}
		        				}
	        				}
	        			}
	        		}
	        	}
			}
			inputfile.close();
			
			/*
			 * SNP(RSnumberList) Extraction
			 */
			String article="";
			String Pmid="";
			HashMap<String,String> RSnumberList = new HashMap<String,String>();
			HashMap<Integer,String> start2id = new HashMap<Integer,String>(); //tmp
			HashMap<String,String> MentionPatternMap2rs = new HashMap<String,String>(); //by pattern
			HashMap<String,String> MentionPatternMap = new HashMap<String,String>(); 
			HashMap<String,String> MentionPatternMap2rs_extend = new HashMap<String,String>(); //by pattern [PMID:24073655] c.169G>C (p.Gly57Arg) -- map to rs764352037
			inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(outputPubTator), "UTF-8"));
			
			HashMap<String,ArrayList<Integer>> PMID2SentenceOffsets = new HashMap<String,ArrayList<Integer>>();
			
			while ((line = inputfile.readLine()) != null)  
			{
				Pattern pat = Pattern.compile("^([^\\|\\t]+)\\|([^\\|\\t]+)\\|(.*)$");
				Matcher mat = pat.matcher(line);
				if(mat.find()) //Title|Abstract
	        	{
					Pmid = mat.group(1);
					String ParagraphContent=mat.group(3);
					article=article+ParagraphContent+"\t";
					
					// split sentences
					iterator.setText(article);
					ArrayList<Integer> Sentence_offsets = new ArrayList<Integer>();
					int Sent_start = iterator.first();
					for (int Sent_last = iterator.next(); Sent_last != BreakIterator.DONE; Sent_start = Sent_last, Sent_last = iterator.next()) 
					{
						Sentence_offsets.add(Sent_start);
					}
					PMID2SentenceOffsets.put(Pmid, Sentence_offsets);
				}
				else if (line.contains("\t")) //Annotation
		    	{
					String anno[]=line.split("\t");
	        		if(anno.length>=6)
	        		{
	        			Pmid = anno[0];
	        			int start = Integer.parseInt(anno[1]);
    	        		int last = Integer.parseInt(anno[2]);
    	        		String mention=anno[3];
    	        		String mentiontype=anno[4];
    	        		String id=anno[5];
    	        		if(anno[5].matches("\\|.*")){anno[5]="p"+anno[5];}
    	        		
    	        		if(mentiontype.equals("DNAMutation"))
	        			{
    	        			if(article.substring(start,last).equals(mention))
        	        		{
        	        			String tmp = mention;
        	        			tmp=tmp.replaceAll(".","D");
        	        			article=article.substring(0,start)+tmp+article.substring(last,article.length());
        	        		}
    	        			start2id.put(start,id);
	        			}
    	        		else if(mentiontype.equals("ProteinMutation"))
	        			{
    	        			if(article.substring(start,last).equals(mention))
        	        		{
        	        			String tmp = mention;
        	        			tmp=tmp.replaceAll(".","P");
        	        			article=article.substring(0,start)+tmp+article.substring(last,article.length());
        	        		}
    	        			start2id.put(start,id);
	        			}
    	        		else if(mentiontype.equals("SNP"))
	        			{
	        				anno[5]=anno[5].replaceAll("^rs","");
	        				RSnumberList.put(anno[0]+"\t"+anno[5],"");
	        				
	        				if(article.substring(start,last).equals(mention))
        	        		{
        	        			String tmp = mention;
        	        			tmp=tmp.replaceAll(".","S");
        	        			article=article.substring(0,start)+tmp+article.substring(last,article.length());
        	        		}
	        				start2id.put(start,id);
	        			}
	        		}
	        	}
				else if(line.length()==0)
				{
					String text=article;
					for(int i=0;i<RS_DNA_Protein.size();i++)
					{
						pat = Pattern.compile("^(.*?)("+RS_DNA_Protein.get(i)+")");
						mat = pat.matcher(text);
						while(mat.find())
						{
							int start=mat.group(1).length();
							int last=start+mat.group(2).length();
							
							ArrayList<String> DP = new ArrayList<String>();
							String S = ""; // the RS# in the pattern
							
							for(int s : start2id.keySet())
							{
								if(s>=start && s<last)
								{
									Pattern pat_rs = Pattern.compile("[Rr][Ss]([0-9]+)$");
									Matcher mat_rs = pat_rs.matcher(start2id.get(s));
									if(mat_rs.find())
									{
										S = mat_rs.group(1);
									}
									else
									{
										DP.add(start2id.get(s));
									}
								}
							}
							if(!S.equals("")) //RS number is in the pattern
							{
								for(int dp=0;dp<DP.size();dp++)
								{
									MentionPatternMap2rs.put(Pmid+"\t"+DP.get(dp),S);
								}
							}
							else 
							{
								if(DP.size()>1)
								{
									MentionPatternMap.put(Pmid+"\t"+DP.get(0),DP.get(1));
									MentionPatternMap.put(Pmid+"\t"+DP.get(1),DP.get(0));
								}
							}
							
							String tmp = mat.group(2);
							tmp=tmp.replaceAll(".", "F");
							text = text.substring(0, start)+tmp+text.substring(last,text.length());
							pat = Pattern.compile("^(.*?)("+RS_DNA_Protein.get(i)+")");
							mat = pat.matcher(text);
						}
					}
					
					article="";
					start2id.clear();
				}
			}
			inputfile.close();
			
			
			/*
			 * Mutation Extraction
			 */
			BufferedWriter outputfile = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(finalPubTator), "UTF-8")); // .location
			inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(outputPubTator), "UTF-8"));
			HashMap<String,String> tmp_mutation = new HashMap<String,String>();
			HashMap<String,String> annotations_mutation = new HashMap<String,String>();
			HashMap<String,String> rs_foundbefore_hash = new HashMap<String,String>();
			HashMap<String,String> rs_foundinText_hash = new HashMap<String,String>();
			String outputSTR="";
			
			while ((line = inputfile.readLine()) != null)  
			{
				Pattern pat = Pattern.compile("^([^\\|\\t]+)\\|([^\\|\\t]+)\\|(.*)$");
				Matcher mat = pat.matcher(line);
				if(mat.find()) //Title|Abstract
	        	{
					outputSTR=outputSTR+line+"\n";
	        	}
				else if (line.contains("\t")) //Annotation
	        	{
					String anno[]=line.split("\t");
	        		if(anno.length>=6)
	        		{
	        			String mentiontype=anno[4];
	        			if(mentiontype.matches("(DNAMutation|ProteinMutation)"))
	        			{
	        				if(anno[5].matches("\\|.*")){anno[5]="p"+anno[5];}
	        				
	        				String component[]=anno[5].split("\\|",-1);
	        				String NormalizedForm="";
	        				String NormalizedForm_reverse="";
	        				String NormalizedForm_plus1="";
	        				String NormalizedForm_minus1="";
	        				String NormalizedForm_protein="";
	        				if(component.length>=3)
	        				{
		        				if(component[1].equals("SUB"))
		        				{
		        					if(component[0].equals("p"))
			        				{
		        						String tmp="";
		        						/*one -> three*/
		        						for(int len=0;len<component[2].length();len++)
		        						{
		        							if(one2three.containsKey(component[2].substring(len, len+1)))
		        							{
		        								if(tmp.equals(""))
		        								{
		        									tmp = one2three.get(component[2].substring(len, len+1));
		        								}
		        								else
		        								{
		        									tmp = tmp +","+ one2three.get(component[2].substring(len, len+1));
		        								}
		        							}
		        						}
		        						component[2]=tmp;
		        						tmp="";
		        						for(int len=0;len<component[4].length();len++)
		        						{
		        							if(one2three.containsKey(component[4].substring(len, len+1)))
		        							{
		        								if(tmp.equals(""))
		        								{
		        									tmp = one2three.get(component[4].substring(len, len+1));	
		        								}
		        								else
		        								{	     
		        									tmp = tmp +","+ one2three.get(component[4].substring(len, len+1));	
		        								}
		        							}
		        						}
		        						component[4]=tmp;
		        						
		        						if(component[2].equals(component[4]))
		        						{
		        							NormalizedForm=component[0]+"."+component[2]+component[3]+"=";
		        						}
		        						else
		        						{
			        						NormalizedForm=component[0]+"."+component[2]+component[3]+component[4];
			        						NormalizedForm_reverse=component[0]+"."+component[4]+component[3]+component[2];
		        						}
			        					String wildtype[]=component[2].split(",");
			        					String mutant[]=component[4].split(",");
			        					String positions[]=component[3].split(",");
			        					
			        					if(wildtype.length == positions.length && wildtype.length>1) //Pair of wildtype&position
			        					{
			        						for (int i=0;i<wildtype.length;i++) //May have more than one pair
				        					{
				        						for (int j=0;j<mutant.length;j++)
					        					{
				        							NormalizedForm=NormalizedForm+"|"+component[0]+"."+wildtype[i]+positions[i]+mutant[j];
					        					}
				        					}
			        					}
			        					else
			        					{
			        						for (int i=0;i<wildtype.length;i++) //May have more than one pair
				        					{
				        						for (int j=0;j<mutant.length;j++)
					        					{
				        							NormalizedForm=NormalizedForm+"|"+component[0]+"."+wildtype[i]+component[3]+mutant[j];
					        					}
				        					}
			        					}
			        				}
		        					else //[rmgc]
		        					{
		        						if(component.length>4)
		    	        				{
			        						component[3]=component[3].replaceAll("^\\+", "");
			        						NormalizedForm=component[0]+"."+component[3]+component[2]+">"+component[4];
			        						NormalizedForm_reverse=component[0]+"."+component[3]+component[4]+">"+component[2];
			        						if(component[3].matches("[\\+\\-]{0,1}[0-9]{1,8}"))
			        						{
				        						NormalizedForm_plus1=component[0]+"."+(Integer.parseInt(component[3])+1)+""+component[2]+">"+component[4];
				        						NormalizedForm_minus1=component[0]+"."+(Integer.parseInt(component[3])-1)+""+component[2]+">"+component[4];
			        						}
				        					String wildtype[]=component[2].split(",");
				        					String mutant[]=component[4].split(",");
				        					String positions[]=component[3].split(",");
				        					
				        					if(wildtype.length == positions.length && wildtype.length>1) //Pair of wildtype&position
				        					{
				        						for (int i=0;i<wildtype.length;i++) //May have more than one pair
					        					{
					        						for (int j=0;j<mutant.length;j++)
						        					{
					        							NormalizedForm=NormalizedForm+"|"+component[0]+"."+positions[i]+wildtype[i]+">"+mutant[j];
						        					}
					        					}
				        					}
				        					else
				        					{
					        					for (int i=0;i<wildtype.length;i++) //May have more than one pair
					        					{
					        						for (int j=0;j<mutant.length;j++)
						        					{
					        							NormalizedForm=NormalizedForm+"|"+component[0]+"."+component[3]+wildtype[i]+">"+mutant[j];
						        					}
					        					}
					        					component[3]=component[3].replaceAll(",", "");
				        					}
				        					if(component[3].matches("[0-9]{5,}"))
				        					{
				        						component[0]="g";
				        						NormalizedForm=NormalizedForm+"|"+component[0]+"."+component[3]+component[2]+">"+component[4];
				        					}
				        					
				        					//protein
			        						{
				        						String tmp="";
				        						/*one -> three*/
				        						for(int len=0;len<component[2].length();len++)
				        						{
				        							if(one2three.containsKey(component[2].substring(len, len+1)))
				        							{
				        								if(tmp.equals(""))
				        								{
				        									tmp = one2three.get(component[2].substring(len, len+1));
				        								}
				        								else
				        								{
				        									tmp = tmp +","+ one2three.get(component[2].substring(len, len+1));
				        								}
				        							}
				        						}
				        						component[2]=tmp;
				        						tmp="";
				        						for(int len=0;len<component[4].length();len++)
				        						{
				        							if(one2three.containsKey(component[4].substring(len, len+1)))
				        							{
				        								if(tmp.equals(""))
				        								{
				        									tmp = one2three.get(component[4].substring(len, len+1));	
				        								}
				        								else
				        								{	     
				        									tmp = tmp +","+ one2three.get(component[4].substring(len, len+1));	
				        								}
				        							}
				        						}
				        						component[4]=tmp;
				        						
				        						if(component[2].equals(component[4]))
				        						{
				        							NormalizedForm_protein="p."+component[2]+component[3]+"=";
				        						}
				        						else
				        						{
				        							NormalizedForm_protein="p."+component[2]+component[3]+component[4];
					        					}
					        					wildtype=component[2].split(",");
					        					mutant=component[4].split(",");
					        					positions=component[3].split(",");
					        					
					        					if(wildtype.length == positions.length && wildtype.length>1) //Pair of wildtype&position
					        					{
					        						for (int i=0;i<wildtype.length;i++) //May have more than one pair
						        					{
						        						for (int j=0;j<mutant.length;j++)
							        					{
						        							NormalizedForm_protein=NormalizedForm_protein+"|"+"p."+wildtype[i]+positions[i]+mutant[j];
							        					}
						        					}
					        					}
					        					else
					        					{
					        						for (int i=0;i<wildtype.length;i++) //May have more than one pair
						        					{
						        						for (int j=0;j<mutant.length;j++)
							        					{
						        							NormalizedForm_protein=NormalizedForm_protein+"|"+"p."+wildtype[i]+component[3]+mutant[j];
							        					}
						        					}
					        					}
				        					}
		    	        				}
			        				}
		        				}
		        				else if(component[1].equals("DEL"))
		        				{
		        					if(component.length>3)
		        					{
		        						if(component[0].equals("p"))
				        				{
			        						String tmp="";
			        						for(int len=0;len<component[3].length();len++)
			        						{
			        							tmp = tmp + one2three.get(component[3].substring(len, len+1));
			        						}
			        						component[3]=tmp;
			        					}
			        					NormalizedForm=component[0]+"."+component[2]+"del"+component[3];
		        					}
		        					NormalizedForm=NormalizedForm+"|"+component[0]+"."+component[2]+"del";
		        				}
		        				else if(component[1].equals("INS"))
		        				{
		        					if(component.length>3)
		        					{
		        						if(component[0].equals("p"))
				        				{
			        						String tmp="";
			        						for(int len=0;len<component[3].length();len++)
			        						{
			        							tmp = tmp + one2three.get(component[3].substring(len, len+1));
			        						}
			        						component[3]=tmp;
			        					}
			        					NormalizedForm=component[0]+"."+component[2]+"ins"+component[3];
		        					}
		        					NormalizedForm=NormalizedForm+"|"+component[0]+"."+component[2]+"ins";
		        				}
		        				else if(component[1].equals("INDEL"))
		        				{
		        					if(component.length>3)
		        					{
		        						//c.*2361_*2362delAAinsA
			        					//c.2153_2155delinsTCCTGGTTTA
			        					if(component[0].equals("p"))
				        				{
			        						String tmp="";
			        						for(int len=0;len<component[3].length();len++)
			        						{
			        							tmp = tmp + one2three.get(component[3].substring(len, len+1));
			        						}
			        						component[3]=tmp;
			        					}
			        					NormalizedForm=component[0]+"."+component[2]+"delins"+component[3];
		        					}
		        				}
		        				else if(component[1].equals("DUP"))
		        				{
		        					if(component.length>3)
		        					{
		        						if(component[0].equals("p"))
				        				{
			        						String tmp="";
			        						for(int len=0;len<component[3].length();len++)
			        						{
			        							tmp = tmp + one2three.get(component[3].substring(len, len+1));
			        						}
			        						component[3]=tmp;
			        					}
			        					NormalizedForm=component[0]+"."+component[2]+"dup"+component[3];
		        					}
		        					NormalizedForm=NormalizedForm+"|"+component[0]+"."+component[2]+"dup";
		        				}
		        				else if(component[1].equals("FS"))
		        				{
		        					if(component[0].equals("p"))
			        				{
		        						String tmp="";
		        						for(int len=0;len<component[2].length();len++)
		        						{
		        							tmp = tmp + one2three.get(component[2].substring(len, len+1));
		        						}
		        						component[2]=tmp;
		        						
		        						tmp="";
		        						if(component.length>=5)
		        						{
			        						for(int len=0;len<component[4].length();len++)
			        						{
			        							tmp = tmp + one2three.get(component[4].substring(len, len+1));
			        						}
			        						component[4]=tmp;
		        						}
		        					}
		        					
		        					if(component.length>=5)
		        					{
		        						NormalizedForm=component[0]+"."+component[2]+component[3]+component[4]+"fs";
		        					}
		        					else if(component.length==4)
		        					{
		        						NormalizedForm=component[0]+"."+component[2]+component[3]+"fs";
		        					}
		        				}
	        				}
	        				String NormalizedForms[]=NormalizedForm.split("\\|");
	        				HashMap<String,String> NormalizedForm_hash = new HashMap<String,String>();
	        				for(int n=0;n<NormalizedForms.length;n++)
	        				{
	        					NormalizedForm_hash.put(NormalizedForms[n], "");
	        				}
	        				NormalizedForm="";
	        				for(String NF : NormalizedForm_hash.keySet())
	        				{
	        					if(NormalizedForm.equals(""))
	        					{
	        						NormalizedForm=NF;
	        					}
	        					else
	        					{
	        						NormalizedForm=NormalizedForm+"|"+NF;
	        					}
	        				}
	        				
	        				if(MentionPatternMap2rs.containsKey(anno[0]+"\t"+anno[5])) //by pattern (RS# in the pattern)
	        				{
	        					outputSTR=outputSTR+line+";Pattern-RS#:"+MentionPatternMap2rs.get(anno[0]+"\t"+anno[5])+"\n";
	        				}
	        				else if(rs_foundinText_hash.containsKey(anno[0]+"\t"+NormalizedForm)) //recognized (RS# in the text)
	        				{
	        					outputSTR=outputSTR+line+";Recognized-RS#:"+rs_foundinText_hash.get(anno[0]+"\t"+NormalizedForm)+"\n";
	        				}
	        				else if(rs_foundbefore_hash.containsKey(anno[0]+"\t"+NormalizedForm)) //found by dictionary-lookup before
	        				{
	        					outputSTR=outputSTR+line+";Foundbefore-RS#:"+rs_foundbefore_hash.get(anno[0]+"\t"+NormalizedForm)+"\n";
	        					MentionPatternMap2rs_extend.put(anno[0]+"\t"+anno[5],rs_foundbefore_hash.get(anno[0]+"\t"+NormalizedForm));
	        				}
	        				else
	        				{
	        					boolean found=false;
		        				String id_tmp=anno[5];
	        					id_tmp=id_tmp.replaceAll("^\\+", "");
	        					for(String pmid_geneid : tmp_gene.keySet())
	        					{
	        						String pmid_geneid_split[] = pmid_geneid.split("\t");
	        						String pmid=pmid_geneid_split[0];
	        						String geneid=pmid_geneid_split[1];
	        						if(pmid.equals(anno[0]))
	        						{
		        						if(Mutation_RS_Geneid_hash.containsKey(id_tmp+"\t"+geneid))
			        					{
		        							outputSTR=outputSTR+line+";Mapping-RS#:"+Mutation_RS_Geneid_hash.get(anno[5]+"\t"+geneid)+"\n"; // by history pairs
		        							MentionPatternMap2rs_extend.put(anno[0]+"\t"+anno[5],Mutation_RS_Geneid_hash.get(anno[5]+"\t"+geneid));
		        							found=true;
		        							break;
		        						}
	        						}
	        					}
	        					if(found == false)
	        					{
			        				/*
			        				 * Search Database - Mutation
			        				 */
		        					if(tmp_mutation.containsKey(anno[5]))
			        				{
			        					annotations_mutation.put(anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+anno[5],tmp_mutation.get(anno[5]));
			        				}
			        				else
			        				{
				        				try {
				        					String DB="var2rs_c";
				        					String Table="var2rs_c";
				        					String DBX="var2rs_Xc";
				        					if(component[0].equals("p"))
					        				{	
				        						DB="var2rs_p";
				        						Table="var2rs_p";
				        						DBX="var2rs_Xp";
					        				}
				        					else if(component[0].equals("c") || component[0].equals("r"))
					        				{	
				        						DB="var2rs_c";
				        						Table="var2rs_c";
				        						DBX="var2rs_Xc";
					        				}
				        					else if(component[0].equals("m"))
					        				{	
				        						DB="var2rs_m";
				        						Table="var2rs_m";
				        						DBX="var2rs_Xm";
					        				}
				        					else if(component[0].equals("n"))
					        				{	
				        						DB="var2rs_n";
				        						Table="var2rs_n";
				        						DBX="var2rs_Xn";
					        				}
				        					else if(component[0].equals("g"))
					        				{	
				        						DB="var2rs_g";
				        						Table="var2rs_g";
				        						DBX="var2rs_Xg";
					        				}
				        					
				        					String  rsNumber="";
				        					
				        					//N prefix
				        					{
				        						c = DriverManager.getConnection("jdbc:sqlite:Database/"+DB+".db");
					        					stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm.split("\\|");
					        					String SQL="SELECT rs FROM "+Table+" WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
				        					}
			        			
				        					//X prefix
				        					{
					        					c = DriverManager.getConnection("jdbc:sqlite:Database/"+DBX+".db");
					        					stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm.split("\\|");
					        					String SQL="SELECT rs FROM "+Table+" WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
				        					}
			        						
			        						//NormalizedForm_reverse
				        					if(!NormalizedForm_reverse.equals(""))
				        					{
					        					c = DriverManager.getConnection("jdbc:sqlite:Database/"+DB+".db");
					        					stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm_reverse.split("\\|");
					        					String SQL="SELECT rs FROM "+Table+" WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
				        					}
				        					
			        						//NormalizedForm - c.-844A>G --> n.-844A>G
				        					/*
			        						if(NormalizedForm.matches(".*[\\-\\_\\+\\*].*") && NormalizedForm.matches("[crmgn]\\..*"))
			        						{
			        							String NormalizedForm_rev = NormalizedForm.replaceAll("[crmg]\\.","n\\.");
			        							c = DriverManager.getConnection("jdbc:sqlite:Database/var2rs_n.db");
					        					stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm_rev.split("\\|");
					        					String SQL="SELECT rs FROM var2rs_n WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
			        						}
			        						*/
			        						//NormalizedForm_reverse - c.-844A>G --> n.-844A>G
			        						/*
			        						if((!NormalizedForm_reverse.equals("")) && NormalizedForm_reverse.matches(".*[\\-\\_\\+\\*].*") && NormalizedForm_reverse.matches("[crmgn]\\..*"))
			        						{
			        							String NormalizedForm_rev = NormalizedForm_reverse.replaceAll("[crmg]\\.","n\\.");
			        							c = DriverManager.getConnection("jdbc:sqlite:Database/var2rs_n.db");
					        					stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm_rev.split("\\|");
					        					String SQL="SELECT rs FROM var2rs_n WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
			        						}
			        						*/
			        						
			        						//NormalizedForm - proteinchange
			        						if(NormalizedForm.matches("p\\..*"))
			        						{
			        							c = DriverManager.getConnection("jdbc:sqlite:Database/var2rs_clinvar_proteinchange.db");
					        					stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm.split("\\|");
					        					String SQL="SELECT rs FROM var2rs_clinvar_proteinchange WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
			        						}
			        						
			        						//NormalizedForm - reverse_proteinchange
			        						if(NormalizedForm.matches("p\\..*"))
			        						{
			        							c = DriverManager.getConnection("jdbc:sqlite:Database/var2rs_clinvar_proteinchange.db");
				        						stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm_reverse.split("\\|");
					        					String SQL="SELECT rs FROM var2rs_clinvar_proteinchange WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
			        						}
			        						
			        						//NormalizedForm_plus1
			        						/*
			        						if(!NormalizedForm_plus1.equals(""))
			        						{
				        						c = DriverManager.getConnection("jdbc:sqlite:Database/"+DB+".db");
					        					stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm_plus1.split("\\|");
					        					String SQL="SELECT rs FROM "+Table+" WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
			        						}
			        						*/
			        						
			        						//NormalizedForm_minus1
			        						/*
			        						if(!NormalizedForm_minus1.equals(""))
			        						{
				        						c = DriverManager.getConnection("jdbc:sqlite:Database/"+DB+".db");
					        					stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm_minus1.split("\\|");
					        					String SQL="SELECT rs FROM "+Table+" WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
				        					}
											*/
			        						
			        						//if DNA cannot found, try protein
			        						if(!NormalizedForm_protein.equals(""))
			        						{
			        							DB="var2rs_p";
				        						Table="var2rs_p";
				        						DBX="var2rs_Xp";
				        						c = DriverManager.getConnection("jdbc:sqlite:Database/"+DB+".db");
					        					stmt = c.createStatement();
					        					String NormalizedForm_arr[]=NormalizedForm_protein.split("\\|");
					        					String SQL="SELECT rs FROM "+Table+" WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					ResultSet rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
					        					
				        						//X prefix
					        					c = DriverManager.getConnection("jdbc:sqlite:Database/"+DBX+".db");
					        					stmt = c.createStatement();
					        					NormalizedForm_arr=NormalizedForm_protein.split("\\|");
					        					SQL="SELECT rs FROM "+Table+" WHERE ";
							        			for(int nfa=0;nfa<NormalizedForm_arr.length;nfa++)
					        					{
							        				SQL=SQL+"var='"+NormalizedForm_arr[nfa]+"' or ";
					        					}
					        					SQL=SQL.replaceAll(" or $", "");
					        					rs = stmt.executeQuery(SQL);
					        					while ( rs.next() ) 
					        					{
					        						if(rsNumber.equals(""))
					        						{
					        							rsNumber = rs.getString("rs");
					        						}
					        						else
					        						{
					        							rsNumber = rsNumber+"|"+rs.getString("rs");
					        						}
					        					}
				        						stmt.close();
				        						c.close();
				        						rsNumber=rsNumber+"|end";
			        						}
			        						
			        						annotations_mutation.put(anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+anno[5],rsNumber);
			        				        tmp_mutation.put(anno[5],rsNumber);
			        					} 
			        					catch ( SQLException e ) 
			        					{
			        						System.err.println( e.getClass().getName() + ": " + e.getMessage() );
			        					}
			        				}
			        				
			        				/*
			        				 * normalization by snp mention in the text
			        				 */
			        				if(annotations_mutation.containsKey(anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+anno[5]))
			        				{
				        				String rsNumbers_arr[] = annotations_mutation.get(anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+anno[5]).split("\\|");
				        				for(int ra=0;ra<rsNumbers_arr.length;ra++)
				        				{
				        					if(RSnumberList.containsKey(anno[0]+"\t"+rsNumbers_arr[ra]))
					        				{
				        						rs_foundinText_hash.put(anno[0]+"\t"+NormalizedForm, rsNumbers_arr[ra]);
				        						outputSTR=outputSTR+line+";SNPinText-RS#:"+rsNumbers_arr[ra]+"\n";
				        						MentionPatternMap2rs_extend.put(anno[0]+"\t"+anno[5],rsNumbers_arr[ra]);
					        					ra = rsNumbers_arr.length; //last
					        					found=true;
					        				}
				        				}
				        			}
			        				if(found == false)
			        				{
				        				/*
				        				 * normalization by gene
				        				 */
			        					HashMap<String,Integer> geneRS2distance_hash= new HashMap<String,Integer>();
			        					HashMap<Integer,String>  sentencelocation_gene_hash= new HashMap<Integer,String>();
			        					int sentencelocation_rs=0;
			        					int sentencelocation_gene=0;
		        						for(String gene : annotations_gene.keySet())
				        				{
				        					String gene_info[] = gene.split("\t"); //pmid,start,last,mention,geneid
				        					
				        					if(anno[0].equals(gene_info[0]) && gene_info.length>=5)
				        					{
				        						/*
				        						 * detect the sentence boundaries for gene and RS#
				        						 */
				        						ArrayList<Integer> SentenceOffsets = PMID2SentenceOffsets.get(anno[0]);
				        						for(int si=0;si<SentenceOffsets.size();si++) // find sentence location for gene mention
				        						{
				        							if(Integer.parseInt(gene_info[1])<=SentenceOffsets.get(si))
				        							{
				        								sentencelocation_gene_hash.put(si-1,gene_info[4]); //gene_info[1] : gene identifier
				        								sentencelocation_gene=si-1;
				        								break;
				        							}
				        						}
				        						for(int si=0;si<SentenceOffsets.size();si++) // find sentence location for rs mention 
				        						{
				        							if(Integer.parseInt(anno[1])<=SentenceOffsets.get(si))
				        							{
				        								sentencelocation_rs=si-1;
				        								break;
				        							}
				        						}
				        										        						
				        						if(Integer.parseInt(gene_info[2])<=Integer.parseInt(anno[1])) // gene --- mutation
				        						{
				        							int distance = Integer.parseInt(anno[1])-Integer.parseInt(gene_info[2]);
				        							if(sentencelocation_gene != sentencelocation_rs){distance=distance+1000;}
					        						if(!geneRS2distance_hash.containsKey(annotations_gene.get(gene)))
				        							{
					        							geneRS2distance_hash.put(annotations_gene.get(gene),distance);
				        							}
				        							else if(distance<geneRS2distance_hash.get(annotations_gene.get(gene)))
				        							{
				        								geneRS2distance_hash.put(annotations_gene.get(gene),distance);
				        							}
				        						}
				        						else if(Integer.parseInt(gene_info[1])>=Integer.parseInt(anno[2])) // mutation --- gene
				        						{
				        							int distance = Integer.parseInt(gene_info[1])-Integer.parseInt(anno[2]);
				        							if(sentencelocation_gene != sentencelocation_rs){distance=distance+1000;}
				        							if(!geneRS2distance_hash.containsKey(annotations_gene.get(gene)))
				        							{
				        								geneRS2distance_hash.put(annotations_gene.get(gene),distance);
				        							}
				        							else if(distance<geneRS2distance_hash.get(annotations_gene.get(gene)))
				        							{
				        								geneRS2distance_hash.put(annotations_gene.get(gene),distance);
				        							}
				        						}
				        						else // mutation & gene may overlap
				        						{
				        							geneRS2distance_hash.put(annotations_gene.get(gene),0);
				        						}
				        					}
				        				}
				        				
				        				ArrayList <String> geneRS_ranked = new ArrayList <String>();
				        				while(!geneRS2distance_hash.isEmpty())
				        				{
				        					int closet_distance=10000;
					        				String closet_geneRS="";
					        				for (String geneRS : geneRS2distance_hash.keySet())
					        				{
					        					if(geneRS2distance_hash.get(geneRS)<closet_distance)
					        					{
					        						closet_distance=geneRS2distance_hash.get(geneRS);
					        						closet_geneRS=geneRS;
					        					}
					        				}
					        				if(closet_geneRS.equals(""))
					        				{	
						        				break;
					        				}
					        				geneRS_ranked.add(closet_geneRS);
					        				geneRS2distance_hash.remove(closet_geneRS);
				        				}
				        				
				        				found = false;
				        				
				        				int number_geneRS=geneRS_ranked.size();
				        				if(number_geneRS>0 && sentencelocation_gene_hash.containsKey(sentencelocation_rs)) // if the closet gene and RS# in the same sentence, only look at the closet gene.
				        				{
				        					number_geneRS=1;
				        				}
				        				
				        				for(int rsg=0;rsg<number_geneRS;rsg++)
			        					{
				        					String target_gene_rs=geneRS_ranked.get(rsg);
					        				String geneRSs[]=target_gene_rs.split("\\|");
					        				if(annotations_mutation.containsKey(anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+anno[5]))
					        				{
					        					String mutationRSs[]=annotations_mutation.get(anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+anno[5]).split("\\|");
						        				String found_rs = "";
						        				
						        				for(int m=0;m<mutationRSs.length;m++)
						        				{
						        					for(int g=0;g<geneRSs.length;g++)
						        					{
						        						if(found==true && mutationRSs[m].equals("end"))
						        						{
						        							m=mutationRSs.length;
						        							g=geneRSs.length;
						        						}
						        						else if(geneRSs[g].equals(mutationRSs[m]) && !geneRSs[g].equals(""))
						        						{
						        							found=true;
						        							if(found_rs.equals(""))
						        							{
						        								found_rs = geneRSs[g];
						        							}
						        							else
						        							{
						        								found_rs = found_rs+"|"+ geneRSs[g];
						        							}
						        						}
							        				}
						        				}
						        				if(!found_rs.equals(""))
						        				{
						        					rs_foundbefore_hash.put(anno[0]+"\t"+NormalizedForm, found_rs);
						        				}
						        				if(found == true)
						        				{
						        					outputSTR=outputSTR+line+";RS#:"+found_rs+"\n";
						        					MentionPatternMap2rs_extend.put(anno[0]+"\t"+anno[5],found_rs);
						        					rsg=geneRS_ranked.size();
						        				}
					        				}
				        				}
				        				if(found == false)
				        				{
				        					if(MentionPatternMap2rs_extend.containsKey(anno[0]+"\t"+anno[5])) //by pattern [PMID:24073655] c.169G>C (p.Gly57Arg) -- map to rs764352037
					        				{
					        					outputSTR=outputSTR+line+";Extended-RS#:"+MentionPatternMap2rs_extend.get(anno[0]+"\t"+anno[5])+"\n";
					        				}
				        					else
				        					{
				        						outputSTR=outputSTR+line+"\n";
				        					}
				        				}
			        				}
	        					}
	        				}
	        				if(MentionPatternMap2rs_extend.containsKey(anno[0]+"\t"+anno[5]))
	        				{
	        					if(MentionPatternMap.containsKey(anno[0]+"\t"+anno[5]))
	        					{
	        						MentionPatternMap2rs_extend.put(anno[0]+"\t"+MentionPatternMap.get(anno[0]+"\t"+anno[5]), MentionPatternMap2rs_extend.get(anno[0]+"\t"+anno[5]));
	        						String tmp_anno5 = MentionPatternMap.get(anno[0]+"\t"+anno[5]).replaceAll("([^A-Za-z0-9])","\\\\$1");
	        						outputSTR=outputSTR.replaceAll(tmp_anno5+"\n", anno[5]+";Extended-RS#:"+MentionPatternMap2rs_extend.get(anno[0]+"\t"+anno[5])+"\n");
	        					}
	        					else
	        					{
	        						MentionPatternMap2rs_extend.put(anno[0]+"\t"+anno[5], MentionPatternMap2rs_extend.get(anno[0]+"\t"+anno[5]));
	        						String tmp_anno5 = anno[5].replaceAll("([^A-Za-z0-9])","\\\\$1");
	        						outputSTR=outputSTR.replaceAll(tmp_anno5+"\n", anno[5]+";Extended-RS#:"+MentionPatternMap2rs_extend.get(anno[0]+"\t"+anno[5])+"\n");
	        					}	
	        				}
	        			}
	        			else
		        		{
	        				outputSTR=outputSTR+line+"\n";
		        		}
	        		}
	        		else
	        		{
	        			outputSTR=outputSTR+line+"\n";
	        		}
	        	}
				else if (line.equals(""))
				{
					if(DisplayRSnumOnly.equals("True"))
					{
						outputSTR=outputSTR.replaceAll(";(|Pattern-|Extended-|Recognized-|Mapping-|SNPinText-|Foundbefore-)RS#:", ";RS#:");
					}
					outputfile.write(outputSTR+"\n");
					outputSTR="";
				}
				else
				{
					outputSTR=outputSTR+line+"\n";
				}
			}
			inputfile.close();
			outputfile.close();
		}
		catch(IOException e1){ System.out.println("[normalization]: "+e1+" Input file is not exist.");}
	}	
}
