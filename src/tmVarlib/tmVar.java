package tmVarlib;
//
// tmVar - Java version
//

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamException;

import org.tartarus.snowball.SnowballStemmer;
import org.tartarus.snowball.ext.englishStemmer;

import edu.stanford.nlp.tagger.maxent.MaxentTagger;

public class tmVar
{
	
	public static MaxentTagger tagger = new MaxentTagger();
	public static SnowballStemmer stemmer = new englishStemmer();
	public static ArrayList<String> RegEx_DNAMutation_STR=new ArrayList<String>(); 
	public static ArrayList<String> RegEx_ProteinMutation_STR=new ArrayList<String>(); 
	public static ArrayList<String> RegEx_SNP_STR=new ArrayList<String>();
	public static ArrayList<String> PAM_lowerScorePair = new ArrayList<String>();
	public static Pattern Pattern_Component_1;
	public static Pattern Pattern_Component_2;
	public static Pattern Pattern_Component_3;
	public static Pattern Pattern_Component_4;
	public static Pattern Pattern_Component_5;
	public static Pattern Pattern_Component_6;
	public static boolean GeneMention = false; // will be turn to true if "ExtractFeature" can find gene mention
	
	public static void main(String [] args) throws IOException, InterruptedException, XMLStreamException 
	{
		/*
		 * Parameters
		 */
		String InputFolder="input";
		String OutputFolder="output";
		String TrainTest="Test"; //Train|Train_Mention|Test
		String DeleteTmp="True";
		String DisplayRSnumOnly="True";
		if(args.length<2)
		{
			System.out.println("\n$ java -Xmx5G -Xms5G -jar tmVar.jar [InputFolder] [OutputFolder]");
			System.out.println("[InputFolder] Default : input");
			System.out.println("[OutputFolder] Default : output\n\n");
		}
		else
		{
			InputFolder=args [0];
			OutputFolder=args [1];
			
			if(args.length>2 && args[2].matches("(Train|Train_Mention|Test|Test_FullText)"))
			{
				TrainTest=args [2];
				if(args[2].matches("(Train|Train_Mention)"))
				{
					DeleteTmp="False";
				}
			}
			if(args.length>3 && args[3].toLowerCase().matches("(true|false)"))
			{
				DeleteTmp=args [3];
			}
			if(args.length>4 && args[4].toLowerCase().matches("(true|false)"))
			{
				DisplayRSnumOnly=args [4];
			}
		}
		
		double startTime,endTime,totTime;
		startTime = System.currentTimeMillis();//start time
		BioCConverter BC= new BioCConverter();
		
		/**
		 * Import models
		 */
		{
			/*
			 * POSTagging: loading model
			 */
			tagger = new MaxentTagger("lib/taggers/english-left3words-distsim.tagger");
			
			/*
			 * Stemming : using Snowball
			 */
			stemmer = new englishStemmer();
			
			/*
			 * PAM 140 : <=-6 pairs 
			 */
			BufferedReader PAM = new BufferedReader(new InputStreamReader(new FileInputStream("lib/PAM140-6.txt"), "UTF-8"));
			String line="";
			while ((line = PAM.readLine()) != null)  
			{
				String nt[]=line.split("\t");
				PAM_lowerScorePair.add(nt[0]+"\t"+nt[1]);
				PAM_lowerScorePair.add(nt[1]+"\t"+nt[0]);
			}
			PAM.close();
			
			/*
			 * HGVs nomenclature lookup - RegEx : DNAMutation
			 */
			BufferedReader RegEx_DNAMutation = new BufferedReader(new InputStreamReader(new FileInputStream("lib/RegEx/DNAMutation.RegEx.txt"), "UTF-8"));
			line="";
			while ((line = RegEx_DNAMutation.readLine()) != null)  
			{
				RegEx_DNAMutation_STR.add(line);
			}
			RegEx_DNAMutation.close();
			
			/*
			 * HGVs nomenclature lookup - RegEx : ProteinMutation
			 */
			BufferedReader RegEx_ProteinMutation = new BufferedReader(new InputStreamReader(new FileInputStream("lib/RegEx/ProteinMutation.RegEx.txt"), "UTF-8"));
			line="";
			while ((line = RegEx_ProteinMutation.readLine()) != null)  
			{
				RegEx_ProteinMutation_STR.add(line);
			}
			RegEx_ProteinMutation.close();
			
			/*
			 * HGVs nomenclature lookup - RegEx : SNP
			 */
			BufferedReader RegEx_SNP = new BufferedReader(new InputStreamReader(new FileInputStream("lib/RegEx/SNP.RegEx.txt"), "UTF-8"));
			line="";
			while ((line = RegEx_SNP.readLine()) != null)  
			{
				RegEx_SNP_STR.add(line);
			}
			RegEx_SNP.close();
			
			//RegEx of component recognition
			Pattern_Component_1 = Pattern.compile("^([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)\\|(fs[^|]*)\\|([^|]*)$");
			Pattern_Component_2 = Pattern.compile("^([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)\\|(fs[^|]*)$");
			Pattern_Component_3 = Pattern.compile("^([^|]*)\\|([^|]*(ins|del|Del|dup|-)[^|]*)\\|([^|]*)\\|([^|]*)$");
			Pattern_Component_4 = Pattern.compile("^([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)$");
			Pattern_Component_5 = Pattern.compile("^([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)\\|([^|]*)$");
			Pattern_Component_6 = Pattern.compile("^((\\[rs\\]|rs|RS|Rs|reference SNP no[.] )[0-9]+)$");
		}
		
		File folder = new File(InputFolder);
		File[] listOfFiles = folder.listFiles();
		for (int i = 0; i < listOfFiles.length; i++)
		{
			if (listOfFiles[i].isFile()) 
			{
				String InputFile = listOfFiles[i].getName();
				
				File f = new File(OutputFolder+"/"+InputFile+".PubTator");
				if(f.exists() && !f.isDirectory()) 
				{ 
					System.out.println(InputFolder+"/"+InputFile+" - Done. (The output file exists in output folder)");
				}
				else
				{
					/*
					 * Mention recognition by CRF++
					 */
					if(TrainTest.equals("Test") || TrainTest.equals("Test_FullText") || TrainTest.equals("Train"))
					{
						/*
						 * Format Check 
						 */
						String Format = "";
						String checkR=BC.BioCFormatCheck(InputFolder+"/"+InputFile);
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
							System.out.println(checkR);
							System.exit(0);
						}
						
						System.out.print(InputFolder+"/"+InputFile+" - ("+Format+" format) : Processing ... \r");
						 
						/*
						 * Pre-processing
						 */
						MentionRecognition MR= new MentionRecognition();
						if(Format.equals("BioC"))
						{
							BC.BioC2PubTator(InputFolder+"/"+InputFile,"tmp/"+InputFile);
							MR.FeatureExtraction("tmp/"+InputFile,"tmp/"+InputFile+".data","tmp/"+InputFile+".location",TrainTest);
						}
						else if(Format.equals("PubTator"))
						{
							MR.FeatureExtraction(InputFolder+"/"+InputFile,"tmp/"+InputFile+".data","tmp/"+InputFile+".location",TrainTest);
						}
						if(TrainTest.equals("Test") || TrainTest.equals("Test_FullText"))
						{
							MR.CRF_test("tmp/"+InputFile+".data","tmp/"+InputFile+".output",TrainTest);
						}
						
						/*
						 * CRF++ output --> PubTator
						 */
						PostProcessing PP = new PostProcessing();
						{
							if(Format.equals("BioC"))
							{
								PP.toME("tmp/"+InputFile,"tmp/"+InputFile+".output","tmp/"+InputFile+".location","tmp/"+InputFile+".ME");
								PP.toPostME("tmp/"+InputFile+".ME","tmp/"+InputFile+".PostME");
								PP.toPostMEData("tmp/"+InputFile,"tmp/"+InputFile+".PostME","tmp/"+InputFile+".PostME.ml","tmp/"+InputFile+".PostME.data",TrainTest);
							}
							else if(Format.equals("PubTator"))
							{
								PP.toME(InputFolder+"/"+InputFile,"tmp/"+InputFile+".output","tmp/"+InputFile+".location","tmp/"+InputFile+".ME");
								PP.toPostME("tmp/"+InputFile+".ME","tmp/"+InputFile+".PostME");
								PP.toPostMEData(InputFolder+"/"+InputFile,"tmp/"+InputFile+".PostME","tmp/"+InputFile+".PostME.ml","tmp/"+InputFile+".PostME.data",TrainTest);
							}
							if(TrainTest.equals("Test") || TrainTest.equals("Test_FullText"))
							{
								PP.toPostMEoutput("tmp/"+InputFile+".PostME.data","tmp/"+InputFile+".PostME.output");
							}
							
							else if(TrainTest.equals("Train"))
							{
								PP.toPostMEModel("tmp/"+InputFile+".PostME.data");
							}
							
							
							/*
							 * Post-processing
							 */
							if(TrainTest.equals("Test") || TrainTest.equals("Test_FullText"))
							{
								//GeneMention = true;
								if(GeneMention == true) // MentionRecognition detect Gene mentions
								{
									PP.output2PubTator("tmp/"+InputFile+".PostME.ml","tmp/"+InputFile+".PostME.output","tmp/"+InputFile+".PostME","tmp/"+InputFile+".PubTator");
									
									if(Format.equals("BioC"))
									{
										PP.Normalization(DisplayRSnumOnly,"tmp/"+InputFile,"tmp/"+InputFile+".PubTator",OutputFolder+"/"+InputFile+".PubTator");
										
									}
									else if(Format.equals("PubTator"))
									{
										PP.Normalization(DisplayRSnumOnly,InputFolder+"/"+InputFile,"tmp/"+InputFile+".PubTator",OutputFolder+"/"+InputFile+".PubTator");
									}
								}
								else
								{
									PP.output2PubTator("tmp/"+InputFile+".PostME.ml","tmp/"+InputFile+".PostME.output","tmp/"+InputFile+".PostME",OutputFolder+"/"+InputFile+".PubTator");
								}
								
								if(Format.equals("BioC"))
								{
									BC.PubTator2BioC_AppendAnnotation(OutputFolder+"/"+InputFile+".PubTator",InputFolder+"/"+InputFile,OutputFolder+"/"+InputFile+".BioC.XML");
								}			
							}
						}
						
						/*
						 * Time stamp - last
						 */
						endTime = System.currentTimeMillis();//ending time
						totTime = endTime - startTime;
						System.out.println(InputFolder+"/"+InputFile+" - ("+Format+" format) : Processing Time:"+totTime/1000+"sec");
						
						/*
						 * remove tmp files
						 */
						if(DeleteTmp.toLowerCase().equals("true"))
						{
							String path="tmp"; 
					        File file = new File(path);
					        File[] files = file.listFiles(); 
					        for (File ftmp:files) 
					        {
					        	if (ftmp.isFile() && ftmp.exists()) 
					            {
					        		if(ftmp.toString().matches("tmp."+InputFile+".*"))
						        	{
					        			ftmp.delete();
						        	}
					        	}
					        }
						}
					}
					else if(TrainTest.equals("Train_Mention"))
					{
						System.out.print(InputFolder+"/"+InputFile+" - Processing ... \r");
						 
						PostProcessing PP = new PostProcessing();
						PP.toPostMEData(InputFolder+"/"+InputFile,"tmp/"+InputFile+".PostME","tmp/"+InputFile+".PostME.ml","tmp/"+InputFile+".PostME.data","Train");
						
						/*
						 * Time stamp - last
						 */
						endTime = System.currentTimeMillis();//ending time
						totTime = endTime - startTime;
						System.out.println(InputFolder+"/"+InputFile+" - Processing Time:"+totTime/1000+"sec");
					}
				}
			}
		}
	}
}
