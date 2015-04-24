
/**
 * Title:        Project of visualizing codon frequencies
 * Description:
 * Copyright:    Copyright (c) 2002
 * Company:      IHES, France
 * @author Andrey Zinovyev, Alessandra Carbone
 * @version 1.0
 * fisa comment: 1.8.2 (exactly) biojava version needed
 */

import java.io.*;
import java.util.*;
import java.lang.Integer;
import java.lang.System;


import org.biojava.bio.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.io.*;


public class CAIJava{

      static int taskModifier = 0;
      static String addFeature = null;
      static boolean calcGCContent = false;
      static boolean fastaFormat = false;
      static int calcCodonAnalysis = 0;
      static StringBuffer AdditionalInfoFile = new StringBuffer("");
      static float GTFAverLen[] = new float[64];
      static Vector TranslationalTable;
      static String externalWValues = null;
      static boolean checkForStability = false;
      static int classNumber = -1;


  public static void main(String [] args) {
	  String[] gbs = {"NC_003028"};
	 

		  args = "C:/Users/weidewind/workspace/asymmetry/NC_004347.gb -f C:/Users/weidewind/Documents/Asymmetry/2015/NC_004347_gcai_out.dat -t product -k 3 -i 15 >burk_out.veo".split("\\s++");
		  for (String mygbfile: gbs){
		  args[0] = "C:/Users/weidewind/workspace/asymmetry/" + mygbfile + ".gb";
		  args[2] = "C:/Users/weidewind/Documents/Asymmetry/2015/DataSetsProcessed/CAI/" + mygbfile + "_gcai_out2.dat";
	  try {
      if(args.length == 0) {
        throw new Exception("Use: CalcCodonFreq GenBankFile");
      }

      // test
      /*float tt[] = new float[11];
      tt[0]=(float)0.5;
      tt[1]=(float)0.1;
      tt[2]=(float)0.1;
      tt[3]=(float)0.3;
      tt[4]=(float)0.2;
      tt[5]=(float)0.4;
      tt[6]=(float)0.7;
      tt[7]=(float)0.9;
      tt[8]=(float)0.65;
      tt[9]=(float)0.05;
      tt[10]=(float)0.08;
      int nn[] = SortCais(tt);
      for (int i = 0; i < nn.length; i++) System.out.print(nn[i]+" ");
      System.out.println();
      for (int i = 0; i < nn.length; i++) System.out.print(tt[nn[i]]+" ");
      System.out.println();*/
      /*StringBuffer s = new StringBuffer();
      for (int i = 0; i < 64; i++) {
      s.append(Utils.TripletName(i));
      }
      calcAminoacidFreq(s.toString(),1,s.length());*/


      Vector Acids = new Vector();
      for (int i = 0; i < 64; i++){
      SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
      SymbolList Amin = RNATools.translate(RNATools.transcribe(sDNA));
      String sss = new String(Amin.seqString());
      if ((!Acids.contains(sss))&&(!sss.equals("*")))
          Acids.addElement(sss);
      }

      TranslationalTable = new Vector();
      for(int i=0;i<64;i++){
      SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
      SymbolList Amin = RNATools.translate(RNATools.transcribe(sDNA));
       if (!Amin.seqString().equals("*")) {
          int n = Acids.indexOf(Amin.seqString());
          TranslationalTable.addElement(new Integer(n));
          }
      else TranslationalTable.addElement(new Integer(-1));
      }


      int CodonsOrder[] = orderCodonUsageByAminoacid(Acids);
      Vector CodonUsages = new Vector();

      int windowSize = 120;
      int evPosition = 21, sAnnotbeg = -1, sAnnotend = -1;
      int numberOfIterations = 0;
      int maximumGeneLength = -1;
      int BayezRatio = 0;
      String fileName = "";

      for(int i=1;i<=args.length;i++)
      if (i<args.length)
        {
         if (args[i].length()>1)
         if(args[i].equals("-ew"))
             {
             externalWValues = args[i+1];
             continue;
             }
         if ((args[i].length()>2)&&(args[i].substring(0,3).equals("-cs"))) {checkForStability = true;
         if (args[i].length()>3) classNumber = Integer.parseInt(args[i].substring(3));;
         continue; }
         if(args[i].equals("-f")) fileName = new String(args[i+1]);
         if(args[i].equals("-i")) numberOfIterations = Integer.parseInt(args[i+1]);
         if(args[i].equals("-k")) taskModifier = Integer.parseInt(args[i+1]);
         if(args[i].equals("-a")) AdditionalInfoFile = new StringBuffer(args[i+1]);
         if(args[i].equals("-m")) maximumGeneLength = Integer.parseInt(args[i+1]);
         if(args[i].equals("-c")) calcCodonAnalysis = Integer.parseInt(args[i+1]);
         if(args[i].equals("-t")) addFeature = new String(args[i+1]);;
         if(args[i].equals("-s")) fastaFormat = true;
         if(args[i].equals("-g")) calcGCContent = true;
         }


//  Prepare gene file

  HashMap GeneAddInfo = new HashMap();
  String StringDelimiters = new String("\t");
  int GeneAddInfoCount=0;
  System.out.println("AdditionalInfoFile = \""+AdditionalInfoFile+"\"");
  if (!AdditionalInfoFile.toString().equals("")){
  System.out.println("Trying to read...");
  LineNumberReader lr = new LineNumberReader(new FileReader(AdditionalInfoFile.toString()));
  String s;
  s = lr.readLine();
  lr.close();
  lr = new LineNumberReader(new FileReader(AdditionalInfoFile.toString()));
  StringTokenizer st = new StringTokenizer(s,StringDelimiters);
  int j=0;
  while ((st.hasMoreTokens())) {st.nextToken(); j++; }
  GeneAddInfoCount = j;


  //lr.reset();
  while ( ((s=lr.readLine()) != null))
     {
     st = new StringTokenizer(s,StringDelimiters);
     //System.out.println("AddInfofile: "+s.substring(1,10)+"\t"+st.countTokens());
     String[] sInfo = new String[j];
     int k = 0;
     while ((st.hasMoreTokens())) {
        String sss = st.nextToken();
        if((k==1)&&(sss.length()==1)) sss="0"+sss;
        if(k<GeneAddInfoCount)
        sInfo[k] = sss.trim().replace(' ','_');
        else
        System.out.println("AddInfofile: out of bound "+s);
        k++;
        }
     if(GeneAddInfo.containsKey(sInfo[0]))
     { String sm[] = (String[])GeneAddInfo.get(sInfo[0]);
       for(int i=1;i<sm.length;i++) sm[i]=sm[i]+";"+sInfo[i]; }
     else
     GeneAddInfo.put(sInfo[0],sInfo);
     }
  }

  //for(Iterator i=GeneAddInfo.keySet().iterator();i.hasNext();){
  //String sm[] = (String[])GeneAddInfo.get((String)(i.next()));
  //for (int ii = 0; ii < sm.length; ii++) System.out.print(sm[ii]+" ");
  //System.out.println();
  //}
//


        // read external WValues
        float extWValues[] = new float[64];
        if(externalWValues!=null){
        LineNumberReader lr = new LineNumberReader(new FileReader(externalWValues));
        int k = 0;
         while(k<64){
         String s = lr.readLine();
         StringTokenizer st = new StringTokenizer(s,"\t ,");
         while ((st.hasMoreTokens()))
            {String ss = st.nextToken();
             String tn = Utils.TripletName(k,"t","c","a","g");
             int i =Utils.TripletNumber(tn);
             extWValues[i] = Float.parseFloat(ss);
             k++;
            }
         }
        }
        Utils.printDistribution(extWValues,"t","c","a","g");

// Random class file
int randomclassnumber = 0;
String randomclass[][] = new String[1000][4];
System.out.println("Looking for randomclass.txt...");
File ff = new File("randomclass.txt");
if(ff.exists()){
LineNumberReader lr = new LineNumberReader(new FileReader("randomclass.txt"));
String sr = null;
while((sr=lr.readLine())!=null){
  StringTokenizer st = new StringTokenizer(sr,"\t");
  randomclass[randomclassnumber][0] = st.nextToken();
  randomclass[randomclassnumber][1] = st.nextToken();
  randomclass[randomclassnumber][2] = st.nextToken();
  randomclass[randomclassnumber][3] = st.nextToken();
  randomclassnumber++;
}

}



HashMap GeneNoteTable = new HashMap();

      PrintStream fout;
      if (fileName.equals(""))
      fout = System.out;
      else
      fout = new PrintStream(new FileOutputStream(fileName));

      File GenBankFile = new File(args[0]);
      System.out.println("Loading sequence...");
      BufferedReader eReader = new BufferedReader(
        new InputStreamReader(new FileInputStream(GenBankFile)));
      SequenceIterator seqI = null;
      if(!fastaFormat)
      seqI = SeqIOTools.readGenbank(eReader);
      else
      seqI = SeqIOTools.readFastaDNA(eReader);
      System.out.println("Loaded...");

      int NumberOfPoints = 0;

      int sAnnBegM = sAnnotbeg;
      int sAnnEndM = sAnnotend;
      int TF[];
      float GTF[] = new float[64];

      // Calculating number of points and CAIW tables
      Vector ListOfGenes = new Vector();
      Vector CAIValues = new Vector();

      if(!fastaFormat)
      while(seqI.hasNext()) {
        System.out.println("Getting seq...");
        Sequence seq = seqI.nextSequence();
        System.out.println("Got "+seq.getName());
        int numCDSFeat = 0;
        for(Iterator i=seq.features(); i.hasNext(); )
          {
          Feature Ft = (Feature) i.next();
          if ((Ft.getType().equals("CDS"))&&(Ft.getSymbols().seqString().length()>10))
             {
             numCDSFeat++;
             try{
             String ss = new String(Ft.getSymbols().seqString());
             if(maximumGeneLength!=-1)
               if(ss.length()>maximumGeneLength)
                  continue;
             TF = Utils.calcTripletFreq(ss,1,ss.length());
    //         k=CalcCodonUsageForSeq(Ft.getSymbols(),fr0,fr1,fr2);
             //ss = new String(Ft.getSymbols().seqString());
             //if (ss.length()/3.0!=(int)(ss.length()/3.0))
             //    System.out.println("Mistake!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
             //System.out.println(Ft.getAnnotation().getProperty(new String("gene")).toString()+"; "+ss.length()+" symbols");
             //String check = new String(Ft.getAnnotation().getProperty("gene").toString());
             //GeneticCodes GCode = new GeneticCodes();
             //System.out.println(GCode.translate(GCode.transcribe(Ft.getSymbols())).seqString());
             //System.out.println(ss);
             for(int j=0;j<64;j++)
               GTF[j]+=TF[j];
             NumberOfPoints++;
             ListOfGenes.addElement(Ft);
             float CAIV[] = new float[numberOfIterations];
             CAIValues.addElement(CAIV);
              } catch (Exception e)
            {System.out.print("\n"+e+" : ");
             System.out.println(Ft);}
            }
          if ((Ft.getType().equals("gene")))
            {
            if((Ft.getAnnotation().containsProperty("note"))&&(Ft.getAnnotation().containsProperty("gene"))){
              String geneN = Ft.getAnnotation().getProperty("gene").toString();
              String geneNote = Ft.getAnnotation().getProperty("note").toString();
              //System.out.println(geneN+"\t"+geneNote);
               GeneNoteTable.put(geneN,geneNote);
            }
            }
         }
        System.out.println(numCDSFeat+" CDS features found");
        }

        System.out.println("Gene Note table: ");
        Set se = GeneNoteTable.keySet();
        for(Iterator i=se.iterator();i.hasNext();)
        {
          String ss = (String)i.next();
          System.out.println(GeneNoteTable.get(ss));
        }

        if(fastaFormat) while(seqI.hasNext()) {
        Sequence seq = seqI.nextSequence();
        System.out.println("Got"+seq.getName());
        String ss = seq.seqString();
        if(ss.indexOf("N")==-1){
             if(maximumGeneLength!=-1)
               if(ss.length()>maximumGeneLength)
                  continue;
             TF = Utils.calcTripletFreq(ss,1,ss.length());
             for(int j=0;j<64;j++)
               GTF[j]+=TF[j];
             NumberOfPoints++;
             ListOfGenes.addElement(seq);
             float CAIV[] = new float[numberOfIterations];
             CAIValues.addElement(CAIV);
        }
        }


        int sum = 0;
        for(int i=0;i<64;i++) sum+=GTF[i];
        for(int i=0;i<64;i++) GTF[i]=GTF[i]/sum;

        System.out.print("CodonUsage"+" sphere ");
        for(int ii=0;ii<64;ii++)
           System.out.print(GTF[ii]+" ");
        System.out.println("0.4 255 0 0");

        float GTFC[] = new float[64];
        System.arraycopy(GTF,0,GTFC,0,64);
        CodonUsages.addElement(GTFC);

        Utils.printDistribution(GTF);
        float GTFAver[] = new float[64];
        for (int i = 0; i < 64; i++) GTFAver[i]=GTF[i];
        GenerateCAIWValuesTable(GTF,"UNIVERSAL");
        GenerateCAIWValuesTableMethodAverage(GTFAver,"UNIVERSAL");
        GenerateCAIWValuesSuppressingRareCodons(GTFAver,"UNIVERSAL",ListOfGenes);
        System.out.println();
        //Utils.printDistribution(GTF);
        //Utils.printDistribution(GTF,"t","c","g","t");
        System.out.println();
        //Utils.printDistribution(GTFAver);

        System.out.println("Number of points = "+NumberOfPoints);

        // --- Calculating all iterations
        Random r = new Random();
        float FirstCAI[] = new float[ListOfGenes.size()];
        float first_p=100.0f;
        float cais[] = new float[ListOfGenes.size()];
        for (int i = 0; i < numberOfIterations; i++) {
        if(i==0){
        for (int j = 0; j < ListOfGenes.size(); j++) {
        if(!fastaFormat){
        if(!checkForStability)
        cais[j]=CalculateCAIValue(((Feature)ListOfGenes.elementAt(j)).getSymbols().seqString(),GTF);
        else
        {    cais[j]=r.nextFloat(); first_p=1f;
        if(classNumber!=-1){
        String name = "**";
        try{
        name = ((Feature)ListOfGenes.elementAt(j)).getAnnotation().getProperty("gene").toString();
        }catch(Exception e){};
        cais[j]=0;
        for(int iii=0;iii<randomclassnumber;iii++)
          if(randomclass[iii][classNumber].equals(name))
            cais[j]=1;
        }
        }
        }else{
        if(!checkForStability)
        cais[j]=CalculateCAIValue(((Sequence)ListOfGenes.elementAt(j)).seqString(),GTF);
        else
        {    cais[j]=r.nextFloat(); first_p=1f; }
        }}
        //for (int ii = 0; ii < cais.length; ii++) System.out.print(cais[ii]+" ");
        //System.out.println();
        for(int hh=0;hh<ListOfGenes.size();hh++) FirstCAI[hh] = cais[hh];
        }

        int prop = (int)(first_p/Math.pow(2,i+1));
        if (prop<1) prop=1;

        System.out.println("Iteration "+i);
        float CodUs[] = new float[64];
        MakeCAIIteration(ListOfGenes,cais,prop,CodUs);

        System.out.print("CodonUsageIt"+i+" sphere ");
        for(int ii=0;ii<64;ii++)
           System.out.print(CodUs[ii]+" ");
        System.out.println("0.2 0 255 0");

        CodonUsages.addElement(CodUs);

        //for (int ii = 0; ii < cais.length; ii++) System.out.print(cais[ii]+" ");
        //System.out.println();


        for (int j = 0; j < ListOfGenes.size(); j++) {
        float mas[] = (float[])CAIValues.elementAt(j);
        mas[i]=cais[j];
        }

        }// ---

        // printing codon usage
        PrintStream cuout = new PrintStream(new FileOutputStream("codonusage"));
        for (int j = 0; j < 61; j++) {
        SymbolList sDNA = DNATools.createDNA(Utils.TripletName(CodonsOrder[j]));
        String trAcid = (RNATools.translate(RNATools.transcribe(sDNA))).seqString();
        cuout.print(trAcid+"\t"+Utils.TripletName(CodonsOrder[j])+"\t");
        for (int i = 0; i < CodonUsages.size(); i++) {
        float cU[] = (float[])(CodonUsages.elementAt(i));
        cuout.print(cU[CodonsOrder[j]]+"\t");
        }
        cuout.println();
        }

        // table for final W-values

        //cuout = new PrintStream(new FileOutputStream("wvalues"));
        int ncu = CodonUsages.size();
        float GTFFinal[] = new float[64];
        //float GTFFinalForClassicCAI[] = new float[64];
        float CodonUsageFinal[] = (float[])(CodonUsages.elementAt(ncu-1));
        for(int i=0;i<64;i++){
        GTFFinal[i]=CodonUsageFinal[i];
        }
        GenerateCAIWValuesTable(GTFFinal,"UNIVERSAL");
        Utils.printDistribution(GTFFinal);
        System.out.println();
        Utils.printDistribution(GTFFinal,"t","c","g","a");
        if(externalWValues!=null){
        System.out.println("External W-table");
        Utils.printDistribution(extWValues,"t","c","g","a");
        }

        // main table
        int addf = 0;
        if(GeneAddInfoCount!=0) addf+=GeneAddInfoCount-1;
        if(addFeature!=null) addf++;
        if(calcGCContent) addf+=2;
        if(externalWValues!=null) addf++;
        //fisa commented out
      //  fout.println((91+(numberOfIterations+addf))+" "+(int)NumberOfPoints);
       
   //     for(int i=0; i<64;i++) fout.println(Utils.TripletName(i)+" FLOAT");
   //     for(int i=0; i<20;i++) fout.println("AA_"+(String)Acids.elementAt(i)+" FLOAT");
   //     fout.println("GENOME"+" STRING");
   //     fout.println("GENENAME"+" STRING");
    //    if(addFeature!=null) fout.println("ADDF"+" STRING");
   //     if(calcGCContent) {
   //       fout.println("GC_CONT"+" FLOAT");
   //       fout.println("GC_CONT3"+" FLOAT");
    //    }
    //    if(externalWValues!=null) fout.println("CAIEXT"+" FLOAT");
    //    fout.println("CAI"+" FLOAT");
   //     fout.println("CAICLASS"+" FLOAT");
   //     fout.println("CAI1"+" FLOAT");
   //     fout.println("CAI2"+" FLOAT");
   //     fout.println("GENELENGTH"+" FLOAT");
   //     for (int i = 0; i < numberOfIterations; i++) fout.println("CAI_IT"+(i+1)+" FLOAT");
   //     for (int i = 1; i < GeneAddInfoCount; i++) fout.println("ADD"+(i)+" STRING");
      // \fisa commented out

/*      eReader = new BufferedReader(
        new InputStreamReader(new FileInputStream(GenBankFile)));
      System.out.println("Loading sequence...");
      seqI = SeqIOTools.readGenbank(eReader);
      System.out.println("Loaded...");*/

     // Here fill table for all genomes
//     while(seqI.hasNext()) {
//        Sequence seq = seqI.nextSequence();

//        for(Iterator i=seq.features(); i.hasNext(); )

          Annotatable Ft = null;

          Date SO = null;
          Date EO = null;

          StringBuffer sOut = new StringBuffer();

          for(int i=0;i<ListOfGenes.size();i++)
          {
          //Date SO1 = new Date();
          //Feature Ft = (Feature) i.next();
          Ft = (Annotatable)ListOfGenes.elementAt(i);
          //if ((Ft.getType().equals("CDS")))
             {
             try{
    //         k=CalcCodonUsageForSeq(Ft.getSymbols(),fr0,fr1,fr2);
             String ss = null;
             if(!fastaFormat)
             ss = ((Feature)Ft).getSymbols().seqString();
             else
             ss = ((Sequence)Ft).seqString();
             //SO = new Date();
             float caiv = FirstCAI[i]; //CalculateCAIValue(ss,GTF);
             float caivclassic = CalculateCAIValue(ss,GTFFinal);
             float caiv1 = CalculateCAIValue(ss,GTFAver);
             float caiv2 = CalculateCAIValue(ss,GTFAverLen);
             //EO = new Date();
             //System.out.println("CAI calc: "+(EO.getTime()-SO.getTime()));
             sOut = new StringBuffer("");

             TF = Utils.calcTripletFreq(ss,1,ss.length());
          //   for (int j=0; j<=63; j++) sOut.append((float)((float)TF[j]/(ss.length()/3))+" "); //fisa commented
             //SO = new Date();
             TF = calcAminoacidFreq(ss,Acids,1,ss.length());
             //EO = new Date();
             //System.out.println("Calc Freq: "+(EO.getTime()-SO.getTime()));

       //      for (int j=0; j<20; j++) sOut.append((float)((float)TF[j]/(ss.length()/3))+" "); //fisa commented

//             for (int j=0; j<=63; j++) sOut.append(TF[j]+" ");
           //  sOut.append("S "); //fisa commented
          //   if (Ft.getAnnotation().containsProperty("gene")){ //fisa commented
          //   sOut.append(" \""+Ft.getAnnotation().getProperty(new String("gene")).toString()+"\" "); //fisa commented
          //    
          //   } //fisa commented
          //   else{ //fisa commented
          //   try{ //fisa commented
          //   sOut.append(" \""+((Sequence)Ft).getName()+"\" ");} //fisa commented
          //   catch(Exception e){ sOut.append("\"Unknown\""); } //fisa commented
          //   }
             try{
             sOut.append(" \""+Ft.getAnnotation().getProperty(new String("product")).toString()+"\" "); //fisa
             }catch(Exception e){ sOut.append("\"Unknown\""); }
             
             try{
             sOut.append(" \""+Ft.getAnnotation().getProperty(new String("old_locus_tag")).toString()+"\" "); //fisa
             }catch(Exception e){ sOut.append("\"Unknown\""); }
             sOut.append(" \""+Ft.getAnnotation().getProperty(new String("db_xref")).toString()+"\" "); //fisa
             if(addFeature!=null){
             String af = new String("");
             if (Ft.getAnnotation().containsProperty("gene"))
             if (Ft.getAnnotation().containsProperty(addFeature))
                af = new String(Ft.getAnnotation().getProperty(addFeature).toString());
             af = af.replace(' ','_');
             if(af.length()>100) af = af.substring(0,100);
             af = " \""+af+"\" ";
           //  sOut.append(af); //fisa commented
             }
             if(calcGCContent){
             //SO = new Date();
             float gc = Utils.GCContent(ss);
             float gc3 = Utils.GCContent3(ss);
             //EO = new Date();
             //System.out.println("GC Cont: "+(EO.getTime()-SO.getTime()));

           //  sOut.append(" "+gc+" "+gc3+" "); //fisa commented
             }
             if(externalWValues!=null) {
             float caivext = CalculateCAIValue(ss,extWValues);
          //   sOut.append(caivext+" ");//fisa commented
             }
          //   sOut.append(caiv+" ");//fisa commented
          //   sOut.append(caivclassic+" ");//fisa commented
         //    sOut.append(caiv1+" ");//fisa commented
         //    sOut.append(caiv2+" ");//fisa commented
        //     sOut.append(ss.length()+" ");//fisa commented
             for (int ii = numberOfIterations-3; ii < numberOfIterations; ii++) {//fisa changed 0 to numberOfIterations-3
                  float mas[] = (float[])CAIValues.elementAt(i);
                  sOut.append(mas[ii]+" "); }
             String geneName = new String("Unknown");
             String geneNote = new String("Unknown");
             if (Ft.getAnnotation().containsProperty("gene"))
               geneName = new String(Ft.getAnnotation().getProperty(new String("gene")).toString());
//             if (Ft.getAnnotation().containsProperty("note"))
//               {
//               geneNote = new String(Ft.getAnnotation().getProperty(new String("note")).toString());
//               System.out.println(geneName+" "+geneNote);
//               }
             if(GeneNoteTable.containsKey(geneName)) geneNote = (String)(GeneNoteTable.get(geneName));
             else
             geneNote = new String("");

             if ((GeneAddInfo.containsKey(geneName))||(GeneAddInfo.containsKey(geneNote))){
             String sm[];
             if (GeneAddInfo.containsKey(geneName))
             sm = (String [])(GeneAddInfo.get(geneName));
             else{
             sm = (String [])(GeneAddInfo.get(geneNote));
             //System.out.println("!!!"+geneNote);
             }
             //GeneAddInfo.get()
             for (int ii = 1; ii < GeneAddInfoCount; ii++) {
                  //System.out.println(sm[ii].substring(1,2)+"//"+sm[ii].substring(0,1));
                  if(sm[ii]==null) {sOut.append("\" \""); continue; }
                  if((sm[ii].length()!=0)&&(sm[ii].substring(0,1).equals("\"")))
                      sOut.append(sm[ii]+" ");
                  else
                      sOut.append("\""+sm[ii]+"\" ");
                  }}
             else
             for (int ii = 1; ii < GeneAddInfoCount; ii++)
                  sOut.append("\"\" ");
             fout.println(sOut.toString());
             //Date EO1 = new Date();
             //System.out.println("Total("+ss.length()+"): "+(EO1.getTime()-SO1.getTime()));

              } catch (Exception e)
            {/**System.out.print("\n"+e+" : ");
             System.out.println(Ft);**/
              System.out.println(sOut.toString());
              e.printStackTrace();}
            }
         }


//     }

    } catch (Throwable t) {
      t.printStackTrace();
      System.exit(1);
    }
  }
  }


public static void GenerateCAIWValuesTable(float WW[], String CodeName){
TranslationTable TranTable = RNATools.getGeneticCode(CodeName);
Alphabet Codons = TranTable.getSourceAlphabet();

HashMap Mp = new HashMap();

try{
for(int i=0;i<64;i++){
   SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
   //System.out.println(sDNA.seqString());
   Symbol sb = ((FiniteAlphabet)Codons).getSymbol((RNATools.transcribe(sDNA)).toList());
   //System.out.println(sb.toString()+" => "+TranTable.translate(sb));
   float f=-1;
   if (Mp.containsKey(TranTable.translate(sb)))
      f = ((Float)Mp.get(TranTable.translate(sb))).floatValue();
   if (WW[i]>f)
      Mp.put(TranTable.translate(sb),new Float(WW[i]));
   }
for(int i=0;i<64;i++){
   SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
   Symbol sb = ((FiniteAlphabet)Codons).getSymbol((RNATools.transcribe(sDNA)).toList());
   if (((Float)Mp.get(TranTable.translate(sb))).floatValue()!=0)
   WW[i]/=((Float)Mp.get(TranTable.translate(sb))).floatValue();
   else WW[i]=0;
}

}catch(Exception e){ System.out.println(e); }
}

public static void GenerateCAIWValuesTableRegardingGeneLength(float WW[], String CodeName,Vector V){
try{
float WWW[][] = new float[V.size()][64];

for(int i=0;i<V.size();i++){
Feature f = (Feature)(V.elementAt(i));
String gs = f.getSymbols().seqString();
int Tc[] = Utils.calcTripletFreq(gs,1,gs.length());
float Tcf[] = new float[64];
for(int j=0;j<64;j++) Tcf[j]=Tc[j]/(gs.length()/3f);
GenerateCAIWValuesTable(Tcf,CodeName);
for(int j=0;j<64;j++) WWW[i][j]=Tcf[j];
}
for(int j=0;j<64;j++) WW[j]=0;
for(int i=0;i<V.size();i++)
  for(int j=0;j<64;j++)
      if(WWW[i][j]!=0) WW[j]+=Math.log(WWW[i][j]);
for(int j=0;j<64;j++) WW[j]/=V.size();
for(int j=0;j<64;j++) WW[j]=(float)Math.exp(WW[j]);
}catch(Exception e){System.out.println(e);}
}

public static void GenerateCAIWValuesSuppressingRareCodons(float WW[], String CodeName,Vector V){
int occ[] = new int[64];
try{
for(int i=0;i<V.size();i++){
Feature f = null;
String gs = null;
if(!fastaFormat){
f = (Feature)(V.elementAt(i));
gs = f.getSymbols().seqString();}
else
gs = ((Sequence)(V.elementAt(i))).seqString();
int Tc[] = Utils.calcTripletFreq(gs,1,gs.length());
for(int j=0;j<64;j++)
  if(Tc[j]!=0) occ[j]++;
}
GenerateCAIWValuesTable(WW,CodeName);
for(int j=0;j<64;j++){
   //WW[j]/=V.size()-occ[j]+1;
   WW[j]*=(float)occ[j]/(float)V.size();
   GTFAverLen[j]=WW[j];
   }
}catch(Exception e){System.out.println(e);}
}


public static void GenerateCAIWValuesTableMethodAverage(float WW[], String CodeName){
TranslationTable TranTable = RNATools.getGeneticCode(CodeName);
Alphabet Codons = TranTable.getSourceAlphabet();

HashMap Mp = new HashMap();
HashMap MpN = new HashMap();

try{
for(int i=0;i<64;i++){
   SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
   //System.out.println(sDNA.seqString());
   Symbol sb = ((FiniteAlphabet)Codons).getSymbol((RNATools.transcribe(sDNA)).toList());
   //System.out.println(sb.toString()+" => "+TranTable.translate(sb));
   float f=0;
   int nn=0;
   if (Mp.containsKey(TranTable.translate(sb)))
      f = ((Float)Mp.get(TranTable.translate(sb))).floatValue();
   if (MpN.containsKey(TranTable.translate(sb)))
      nn = ((Integer)MpN.get(TranTable.translate(sb))).intValue();
   Mp.put(TranTable.translate(sb),new Float(WW[i]+f));
   MpN.put(TranTable.translate(sb),new Integer(nn+1));
//   System.out.println(sb.getName()+" "+(WW[i]+f)+" "+(nn+1));
   }
for(int i=0;i<64;i++){
   SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
   Symbol sb = ((FiniteAlphabet)Codons).getSymbol((RNATools.transcribe(sDNA)).toList());
   if (((Float)Mp.get(TranTable.translate(sb))).floatValue()!=0)
      WW[i]/=((Float)Mp.get(TranTable.translate(sb))).floatValue();
   else WW[i]=0;
}

}catch(Exception e){ System.out.println(e); }
}

public static int[] calcAminoacidFreq(String s, Vector Acids, int st, int en){
    int[] Res = new int[20];

try{

    String ss = s.substring(st-1,en-1);
    int len = ss.length();
    SymbolList Tripl;
    for(int i=0; i<len/3; i++)
       {
       String sTr = ss.substring(i*3,i*3+3);
       int n = 0;
       //System.out.println(Amin.seqString()+" "+sTr);
       int k = Utils.TripletNumber(sTr);
       k = ((Integer)TranslationalTable.elementAt(k)).intValue();
       if (k!=-1) {
          Res[k]+=1;
          }
       }
} catch(Exception e){System.out.println(e);};
    return Res;
}

public static int[] orderCodonUsageByAminoacid(Vector Acids){
int Res[] = new int[64];
try{
int k=0;

for(int j=0;j<Acids.size();j++){

String acName = new String((String)Acids.elementAt(j));

for(int i=0;i<64;i++){
SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
String trAcid = (RNATools.translate(RNATools.transcribe(sDNA))).seqString();
if(trAcid.equals(acName)){
   Res[k] = i;
   k++;
   }
}
}
}catch(Exception e){System.out.println(e+" in  orderCodonUsageByAminoacid"); }
return Res;
}

public static float CalculateCAIValue(String s,float CAIWValuesTable[]){
int len = s.length();
float res = 0;
int L=0;
for(int i=0; i<len/3; i++)
  {
  String sTr = s.substring(i*3,i*3+3);
  int n = Utils.TripletNumber(sTr);
  if (CAIWValuesTable[n]!=0)
    {
    res+=Math.log(CAIWValuesTable[n]);
    L++;
    }
  else res+=Math.log(0.01);
  }
return (float)Math.exp((double)(1.0/(float)L*res));
}

public static void MakeCAIIteration(Vector ListOfGenes,float cais[],int prop,float CodUs[]){
try{

float memcais[] = new float[cais.length];
for (int i = 0; i < cais.length; i++) {
memcais[i]=cais[i];
}

int geneNumber = (int)(0.01*prop*ListOfGenes.size());
System.out.println("First "+geneNumber+" genes ("+prop+"%)");
int nums[] = SortCais(cais);
//for (int i = 0; i < nums.length; i++) System.out.print(cais[nums[i]]+" ");
//System.out.println();
Vector BiggestGenes = new Vector();
for (int i = 0; i < geneNumber; i++)
  BiggestGenes.addElement(ListOfGenes.elementAt(nums[i]));

// Calculate new codon usage

      int TF[];
      float GTF[] = new float[64];

for(int i=0;i<BiggestGenes.size();i++)
{
Feature Ft = null;
String ss = null;
if(!fastaFormat){
Ft = (Feature)BiggestGenes.elementAt(i);
ss = new String(Ft.getSymbols().seqString());}
else
ss = ((Sequence)(BiggestGenes.elementAt(i))).seqString();
TF = Utils.calcTripletFreq(ss,1,ss.length());
for(int j=0;j<64;j++)
   GTF[j]+=TF[j];
}

int sum = 0;
for(int i=0;i<64;i++) sum+=GTF[i];
for(int i=0;i<64;i++) GTF[i]=GTF[i]/sum;
for(int i=0;i<64;i++) CodUs[i]=GTF[i];

// Calculate new cai table

Utils.printDistribution(GTF);
if (taskModifier==0)
GenerateCAIWValuesTable(GTF,"UNIVERSAL");
if (taskModifier==1)
GenerateCAIWValuesTableMethodAverage(GTF,"UNIVERSAL");
if (taskModifier==2)
GenerateCAIWValuesTableRegardingGeneLength(GTF,"UNIVERSAL",BiggestGenes);
if (taskModifier==3)
GenerateCAIWValuesSuppressingRareCodons(GTF,"UNIVERSAL",BiggestGenes);

// Calculate new cais

for(int i=0;i<ListOfGenes.size();i++)
{
Feature Ft = null;
String ss = null;
if(!fastaFormat){
Ft = (Feature)ListOfGenes.elementAt(i);
ss = new String(Ft.getSymbols().seqString());}
else
ss = ((Sequence)ListOfGenes.elementAt(i)).seqString();
    float caiv = CalculateCAIValue(ss,GTF);
    cais[i]=caiv;
}

} catch (Exception e){System.out.print("\n"+e);}

}

public static int[] SortCais(float cais[]){
int res[]=new int[cais.length];
for (int i = 0; i < res.length; i++) res[i]=i;

int i,j,k,inc,n=cais.length;
float v;

inc=1;
do {
	inc *= 3;
	inc++;
} while (inc <= n);

do {
	inc /= 3;
	for (i=inc+1;i<=n;i++) {
		v=cais[res[i-1]];
		j=i;
                k=res[i-1];
		while (cais[res[j-inc-1]]<v) {
			//cais[j]=cais[j-inc];
                        res[j-1]=res[j-inc-1];
			j -= inc;
			if (j <= inc) break;
		}
		//cais[j]=v;
                res[j-1]=k;
	}
} while (inc > 0);

return res;
}

}
