import java.io.*;
import java.util.*;
import java.text.*;

import org.biojava.bio.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.io.*;


public class Utils {

  public Utils() {
  }
    public static int[] calcTripletFreq(SymbolList seq, int st, int en) throws Exception
    {
    int[] Res = new int[64];
    SymbolList LS = seq.subList(st,en);
    int len = LS.length();
    SymbolList Tripl;
    for(int i=0; i<len/3; i++)
       {
       Tripl = LS.subList(i*3+1,i*3+3);
       String sTr = Tripl.seqString();
       int n = TripletNumber(sTr.toString());
       Res[n]+=1;
       }
    return Res;
    }

    public static int[] calcTripletFreq(String s, int st, int en) throws Exception
    {
    int[] Res = new int[64];
    //SymbolList LS = seq.subList(st,en);
    String ss = s.substring(st-1,en-1);
    int len = ss.length();
    SymbolList Tripl;
    for(int i=0; i<len/3; i++)
       {
       //Tripl = LS.subList(i*3+1,i*3+3);
       String sTr = ss.substring(i*3,i*3+3);
       int n = TripletNumber(sTr.toString());
       Res[n]+=1;
       }
    return Res;
    }

    public static float[] calcTripletFreqF(SymbolList seq, int pos, int wSize, int Bayez) throws Exception
    {
    float[] Res = new float[64];
    int nt = (int)(wSize/6);
    int ntb = 0;
    if (Bayez>1)
       ntb = (int)(wSize/6/Bayez);
    nt = (int)(nt/2)*2+1;
    ntb = (int)(ntb/2)*2+1;
    int st = pos-1-(nt-1)*3, en = pos+1+(nt-1)*3;
    int stb = pos-1-(ntb-1)*3, enb = pos+1+(ntb-1)*3;
    int fi[] = calcTripletFreq(seq,st,en);
    int fib[] = new int[64];
    if (Bayez>1)
      fib = calcTripletFreq(seq,stb,enb);
    for(int i=0; i<64; i++)
       {
       if (Bayez>1)
       Res[i]=(((float)fi[i]+(float)fib[i])/(nt+ntb));
       else
       Res[i]=(float)fi[i]/nt;
       }
    return Res;
    }

    public static float[] calcTripletFreqF(String s, int pos, int wSize, int Bayez) throws Exception
    {
    float[] Res = new float[64];
    int nt = (int)(wSize/6);
    int ntb = 0;
    if (Bayez>1)
       ntb = (int)(wSize/6/Bayez);
    nt = (int)(nt/2)*2+1;
    ntb = (int)(ntb/2)*2+1;
    int st = pos-1-(nt-1)*3, en = pos+1+(nt-1)*3;
    int stb = pos-1-(ntb-1)*3, enb = pos+1+(ntb-1)*3;
    int fi[] = calcTripletFreq(s,st,en);
    int fib[] = new int[64];
    if (Bayez>1)
      fib = calcTripletFreq(s,stb,enb);
    for(int i=0; i<64; i++)
       {
       if (Bayez>1)
       Res[i]=(((float)fi[i]+(float)fib[i])/(nt+ntb));
       else
       Res[i]=(float)fi[i]/nt;
       }
    return Res;
    }


/*    public static float[] calcTripletFreqF(String s, int st, int en) throws Exception
    {
    int[] Res = new int[64];
    //SymbolList LS = seq.subList(st,en);
    String ss = s.substring(st-1,en-1);
    int len = ss.length();
    SymbolList Tripl;
    for(int i=0; i<len/3; i++)
       {
       //Tripl = LS.subList(i*3+1,i*3+3);
       String sTr = ss.substring(i*3,i*3+3);
       int n = TripletNumber(sTr.toString());
       Res[n]+=1;
       }
    return Res;
    }*/

    public static int[] calcTripletFreqComplementary(SymbolList seq, int st, int en) throws Exception
    {
    int[] Res = new int[64];
    SymbolList LS = seq.subList(st,en);
    int len = LS.length();
    SymbolList Tripl;
    for(int i=0; i<len/3; i++)
       {
       Tripl = DNATools.complement(LS.subList(i*3+1,i*3+3));
       StringBuffer sTr = new StringBuffer(Tripl.seqString());
       int n = TripletNumber(sTr.reverse().toString());
       Res[n]+=1;
       }
    return Res;
    }

    public static int TripletNumber(String st)
    {
    char[] cA = new char[3];
    st.getChars(0,3,cA,0);
    return BaseNumber(new String(cA,0,1))*16+BaseNumber(new String(cA,1,1))*4+BaseNumber(new String(cA,2,1));
    }

    public static int TripletNumber(String st, String let1, String let2, String let3, String let4)
    {
    char[] cA = new char[3];
    st.getChars(0,3,cA,0);
    return BaseNumber(new String(cA,0,1),let1,let2,let3,let4)*16+
           BaseNumber(new String(cA,1,1),let1,let2,let3,let4)*4+
           BaseNumber(new String(cA,2,1),let1,let2,let3,let4);
    }

    public static String BaseName(int n){
    String s = new String("");
    switch(n){
    case 0: s = new String("a"); break;
    case 1: s = new String("c"); break;
    case 2: s = new String("g"); break;
    case 3: s = new String("t"); break;
    }
    return s;
    }

    public static String BaseName(int n, String let1, String let2, String let3, String let4){
    String s = new String("");
    switch(n){
    case 0: s = new String(let1); break;
    case 1: s = new String(let2); break;
    case 2: s = new String(let3); break;
    case 3: s = new String(let4); break;
    }
    return s;
    }

    public static String TripletName(int n){
    StringBuffer s  = new StringBuffer();
    int i1 = (int)((float)n/16);
    int i2 = (int)((float)(n-i1*16)/4);
    int i3 = n-i1*16-i2*4;
    s.append(BaseName(i1));
    s.append(BaseName(i2));
    s.append(BaseName(i3));
    return s.toString();
    }

    public static String TripletName(int n, String let1, String let2, String let3, String let4){
    StringBuffer s  = new StringBuffer();
    int i1 = (int)((float)n/16);
    int i2 = (int)((float)(n-i1*16)/4);
    int i3 = n-i1*16-i2*4;
    s.append(BaseName(i1,let1,let2,let3,let4));
    s.append(BaseName(i2,let1,let2,let3,let4));
    s.append(BaseName(i3,let1,let2,let3,let4));
    return s.toString();
    }

    public static int BaseNumber(String cn)
    {
    int Res = 0;
    if (cn.equals("a")) Res = 0;
    if (cn.equals("c")) Res = 1;
    if (cn.equals("g")) Res = 2;
    if (cn.equals("t")) Res = 3;
    return Res;
    }

    public static int BaseNumber(String cn, String let1, String let2, String let3, String let4)
    {
    int Res = 0;
    if (cn.equals(let1)) Res = 0;
    if (cn.equals(let2)) Res = 1;
    if (cn.equals(let3)) Res = 2;
    if (cn.equals(let4)) Res = 3;
    return Res;
    }


    public static float GCContent(String cn){
    float gc = 0;
    for(int i=0;i<cn.length();i++){
    if((cn.charAt(i)=='g')||(cn.charAt(i)=='c')) gc+=1;
    }
    return gc/cn.length();
    }

    public static float GCContent3(String cn){
    float gc = 0;
    for(int i=2;i<cn.length();i+=3){
    if((cn.charAt(i)=='g')||(cn.charAt(i)=='c')) gc+=1;
    }
    return gc/(cn.length()/3.0f);
    }

    public static int CalcIsCoding(Sequence seq, int pos)
    {
    int Res = 0;
    /*for(Iterator i=seq.features(); (i.hasNext())&&(Res==0); )
      {
      Feature Ft = (Feature) i.next();
      if ((Ft.getType().equals("CDS"))&&(Ft.getLocation().contains(pos))) Res=1;
      }*/
    String s = GetFeatureName(seq,pos);
    if((s.indexOf("CDS")>-1)&&(s.indexOf("intron")==-1)) Res = 1;
    return Res;
    }

    public static int CalcIsCodingGene(Sequence seq, int pos)
    {
    int Res = 0;
    /*for(Iterator i=seq.features(); (i.hasNext())&&(Res==0); )
      {
      Feature Ft = (Feature) i.next();
      if ((Ft.getType().equals("CDS"))&&(Ft.getLocation().contains(pos))) Res=1;
      }*/
    String s = GetFeatureName(seq,pos);
    if(s.indexOf("gene")>-1) Res = 1;
    return Res;
    }

    public static int CalcIsCoding2(Sequence seq, int pos)
    {
    int Res = 0;
    /*for(Iterator i=seq.features(); (i.hasNext())&&(Res==0); )
      {
      Feature Ft = (Feature) i.next();
      if ((Ft.getType().equals("CDS"))&&(Ft.getLocation().contains(pos))) Res=1;
      }*/
    String s = GetFeatureName(seq,pos);
    if(((s.indexOf("CDS")>-1)||(s.indexOf("tRNA")>-1)||(s.indexOf("rRNA")>-1))&&(s.indexOf("intron")==-1)) Res = 1;
    return Res;
    }

    public static int CalcIsGeneComplement(Sequence seq, int pos)
    {
    int Res = 0;
    for(Iterator i=seq.features(); (i.hasNext())&&(Res==0); )
      {
      Feature Ft = (Feature) i.next();
      if ((Ft.getType().equals("CDS"))&&(Ft.getLocation().contains(pos)))
         {
         Res=1;
         for (Iterator ai = Ft.getAnnotation().asMap().entrySet().iterator(); ai.hasNext(); ) {
	 Map.Entry me = (Map.Entry) ai.next();
         if (me.getKey().toString().equals("internal_data"))
         if (me.getValue().toString().indexOf("complement")>=0)
            Res=-1;
         }
         }
      }
    return Res;
    }

    public static int NumberOfBase(StringBuffer seq, char b){
    int c=0;
    for(int i=0;i<seq.length();i++)
      if (seq.charAt(i)==b) c++;
    return c;
    }

    public static String GetFeatureName(Sequence seq, int pos)
    {
    int Res = 0;
    String sRes = "";
    for(Iterator i=seq.features(); (i.hasNext())&&(Res==0); )
      {
      Feature Ft = (Feature) i.next();
      if (Ft.getLocation().contains(pos))
         {
         //Res=1;
         if (!(Ft.getType().equals("source")))
           {
           if (sRes.equals(""))
           sRes = new String(Ft.getType());
           else
           sRes = new String(sRes+","+Ft.getType());
           }
         }
      }
    return sRes;
    }


    public static String ComplementaryTripletName(int n){
    StringBuffer s  = new StringBuffer();
    int i1 = (int)((float)n/16);
    int i2 = (int)((float)(n-i1*16)/4);
    int i3 = n-i1*16-i2*4;
    s.append(ComplementaryBaseName(i1));
    s.append(ComplementaryBaseName(i2));
    s.append(ComplementaryBaseName(i3));
    return s.toString();
    }

    public static String ComplementaryBaseName(int n){
    String s = new String("");
    switch(n){
    case 0: s = new String("t"); break;
    case 1: s = new String("g"); break;
    case 2: s = new String("c"); break;
    case 3: s = new String("a"); break;
    }
    return s;
    }


public static float FreqReconstruct(float fr[], String codon,int phase){
float f = 0;
//System.out.println(codon);
for(int i=0;i<64;i++)
  for(int j=0;j<64;j++)
   {
   String hexa = new String(TripletName(i)+TripletName(j));
   if (hexa.substring(phase,phase+3).equals(codon)){
          //System.out.println(hexa+" "+fr[i]+" "+fr[j]+" "+f);
          f+=fr[i]*fr[j];
          }
   }
return f;
}

public static float[] ReconstructPhaseDistribution(float distr[], int phase){
float f[] = new float[64];
for(int i=0;i<64;i++)
  f[i]=FreqReconstruct(distr,TripletName(i),phase);
return f;
}

public static float DistributionDistance(float f1[],float f2[]){
double r=0;
for(int i=0;i<64;i++)
  r+=Math.abs(f1[i]-f2[i]);//*(f1[i]-f2[i]);
// r+=(f1[i]-f2[i])*(f1[i]-f2[i]);
//return (float)Math.sqrt(r);
return (float)r;
}

public static float DistributionNorm(float f1[]){
double r=0;
for(int i=0;i<64;i++)
  r+=Math.abs(f1[i]);//*f1[i];
//  r+=f1[i]*f1[i];
//return (float)Math.sqrt(r);
return (float)r;
}


public static float ReconstructionQuality(float distr[],float realdist[],int phase){
float f[] = ReconstructPhaseDistribution(distr,phase);
return DistributionDistance(realdist,f)/DistributionNorm(distr);
}

public static float EvaluateDistributionGoodness(float f[]){
float p[] = ReconstructPhaseDistribution(f,1);
return DistributionDistance(f,p)/DistributionNorm(f)/2;
}

public static float[] ComplementaryDistribution(float fr[]){
float f[] = new float[64];
for(int i=0;i<64;i++){
  StringBuffer s = new StringBuffer(ComplementaryTripletName(i));
  s.reverse();
  f[TripletNumber(s.toString())] = fr[i];
  }
return f;
}

public static void printDistribution(float f[]){
DecimalFormat df = new DecimalFormat();
df.applyPattern("#,###0.000");
for(int i=0;i<8;i++)
 {
 for(int j=0;j<8;j++){
   System.out.print(TripletName(i*8+j)+":"+df.format(f[i*8+j])+"\t");
   }
 System.out.println();
 }
}

public static void printDistribution(float f[], String let1, String let2, String let3, String let4){
DecimalFormat df = new DecimalFormat();
df.applyPattern("#,###0.000");
for(int i=0;i<8;i++)
 {
 for(int j=0;j<8;j++){
   String s = TripletName(i*8+j,let1,let2,let3,let4);
   int k = TripletNumber(s);
   System.out.print(s+":"+df.format(f[k])+"\t");
   }
 System.out.println();
 }
}


public static float ComplementaryAsymmetry(float fr0[],float fr1[], float fr2[], float cfr0[], float cfr1[], float cfr2[]){
float r=0;
float d00 = DistributionDistance(fr0,cfr0);
float d01 = DistributionDistance(fr0,cfr1);
float d02 = DistributionDistance(fr0,cfr2);
float d10 = DistributionDistance(fr1,cfr0);
float d11 = DistributionDistance(fr1,cfr1);
float d12 = DistributionDistance(fr1,cfr2);
float d20 = DistributionDistance(fr2,cfr0);
float d21 = DistributionDistance(fr2,cfr1);
float d22 = DistributionDistance(fr2,cfr2);
float r1 = Math.min(d00,Math.min(d01,d02));
float r2 = Math.min(d10,Math.min(d11,d12));
float r3 = Math.min(d20,Math.min(d21,d22));
r = (r1+r2+r3)/3;
return r;
}

public static float[] CalcJunkDistribution(float fr0[],float fr1[], float fr2[], float cfr0[], float cfr1[], float cfr2[]){
float res[] = new float[64];
for(int i=0;i<64;i++)
  res[i]=(fr0[i]+fr1[i]+fr2[i]+cfr0[i]+cfr1[i]+cfr2[i])/6;
return res;
}

public static String StringSpaces(int n){
StringBuffer ss = new StringBuffer("");
for(int i=0;i<n;i++)
 ss.append(" ");
return ss.toString();
}

public static String getJunkString(Sequence seq){
StringBuffer ss = new StringBuffer(seq.seqString());
    for(Iterator i=seq.features(); i.hasNext(); )
      {
      Feature Ft = (Feature) i.next ();
      int i1 = Ft.getLocation().getMin();
      int i2 = Ft.getLocation().getMax()+1;
         if (!(Ft.getType().equals("source")))
           if (Ft.getType().equals("gene"))
           {
           ss.replace(i1-1,i2-1,StringSpaces(i2-i1));
           }
      }
StringBuffer s1 = new StringBuffer("");
for(int i=0;i<ss.length();i++)
if (ss.charAt(i)!=' ')
  s1.append(ss.charAt(i));
return s1.toString();
}

public static String FormatSeq(StringBuffer s){
StringBuffer ss = new StringBuffer();
int n = (int)((float)s.length()/60);
for(int i=0;i<n;i++){
ss.append("    "+(i*60+1)+" ");
   for(int j=0;j<6;j++)
     ss.append(s.substring(i*60+j*10,i*60+j*10+10)+" ");
ss.append("\n");
}
return ss.toString();
}


}