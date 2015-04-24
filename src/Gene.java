import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class Gene {
	private Chromosome chromosome;
	private int start;
	private int end;
	private boolean straight;

	public Gene (Chromosome chromosome, int start, int end, boolean straight){
		this.chromosome = chromosome;
		this.start = start;
		this.end = end;
		this.straight = straight;
	}

	public Chromosome getChromosome() {
		return chromosome;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public boolean isStraight() {
		return straight;
	}

	/** Returns 1 if the gene lies on the leading strand, 
	 * 0 if it lies on the lagging 
	 * and 2 if its place can't be defined. 
	 * Returns 3 if something is wrong.	
	 * @return
	 */
	public int Strand(){
		int conf = this.chromosome.getConfiguration();
		int ori_start = this.chromosome.getOri_start();
		int ori_end = this.chromosome.getOri_end();
		int ter_start = this.chromosome.getTer_start();
		int ter_end = this.chromosome.getTer_end();
		int chr_size = this.chromosome.getSize();
		int start = this.start;
		int end = this.end;
		boolean straight = this.straight;

		int strand = geneStrand(conf, ori_start, ori_end, ter_start, ter_end, chr_size, start, end, straight);

		return strand;
	}


	public int Strand(int conf, int ori_start, int ori_end, int ter_start, int ter_end, int chr_size){
		int start = this.start;
		int end = this.end;
		boolean straight = this.straight;

		int strand = geneStrand(conf, ori_start, ori_end, ter_start, ter_end, chr_size, start, end, straight);

		return strand;
	}

	public static int geneStrand(int conf, int ori_start, int ori_end, int ter_start, int ter_end, int chr_size, int start, int end, boolean straight){
		int strand = 3;

		if (conf == 1){

			if (start > ori_end & end < ter_start){
				if (straight) strand = 1;
				else strand = 0;		
			}
			else if (start > ter_end & end <= chr_size){
				if (straight) strand = 0;
				else strand = 1;		
			}
			else if (start < ori_start & end < ori_start){
				if (straight) strand = 0;
				else strand = 1;		
			}
			else if (start > ter_end & end < ori_start){
				if (straight) strand = 0;
				else strand = 1;		
			}

			else strand = 2;

		}

		if (conf == 2){

			if (start > ori_end & end < ter_start){
				if (straight) strand = 1;
				else strand = 0;		
			}
			else if (start > ter_end & end < ori_start){
				if (straight) strand = 0;
				else strand = 1;		
			}

			else strand = 2;

		}

		if (conf == 3){

			if (start < ter_start & end < ter_start){
				if (straight) strand = 1;
				else strand = 0;		
			}
			else if (start > ter_end & end < ori_start){
				if (straight) strand = 0;
				else strand = 1;		
			}
			else if (start > ori_end & end <= chr_size){
				if (straight) strand = 1;
				else strand = 0;		
			}
			else if (start > ori_end & end < ter_start){
				if (straight) strand = 1;
				else strand = 0;		
			}
			else strand = 2;

		}
		return strand;
	}

	public static ArrayList<String> listFilesForFolder(File folder) {
		ArrayList<String> arr = new ArrayList<String>();
		for (String fileEntry : folder.list(new FilenameFilter() {
			public boolean accept(File dir, String filename) {
				if (filename.charAt(0) == '.'){
					return false;
				}
				return !filename.endsWith(".pl");
			}
		})) {

			arr.add(fileEntry);

		}
		return arr;
	}

	public static void main (String[] args){
		// Folder with perl-generated full gb files
		File folder = new File("C:/Users/weidewind/workspace/asymmetry/");
		HashMap<String, Chromosome> hash = new HashMap<String, Chromosome>();
		try{
		// File with ori data from DORIC	
		BufferedReader br = new BufferedReader(new FileReader("C:/Users/weidewind/Documents/Asymmetry/2015/RawDataSets/ORI_DB.txt"));
		br.readLine();
		while(br.ready()){
			String str = br.readLine();
			String[] oriLine = str.split("\t");
			String id = oriLine[1].split("[_\\.]")[1];
			int oriStart = Integer.parseInt(oriLine[4].split("\\.+")[0]);
			int oriEnd = Integer.parseInt(oriLine[4].split("[\\.\\s]+")[1]);
			int ter = Integer.parseInt(oriLine[6].split("\\s++")[3]);
			int length = Integer.parseInt(oriLine[9].split("\\s++")[0]);
			boolean isMultiOri = false;
			if (hash.containsKey(id)){

				isMultiOri = true;
			}
			
			Chromosome c = new Chromosome(length, oriStart, oriEnd, ter, true, isMultiOri);

			hash.put(id, c);
		}
		}catch (FileNotFoundException e) {

			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (String filename: listFilesForFolder(folder)){

			try {
				BufferedReader br = new BufferedReader(new FileReader("C:/Users/weidewind/workspace/asymmetry/" + filename));
				PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("C:/Users/weidewind/Documents/Asymmetry/2015/RawDataSets/" + filename + "_strand.txt")));
				

				String line = "";
				String id = filename.substring(filename.lastIndexOf('_')+1); 
				
				Chromosome chr = hash.get(id);
				if (chr.isMultiOri()){
					System.out.println(id);
				}
				
//				BufferedReader br = new BufferedReader(new FileReader("C:/Users/weidewind/workspace/asymmetry/Bacillus_subtilis_subsp._NC_000964"));
//				PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("C:/Users/weidewind/Documents/Asymmetry/2015/RawDataSets/Bacillus_subtilis_subsp._NC_000964_strand.txt")));
//
//				int oriStart = 4214414;
//				int oriEnd = 1938;
//				int ter =   1941646;
//				int length = 4214630;
//				boolean isMultiOri = false;
//				
//				Chromosome chr = new Chromosome(length, oriStart, oriEnd, ter, true, isMultiOri);
				
				chr.expandOri(5000);
				chr.expandTer(5000);
				while(br.ready()){
					line = br.readLine().trim();
					if (line.equals("")) break;
					String[] info  = line.split("\\s++");
					boolean straight = true;
					if (Integer.parseInt(info[5])==1) straight = true;
					if (Integer.parseInt(info[5])==-1) straight = false;
					Gene gene = new Gene(chr, Integer.parseInt(info[3]), Integer.parseInt(info[4]), straight);
					out.println(line + "\t" + gene.Strand());

				}
				out.close();
				br.close();

			} catch (FileNotFoundException e) {

				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
//	}

}