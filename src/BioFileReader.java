
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.regex.PatternSyntaxException;
import java.sql.*;

public class BioFileReader {
	
	public static void FastaOneLiner(String infile, String outfile){
		try {
			BufferedReader br = new BufferedReader(new FileReader(infile));
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
			while(br.ready()){
				String str = br.readLine();
				if (str.charAt(0)=='>'){
					out.print("\n");
					continue;
				}
				out.print(str.trim());	
			}
			br.close();
			out.close();
		} catch (FileNotFoundException e) {

			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	public static void GBReader(String infile, String outfile){

		try {
			BufferedReader br = new BufferedReader(new FileReader(infile));
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
			String line = "";
			String check = "";

			while(br.ready()){

				line = br.readLine();
		//		line.trim();
				if (line.substring(0,2).equals("//")|line.trim().equals("")){
					break;
				}
				if (line.length()>8){
					check = line.substring(0, 8);
				}
				else check = "";

				while(!(check.equals("     CDS")|check.equals("     tRN")|check.equals("     rRN"))&br.ready()){
		//			line = br.readLine().trim();
					check = line.substring(0, 8);
				}

				int[] position  = parsePosition(line.trim());
				out.print(line.substring(0, 4).trim() +"\t" + position[2] + "\t" + position[0] + "\t" + position[1]+ "\t");

				String[] features = parseFeatures(br);
				out.print(features[0] + "\t" + features[1] + "\t" + features[2] +"\t" + features[3] + "\t" + features[4] + "\t" + features[5] + "\n");


			}
			br.close();
			out.close();

		} catch (FileNotFoundException e) {

			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	
	
	public static void GBReader(String infile){

		try {
			BufferedReader br = new BufferedReader(new FileReader(infile));
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("dfgd"))); 
			
			Class.forName("com.mysql.jdbc.Driver");
			Connection conn;
			conn = DriverManager.getConnection(
			        "jdbc:mysql://localhost:3306/bio",
			        "root", "qwe");
			
			if (conn == null) {
	            System.out.println("Нет соединения с БД!");
	            System.exit(0);
	        }
			
			
			
			
			String line = "";
			String check = "";
			
			String refseq = br.readLine().trim().split("\\s++")[1];
			
			Statement stmt = conn.createStatement();
			ResultSet rs = stmt.executeQuery("SELECT ori_start FROM ori WHERE refseq LIKE '" + refseq + "%';");
			
			while (rs.next()) {
	            System.out.println(rs.getInt("ori_start"));
	        }
			
			while(br.ready()){

				line = br.readLine();
		//		line.trim();
				if (line.substring(0,2).equals("//")|line.trim().equals("")){
					break;
				}
				if (line.length()>8){
					check = line.substring(0, 8);
				}
				else check = "";

				while(!(check.equals("     CDS")|check.equals("     tRN")|check.equals("     rRN"))&br.ready()){
		//			line = br.readLine().trim();
					check = line.substring(0, 8);
				}
				
                String type = line.substring(0, 4).trim();
				int[] position  = parsePosition(line.trim());
				String[] features = parseFeatures(br);

				
				PreparedStatement preparedStatement = conn
				          .prepareStatement("insert into  gene " +
				          		            "(refseq, gi, locus, id, start, end, straight, product, type, name)" +
				          		            " values (?, ?, ?, ?, ?, ?, ?, ?, ?)");
				     
				      preparedStatement.setString(1, refseq);
				      preparedStatement.setString(2, features[4]);
				      preparedStatement.setString(3, features[1]);
				      preparedStatement.setString(4, features[3]);
				      preparedStatement.setInt(5, position[0]);
				      preparedStatement.setInt(6, position[1]);
				      preparedStatement.setInt(7, position[2]);
				      preparedStatement.setString(8, features[2]);
				      preparedStatement.setString(9, type);
				      preparedStatement.setString(10, features[0]);
				      preparedStatement.executeUpdate();
				
			
				stmt.executeQuery("show warnings;");
				System.out.println(preparedStatement.getWarnings());
				
				out.print(line.substring(0, 4).trim() +"\t" + position[2] + "\t" + position[0] + "\t" + position[1]+ "\t");

				
				out.print(features[0] + "\t" + features[1] + "\t" + features[2] +"\t" + features[3] + "\t" + features[4] + "\t" + features[5] + "\n");


			}
			br.close();
			out.close();
			conn.close();

		} catch (FileNotFoundException e) {

			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}	
	
	
	public static int[] parsePosition(String cds){
		int[] position = new int[4];
		try {
			
			String temp_position = cds.split("\\s++")[1];
			
			if (temp_position.charAt(0)=='c'){
				String[] t = temp_position.split("\\W++");
				position[0] = Integer.parseInt(t[2]);
				position[1] = Integer.parseInt(t[1]);
				position[2] = 0;
			}
			else {
				String[] t = temp_position.split("\\W++");
				position[0] = Integer.parseInt(t[0]);
				position[1] = Integer.parseInt(t[1]);
				position[2] = 1;

			}

		} 
		catch (ArrayIndexOutOfBoundsException e){
			System.out.println("Strange CDS: " + cds);
		}
		return position;
	}

	public static String[] parseFeatures(BufferedReader br){
		String[] features = new String[6];

		String line;
		String check = "";
		try {
			
			line = br.readLine();
//					line.trim();
			
			if (line.length()>4){
				check = line.substring(0, 9);
			}
			else check = "";

			while(!check.equals("     gene")&br.ready()&!line.equals("//")&!line.equals("")){
				
				line.trim();
				if (line.charAt(0)=='/'){
					String[] blocks = line.split("[/=\"]++");
					String feature_type = blocks[1];
					if (feature_type.equals("gene")){
						features[0] = blocks[2];
					}
					if (feature_type.equals("locus_tag")){
						features[1] = blocks[2];
					}
					if (feature_type.equals("product")){
						features[2] = blocks[2];
					}
					if (feature_type.equals("protein_id")){
						features[3] = blocks[2];
					}
					if (feature_type.equals("db_xref")){

						String[] id = line.split("\"");
						if (id[1].substring(0, 2).equals("GI")){
							features[4] = id[1];
						}
						if (id[1].substring(0, 2).equals("Ge")){
							features[5] = id[1];
						}

					}
				}


				if (line.equals("//")|line.equals("")){
					break;
				}
				line = br.readLine().trim();
				try {
					check = line.substring(0, 4);
				} catch (IndexOutOfBoundsException e){
					check = "";
				}
			}
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		return features;
	}




	public static void main (String[] args){
//		//GBReader("src/Files/Acin.gb", "src/Files/Acin.xls");
//		try {
//			Class.forName("com.mysql.jdbc.Driver");
//		} catch (ClassNotFoundException e1) {
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		}
//		Connection conn;
//		try {
//			conn = DriverManager.getConnection(
//			        "jdbc:mysql://localhost:3306/bio",
//			        "root", "qwe");
//			if (conn == null) {
//	            System.out.println("Нет соединения с БД!");
//	            System.exit(0);
//	        }
//			
//			Statement stmt = conn.createStatement();
//			ResultSet rs = stmt.executeQuery("SELECT ori_start FROM ori WHERE refseq LIKE 'NC_002932%';");
//			int counter = 0;
//			while (rs.next()) {
//	            System.out.println(rs.getInt("ori_start"));
//	            counter++;
//	        }
//			
//			if (counter == 0){
//				System.out.println("No such chromosome");
//			}
//			else if (counter == 1){
//				//здесь написать
//			}
//			else {
//				System.out.println("Too many origins");
//			}
//			
//			
//			conn.close();
//			
//		} catch (SQLException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//	 
	    
		FastaOneLiner("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/ExternalResults/n1.all.protein.fa","C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/ExternalResults/n1.all.protein.aln");
		
	}
}
