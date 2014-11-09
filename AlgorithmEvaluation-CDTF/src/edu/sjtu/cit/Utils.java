package edu.sjtu.cit;

import java.io.*;
import java.util.StringTokenizer;

import Jama.Matrix;

public class Utils {

  /**
   * Read matrix from file in complete format.
   * 
   * @param file  matrix file.
   * @param x     matrix to-be-load.
   * @throws IOException
   */
  public static void readFullMatrix(String file, Matrix x) throws IOException {
    BufferedReader br = new BufferedReader(new FileReader(new File(file)));
    
    String line = "";
    int row = 0, col = 0;
    while((line = br.readLine()) != null) {
      col = 0;
      StringTokenizer st = new StringTokenizer(line, ",");
      while(st.hasMoreTokens()) {
        x.set(row, col, Double.parseDouble(st.nextToken()));
        ++col;
      }
      ++row;
    }
    
    br.close();
  }
  
  
  /**
   * Write matrix to file in complete format.
   * 
   * @param file  matrix file.
   * @param x     matrix to-be-write.
   * @throws IOException
   */
  public static void writeFullMatrix(String file, Matrix x) throws IOException {
    PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(file)), true);
    int m = x.getRowDimension(), n = x.getColumnDimension();
    
    for(int i=0; i<m; i++) {
      for(int j=0; j<n-1; j++) writer.printf("%f,", x.get(i, j));
      writer.printf("%f\n", x.get(i, n-1));
    }
    
    writer.flush();
    writer.close();
  }
  
  
  public static void readCSVIndex(String file, int[][] index) throws IOException {
    BufferedReader br = new BufferedReader(new FileReader(new File(file)));
    
    String line = "";
    int cnt = 0;
    while((line = br.readLine()) != null) {
      String[] sp = line.split(",");
      int u = Integer.parseInt(sp[0]) - 1,
          i = Integer.parseInt(sp[1]) - 1;
      index[cnt][0] = u;
      index[cnt][1] = i;
      ++cnt;
    }
    
    br.close();
  }
  
  /**
   * Read matrix from file in `<user,item,rating>` format.
   * 
   * @param file    matrix file.
   * @param x       matrix to-be-load.
   * @param binary  if implict or explicit.
   * @param scale   if scale rating.
   * @throws IOException
   */
  public static int readCSVMatrix(String file, Matrix x, boolean binary, double scale) throws IOException {
    BufferedReader br = new BufferedReader(new FileReader(new File(file)));
    
    String line = "";
    int cnt = 0;
    while((line = br.readLine()) != null) {
      String[] sp = line.split(",");
      int u = Integer.parseInt(sp[0]) - 1,
          i = Integer.parseInt(sp[1]) - 1;
      double r = Double.parseDouble(sp[2]);
      
      if(binary) r = r > 0 ? 1:0;
      else r /= scale;
      
      x.set(u, i, r);
      ++cnt;
    }
    
    br.close();
    return cnt;
  }
  
  /**
   * Read matrix from file in `<user,item,rating>` format.
   * 
   * @param file    matrix file.
   * @param x       matrix to-be-load.
   * @param binary  if implict or explicit.
   * @param scale   if scale rating.
   * @param alpha	IF parameter alpha
   * @throws IOException
   */
  public static int readCSVMatrix(String file, Matrix x, boolean binary, double scale, double alpha) throws IOException {
    BufferedReader br = new BufferedReader(new FileReader(new File(file)));
    
    String line = "";
    int cnt = 0;
    while((line = br.readLine()) != null) {
      String[] sp = line.split(",");
      int u = Integer.parseInt(sp[0]) - 1,
          i = Integer.parseInt(sp[1]) - 1;
      double r = Double.parseDouble(sp[2]);
      
      if(binary) r = r > 0 ? 1:0;
      else r /= scale;
      
      x.set(u, i, r*alpha);
      ++cnt;
    }
    
    br.close();
    return cnt;
  }
}
