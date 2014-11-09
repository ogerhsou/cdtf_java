package edu.sjtu.cit.cdtf;

import java.io.*;
import java.util.*;

import org.apache.commons.cli.*;

import Jama.Matrix;
import edu.sjtu.cit.*;

public class Main {
	static final int RANDOM = 0;
	static final int SVD = 1;
	static final PrintStream out = System.out;
  @SuppressWarnings("static-access")
  public static void main(String[] args) throws Exception {
	long tic = System.currentTimeMillis();
    long initS = System.currentTimeMillis();
    
    Options options = new Options();
    
    Option help = OptionBuilder.withLongOpt("help")
                               .withDescription("Show this help.")
                               .create("h");
    
    Option initi = OptionBuilder.withArgName("domain initialization method")
                                .isRequired()
                                .hasArgs()
                                .withLongOpt("initi")
                                .withDescription("Domain initialization configuration.")
                                .create("i");  
    
    Option assignFactors = OptionBuilder.withArgName("assign original value to factors")
							            .hasArgs()
							            .withLongOpt("assign")
							            .withDescription("Assign original value to factors.")
							            .create("a"); 
    
    Option EMFlag = OptionBuilder.withLongOpt("em")
					            .withDescription("Use EM or not.")
					            .create("e");
    
    Option resultFlag = OptionBuilder.withLongOpt("hide")
            .withDescription("Hide fitting result or not.")
            .create("H");
    
    Option domainNum = OptionBuilder.withArgName("number of domains")
						            .isRequired()
						            .hasArgs()
						            .withLongOpt("num")
						            .withDescription("Number of domains.")
						            .create("N");  
    
    Option factorIndicator = OptionBuilder.withArgName("indicators of three factors")
							            .hasArgs()
							            .withLongOpt("Indicator")
							            .withDescription("Indicators of three factors.")
							            .create("I");
    
    Option conv = OptionBuilder.withArgName("convergent value")
								            .hasArgs()
								            .withLongOpt("conv")
								            .withDescription("Convergent value.")
								            .create("c");
    
    Option facRegularWeight = OptionBuilder.withArgName("factor regularization weights")
                                    .isRequired()
                                    .hasArgs()
                                    .withLongOpt("factor")
                                    .withDescription("Factor regularization weights.")
                                    .create("f");
    
    Option R = OptionBuilder.withArgName("R")
                            .hasArg()
                            .isRequired()
                            .withLongOpt("dim")
                            .withDescription("Factor common dimension.")
                            .create("r");
    
    Option alpha = OptionBuilder.withArgName("IF parameter alpha")
            .hasArg()
            .withLongOpt("alpha")
            .withDescription("IF parameter alpha.")
            .create("A");
    
    Option N = OptionBuilder.withArgName("maxIter")
                            .hasArg()
                            .withLongOpt("niter")
                            .withDescription("Maxium iteration number.")
                            .create("n");
    
    Option target = OptionBuilder.withArgName("targetKey")
                                 .hasArg()
                                 .withLongOpt("target")
                                 .withDescription("Target domain key.")
                                 .create("t");
    
    Option binary = OptionBuilder.withLongOpt("binary")
                                 .withDescription("Implicit or Explicit.")
                                 .create("b");
    
    Option test = OptionBuilder.withArgName("test config")
                               .hasArg()
                               .isRequired()
                               .withLongOpt("test")
                               .withDescription("Test source(tableName if database used, or it'll be filename).")
                               .create("T");
    
    Option atK = OptionBuilder.withArgName("topKs")
                              .hasArgs()
                              .withLongOpt("topk")
                              .withDescription("List of different topKs.")
                              .create("k");
    
    Option output = OptionBuilder.withArgName("path")
					            .hasArg()
					            .withLongOpt("output")
					            .withDescription("Output path of factors.")
					            .create("o");
    
    Option log = OptionBuilder.withArgName("log suffix")
					            .hasArg()
					            .withLongOpt("log")
					            .withDescription("Log files' suffix.")
					            .create("l");
					    
    options.addOption(help);
    options.addOption(alpha);
    options.addOption(log);
    options.addOption(resultFlag);
    options.addOption(domainNum);
    options.addOption(conv);
    options.addOption(EMFlag);
    options.addOption(assignFactors);
    options.addOption(factorIndicator);
    options.addOption(initi);
    options.addOption(facRegularWeight);
    options.addOption(R);
    options.addOption(N);
    options.addOption(target);
    options.addOption(binary);
    options.addOption(test);
    options.addOption(atK);
    options.addOption(output);
    
    
    CommandLineParser parser = new PosixParser();
    CommandLine cmd = null;
    try {
      cmd = parser.parse(options, args);
    } catch(MissingOptionException e) {
      HelpFormatter h = new HelpFormatter();
      h.printHelp("CDTF", options);
      
      System.err.println(e.toString());
      return;
    }
    
    if(cmd.hasOption("h")) {
      HelpFormatter h = new HelpFormatter();
      h.printHelp("CDTF", options);
      return;
    }
    
    if(cmd.hasOption("o"))
    	redirect(cmd.getOptionValue("o"));
    
    CDTF cdtf = new CDTF();
    
    int num = Integer.parseInt(cmd.getOptionValue("N"));
    Matrix[] X = new Matrix[num];
    double[] wX = new double[num];	//weight of each matrix
    
    /* initialization inputs */
    boolean biFlag = false;
    if(cmd.hasOption("b")) biFlag = true;
    
    boolean EM = false;
    if(cmd.hasOption("e")) EM = true;
    
    String[] initiStr = cmd.getOptionValues("i");   
    for(String s : initiStr) {
      InputStream is = new FileInputStream(s);
      Properties prop = new Properties();
      prop.load(is);
      is.close();
      
      int pKey = Integer.parseInt(prop.getProperty("primaryKey"));
      double weight = Double.parseDouble(prop.getProperty("weight"));
      String file = prop.getProperty("file");
      String[] szStr = prop.getProperty("size").split(",");
      int sz[] = new int[] { Integer.parseInt(szStr[0]),
          Integer.parseInt(szStr[1]) };
      
      Matrix data = new Matrix(sz[0], sz[1], 0.0);
      
      double IFPara = 1;
	  	if(cmd.hasOption("A"))
	  		IFPara = Double.parseDouble(cmd.getOptionValue("A"));
      
	  System.out.println("The weight of training elements: " + IFPara);
      Utils.readCSVMatrix(file, data, biFlag, 1, IFPara);
      
      X[pKey] = data;  
      wX[pKey] = weight;
      
      System.out.printf("Domain<%d>: %dx%d\n", pKey, sz[0], sz[1]);
      System.out.println("- weight:\t" + weight);
  }
  System.out.println("");
    
    boolean[] fIndicator = {false,false,false};
    double convCrit = 0;
    int assignType = RANDOM;
    
    
    if(cmd.hasOption("I")) {
    	int cnt=0;
    	 String[] fIndiStr = cmd.getOptionValues("I");   
    	    for(String s : fIndiStr) {
    	    	if(s == "0")
    	    		fIndicator[cnt++] = false;
    	    	else
    	    		fIndicator[cnt++] = true;
    	    }
      }
    
    if(cmd.hasOption("c")) {
    	 String convStr = cmd.getOptionValue("c");   
    	 convCrit = Double.parseDouble(convStr);
      }
    
    if(cmd.hasOption("a")) {
   	 String assignStr = cmd.getOptionValue("a");   
   	 if(assignStr == "SVD" || assignStr == "1")
   		 assignType = SVD;
     }
    
    /* factors */
    String[] facWeightStr = cmd.getOptionValues("f");
    double[] lambda = new double[facWeightStr.length];
    System.out.printf("factor weight:");
    for(int i=0; i<facWeightStr.length; i++) {
      lambda[i] = Double.parseDouble(facWeightStr[i]);
      System.out.printf(" %f", lambda[i]);
    }
    System.out.println("\n");
    
    /* dimension */
    int dim = Integer.parseInt(cmd.getOptionValue("r"));
    System.out.printf("common dimension R: %d\n\n", dim);
    
    // if verbose
    int maxiter=100;
    if(cmd.hasOption("n"))
    	maxiter = Integer.parseInt(cmd.getOptionValue("n"));
    System.out.printf("max iteration: %d\n\n", maxiter);
    
    long initE = System.currentTimeMillis();
    System.out.printf("[INFO] initialized in %d ms.\n", initE - initS);   
    System.out.println("[INFO] begin training ...");
    
    /* train */
    int targetKey = 0;
    if(cmd.hasOption("t"))
    	targetKey = Integer.parseInt(cmd.getOptionValue("t"));
    System.out.printf("target domain: %d\n\n", targetKey);
    
    System.out.println("");
    
    int sFlag = 0;
    if(cmd.hasOption("H"))
    	sFlag = 1;
    
    int[] opt = {maxiter, assignType, sFlag};
    
    Matrix predValue = cdtf.run(X,dim,fIndicator,opt,convCrit,lambda,wX,EM,targetKey);    
    
    if(cmd.hasOption("T")) {
    	 Matrix trainDat;
    	    if(!EM)
    	    	trainDat = X[targetKey];
    	    else
    	    	trainDat = X[targetKey].copy();
    	
    	Matrix testDat = new Matrix(trainDat.getRowDimension(), trainDat.getColumnDimension(), 0.0);
    	
        Evaluator e = new Evaluator();
        if(cmd.hasOption("l")) {
          e.setIfLog(true);
          e.setLogfile(cmd.getOptionValue("l"));
        }

        String testStr = cmd.getOptionValue("T");
        int indexNum = Utils.readCSVMatrix(testStr, testDat, biFlag, 1);
        
        if(biFlag) { // implicit
            int[] atKs;
            if(cmd.hasOption("k")) {
              String[] kStr = cmd.getOptionValues("k");
              atKs = new int[kStr.length];
              for(int i=0; i<kStr.length; i++)
                atKs[i] = Integer.parseInt(kStr[i]);
            }
            else {
              atKs = new int[]{5, 10, 15, 20, 50, 80, 100, 200};
            }

            e.avgPrecisionAndRecallAtK(
                predValue,
                trainDat, 
                testDat, 
                atKs);
            
            System.out.println();
            for(int K : atKs) {
              Evaluator.PR pr = e.getPRatK(K);
              System.out.printf("top[%3d]: hit = %d\tprecTot = %d\trecTot = %d\tavgPrec = %f\tavgRec = %f\n",
                  K, pr.hit, pr.precTot, pr.recTot, pr.avgPrecision, pr.avgRecall);
            }
          }
          else {  // explicit
            
	    	  int[][] testIdx = new int[indexNum][2];
	          Utils.readCSVIndex(testStr, testIdx);
            
            e.rmseAndMae(testIdx, predValue, testDat);
            System.out.printf("rmse:\t%f\nmae:\t%f\n", e.getRmse(), e.getMae());
          }
        
      }
    
    long toc = System.currentTimeMillis();
	
	System.out.println("Elapsed time: " + (toc-tic)*1.0/1000 + "s");	//calculate the time
	
	if(cmd.hasOption("o"))
	{
		System.setOut(out);
		System.out.println("Finished! Please check the output result on file!");
	}
  }
  
  
  /**
	 * Redirect the output to file
	 * @param J The common dimension of factors, which needs to be adjusted 
	 */
	public static void redirect(String Path)
	{
		PrintStream ps=null;
		try {
			ps = new PrintStream(Path);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.setOut(ps);
	}
}
