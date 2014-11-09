package edu.sjtu.cit;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import edu.umbc.cs.maple.utils.MathUtils;
import Jama.Matrix;



public class Evaluator {

  // {{{
  private double mae  = 0.0;
  private double rmse = 0.0;
  private HashMap<Integer, PR> precAndRec = new HashMap<Integer, PR>();
  private boolean ifLog = false;
  private String logfile = "log.txt";
  // }}}
  
  // {{{
  /**
   * Helper class to handle precision and recall, with hit/precTotal/recTotal included.
   * 
   * @author lythesia
   *
   */
  public class PR {
    public int hit = 0;
    public int precTot = 0;
    public int recTot = 0;
    public double avgPrecision = 0.0;
    public double avgRecall = 0.0;
  }
  
  /**
   * Helper class to wrap rating and item-id for sorting.
   * 
   * @author lythesia
   *
   */
  class Score implements Comparable<Score> {

    public double rating;
    public int    id;
    
    public Score(double rating, int id) {
      this.rating = rating;
      this.id = id;
    }

    /* for descend sorting */
    public int compareTo(Score s) {
      if(rating < s.rating) return 1;
      else if(rating > s.rating) return -1;
      else return 0;
    }    
  }
  // }}}
  
  // {{{
  /**
   * Calculate RMSE and MAE with given range to cast predict rating.
   * 
   * @param testIdx index that represents positive rating in test matrix.
   * @param pred    predict matrix.
   * @param test    test matrix.
   */
  public void rmseAndMae(int[][] testIdx, Matrix pred, Matrix test, double l, double r) {
    int missNum = testIdx.length;
    mae = rmse = 0.0;
    for(int i=0; i<missNum; i++) {
      double predVal = pred.get(testIdx[i][0], testIdx[i][1]);
      if(predVal > r) // cast upperbound
        predVal = r;
      if(predVal < l) // cast lowerbound
        predVal = l;
      
      double diff = test.get(testIdx[i][0], testIdx[i][1]) - predVal;
      
      rmse += diff * diff;
      mae  += Math.abs(diff);
    }
    rmse = Math.sqrt(rmse / missNum);
    mae /= missNum;
  }
  
  /**
   * Calculate RMSE and MAE.
   * 
   * @param testIdx index that represents positive rating in test matrix.
   * @param pred    predict matrix.
   * @param test    test matrix.
   */
  public void rmseAndMae(int[][] testIdx, Matrix pred, Matrix test) {
    rmseAndMae(testIdx, pred, test, 0, 1);
  }
  
  /**
   * Get top(K) of predict matrix to recommend with given K and user-id.
   * 
   * @param pred  predict matrix.
   * @param train train matrix.
   * @param test  test matrix.
   * @param uid   user id.
   * @param K     top-K.
   * @return      List<item-id>.
   */
  public ArrayList<Integer> topK(Matrix pred, Matrix train, Matrix test, int uid, int K) {
    int len = pred.getColumnDimension();
    int R = 0;
    
    Score[] rank = new Score[len];
    for(int i=0; i<len; i++) {
      double predVal = pred.get(uid, i);
      if(train.get(uid, i) != 0.0)
        predVal = 0.0; // already in train, do not recommend it
      else
        R++;
      rank[i] = new Score(predVal, i);
    }
    Arrays.sort(rank);
    
    R = Math.min(R, K); // if K is to large, there'll be less than K recommend items.
    ArrayList<Integer> top = new ArrayList<Integer>(R);
    for(int i=0; i<R; i++) {
      top.add(rank[i].id);
    }
    
    return top;
  }
  
  /**
   * Calculate average precision and recall among all users.
   * 
   * @param pred  predict matrix.
   * @param train train matrix.
   * @param test  test matrix.
   * @param atK   list of top-Ks.
   * @throws IOException 
   */
  public void avgPrecisionAndRecallAtK(Matrix pred, Matrix train, Matrix test, int[] atK) throws IOException {
    int nUser = pred.getRowDimension(),
        nItem = pred.getColumnDimension();
    int maxK = MathUtils.maxValue(atK);
    
    // init
    for(int K : atK)
      precAndRec.put(K, new PR());    
    
    int nzNum = 0; // number of valid users, who rated at least one item in test matrix.
    for(int i=0; i<nUser; i++) {
      boolean nzflag = false;
      for(int j=0; j<nItem; j++) {
        if(test.get(i, j) != 0) {
          nzflag = true;
          break;
        }
      }
      if(nzflag) { // valid user
        
        nzNum++;
                
        ArrayList<Integer> topk = topK(pred, train, test, i, maxK); // top-N of current user, N set to nItem to fit different top-K for 1-pass computing
        
        for(int K : atK) {          
          double tmpPrec = 0.0,
                 tmpRec  = 0.0;
          int    tmpHit  = 0,
                 tmpPrecTot = 0,
                 tmpRecTot = 0;
          
          tmpPrecTot = Math.min(K, topk.size());
          
          for(int j=0; j<nItem; j++) {
            if(test.get(i, j) != 0) {
              tmpRecTot++;
              int idx = topk.indexOf(j);
              if(idx != -1 && idx < tmpPrecTot)
                tmpHit++;
            }
          }
          tmpPrec = 1.0 * tmpHit / tmpPrecTot;
          tmpRec  = 1.0 * tmpHit / tmpRecTot;
          
          // logging
          if(ifLog) {
            String fileAtK = "K" + K + "_" + logfile + ".txt";
            PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileAtK, true)), true);
            writer.printf("%d, %d, %d, %d, %f, %f\n", 
                i, tmpHit, tmpPrecTot, tmpRecTot, tmpPrec, tmpRec);
            writer.flush();
            writer.close();
          }

          PR prAtK = precAndRec.get(K);
          
          prAtK.hit += tmpHit;
          prAtK.precTot += tmpPrecTot;
          prAtK.recTot += tmpRecTot;
        }        
      }      
    }
    
    // average precision/recall
    if(nzNum > 0) {
      for(int K : atK) {
        PR prAtK = precAndRec.get(K);
        prAtK.avgPrecision = 1.0 * prAtK.hit / prAtK.precTot;
        prAtK.avgRecall    = 1.0 * prAtK.hit / prAtK.recTot;
      }
    }
  }
  
  /**
   * Get precision and recall with given top-K.
   * 
   * @param K
   * @return
   */
  public PR getPRatK(int K) {
    return precAndRec.get(K);
  }
  // }}}
  
  public void setIfLog(boolean ifLog) {
    this.ifLog = ifLog;
  }
  
  public double getMae() {
    return mae;
  }

  public void setMae(double mae) {
    this.mae = mae;
  }

  public double getRmse() {
    return rmse;
  }

  public void setRmse(double rmse) {
    this.rmse = rmse;
  }

  public HashMap<Integer, PR> getPrecAndRec() {
    return precAndRec;
  }

  public void setPrecAndRec(HashMap<Integer, PR> precAndRec) {
    this.precAndRec = precAndRec;
  }

  public String getLogfile() {
    return logfile;
  }

  public void setLogfile(String logfile) {
    this.logfile = logfile;
  }
}
