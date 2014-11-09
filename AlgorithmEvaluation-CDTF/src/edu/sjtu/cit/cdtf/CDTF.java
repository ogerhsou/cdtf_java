package edu.sjtu.cit.cdtf;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import edu.umbc.cs.maple.utils.JamaUtils;
import Jama.*;

public class CDTF {
	static int RANDOM = 0;	//To initialize factors by RANDOM
	static int SVD = 1;	//To initialize factors by SVD
	static double eps = 2.2204e-16;	//very small number 
	
	Matrix A, B, C;	//the three factors of CDTF
	
	//indicator matrix whose entry has been rated (for observed values) 
	//and zero otherwise (for missing values)
	Matrix M;
	
	Matrix XtX;	//square matrix X'X
	Matrix predValues;	//save the reconstructive matrix
	
	Matrix[] P;	//the column-wise orthonormal matrix
	Matrix[] X;	//store the origin matrices array
	Matrix[] missingOnes;	//matrices array to store the missing indices
	
	double[] wX;	//store the weight of each matrix
	double[] lambda;	//store each factor's regularization
	
	int J;	//the common dimension of three factors
	int I;	//the row of main matrix
	int K;	//the number of input matrices
	int maxIt = 100;	//inner maximum iteration time
	int initi;	//initial model (for further use)
	int missNum, allNum;	//store the existing value of test set, the other store the all 
	
	boolean[] factorIndicator;	//To indicate the factor constant or not
	boolean showOut;	//indication of showing output or not
	boolean testFlag=false;	//use to test with matlab data
	boolean EM;	//EM is true means it will use EM algorithm to update X[k]
	boolean initFacMatrix;	//judge if factors have been initialized or not
	
	double absErr = 1e-6;	//to judge if it is convergent(fit < absErr)
	double convCrit = 1e-3;	//to compare with the tolerance 
	double tol;	//store the tolerance
	
	/**
	 * The main function of CDTF
	 * @param initFacMatrix	The flag to see if it has been initialized
	 * @param X	Input matrices array
	 * @param R	The common dimension of three factors
	 * @param cons	Constraint for further use
	 * @param options	Receiving parameters
	 * @param lambda	Regularization of each factor
	 * @param wX	Weight of each factor
	 * @param EM missingData is true means it will use EM algorithm to update X[k]
	 */
	public Matrix run(Matrix[] X, int R, boolean[] factorIndicator, int[] options, double conv, double[] lambda, double[] wX, boolean EM, int targetKey)
	{
		//===========Initialize the values===========================
		this.X = X;
		this.EM = EM;
		this.lambda = lambda;
		this.factorIndicator = factorIndicator;
		this.wX = wX;
		if(conv!=0)
			this.convCrit = conv;
		if(options[0]!=0)
		{
			maxIt = options[0];
		}
		initi = options[1];
		showOut = options[2]==0;
		
		I = X[0].getRowDimension();
		J = R;
		K = X.length;
		
		if(showOut)
		{
			System.out.println("Convergence criterion        : " + convCrit);
			System.out.println("Maximal number of iterations : " + maxIt);
			System.out.println("Common size of factors       : " + R);
		}
		
		if(!initFacMatrix)
			P = new Matrix[K];
		
		//========================================================
		
		//scale slabs
		for(int k=0;k<K;++k)
		{
			X[k].timesEquals(wX[k]);
			if(initFacMatrix)
			{
				for(int j=0;j<C.getColumnDimension();++j)
				{
					C.set(k, j, C.get(k, j)*wX[k]);
				}
			}
		}
		
		//===============update X[k] by EM algorithm==========================
		
		//find missing ones and replace with average	
		//if missingData is true, X[k] will be updated by EM algorithm, else, it will not
		if(EM)
		{
			if(!initFacMatrix)	//only need to be executed at the first time
			{
				missingOnes = new Matrix[K];
				missNum=0;allNum=0;
				for(int k=0;k<K;++k)
				{
					Matrix x = X[k];
					missingOnes[k] = new Matrix(x.getRowDimension(), x.getColumnDimension());
					double sum=0, cnt=0;
					
					int xrow = x.getRowDimension(), xcol = x.getColumnDimension();
					for(int i=0;i<xrow;++i)
					{
						for(int j=0;j<xcol;++j)
						{
							if(x.get(i, j)!=0)	//if the coordinate has value 
							{
								sum+=x.get(i, j);	//store the sum of the matrix
								++cnt;	//store the number of existing value of training set
							}
							else	//if there is no value
							{
								missingOnes[k].set(i, j, 1);	//store the sub
								++missNum;	//count the missing number, whose value is 0
							}
						}
					}
					
					double mean = sum * 1.0 / cnt;	//calculate the mean value of existing values on this matrix
					
					for(int i=0;i<xrow;++i)
					{
						for(int j=0;j<xcol;++j)
						{
							if(x.get(i, j)==0)
							{
								x.set(i, j, mean);	//set the missing value to be the mean value
							}
						}
					}
						
					X[k] = x;	//update it
					allNum += xrow * xcol;	//calculate all the combinations of all the matrices
				}
				
				if(showOut)
				{
					double percMiss = missNum * 1.0 / allNum;	//the percentage of missing values of training set  
					System.out.println("Miss: " + missNum + " All: " + allNum);
					System.out.println("Missing data handled by EM: " + percMiss);
				}
			}
			else	//if initialized before
			{
				for(int k=0;k<K;++k)
				{
					Matrix x = X[k];
					int xrow = x.getRowDimension(), xcol = x.getColumnDimension();
					Matrix tmpDiag = new Matrix(C.getColumnDimension(), C.getColumnDimension());
					for(int i=0;i<tmpDiag.getColumnDimension();++i)
					{
						tmpDiag.set(i, i, C.get(k, i));	//change the kth row of C to be a diagnal matrix
					}
					M = A.times(tmpDiag).times(P[k].times(B).transpose());	//calculate M matrix
					
					for(int i=0;i<xrow;++i)
					{
						for(int j=0;j<xcol;++j)
						{
							if(x.get(i, j)==0)
							{
								x.set(i, j, M.get(i, j));	//replace the missing value with M
							}
						}
					}
					X[k] = x;	//update it
				}
			}
		}
		
		//==============initialize factors====================
		double fit = 10000000;	//original fit: set it to be large but not infinity;
		if(!initFacMatrix)	//if the factors have not been initialized
		{	
			if(initi==SVD)	//use SVD to initialize the factors
			{
				//get the sum of (User X User) square matrix 
				XtX = X[0].times(X[0].transpose());
				for(int k=1;k<K;++k)
				{
					XtX = XtX.plus(X[k].times(X[k].transpose()));	
				}
				
				System.out.println("SVD based initialization.");
				SingularValueDecomposition svd = new SingularValueDecomposition(XtX);	//SVD of XtX
				A = svd.getU();	//the left matrix of SVD XtX
				
				A = A.getMatrix(0,A.getRowDimension()-1,0,R-1);	//extract first F columns to be A
				C = Matrix.random(K, R).timesEquals(0.1).plus(JamaUtils.ones(K, R));	//initialize C with random value
				B = Matrix.identity(R, R);		//B is identity matrix of F X F
				fit = XtX.trace();	//the fit is the trace of XtX
			}
			else if(initi==RANDOM)	//initialize factors with random values
			{
				System.out.println("Random based initialization.");
				A = Matrix.random(I, R);
				C = Matrix.random(K, R);
				B = Matrix.identity(R, R);
			}
		}
		
		//==================calculate fit (remove to save time)===============
//		if(initFacMatrix || initi!=1)
//		{
//			//get the sum of (User X User) square matrix 
//			XtX = X[0].times(X[0].transpose());
//			for(int k=1;k<K;++k)
//			{
//				XtX = XtX.plus(X[k].times(X[k].transpose()));	
//			}
//		}
		//===================================================================
		
		System.out.println("Fitting model ...");
		
		tol = Double.POSITIVE_INFINITY;
		int it=0;
		while(tol>convCrit && it<maxIt && fit>absErr)	
		{
			++it;
			double oldfit = fit;
			
			Matrix[] Y = new Matrix[K];
			
			//update P
			for(int k=0;k<K;++k)
			{
				Matrix tmpdiagC = new Matrix(C.getColumnDimension(), C.getColumnDimension());
				for(int tra=0;tra<tmpdiagC.getColumnDimension();++tra)
				{
					tmpdiagC.set(tra, tra, C.get(k, tra));
				}
				
				Matrix Qk = X[k].transpose().times(A.times(tmpdiagC).times(B.transpose()));	//Qk = X{k}'*(A*diag(C(k,:))*B');		
				
				P[k] = Qk.times(psqrt(Qk.transpose().times(Qk)));	//P{k} = Qk*psqrt(Qk'*Qk);
				Y[k] = X[k].times(P[k]);	//Y(:,:,k) = X{k}*P{k};
			}
			
			Matrix reshapeY = new Matrix(I,J*K);	//dim: I X J*K
			
			for(int k=0;k<K;++k)
			{
				reshapeY.setMatrix(0, I-1, J*k, J*(k+1)-1, Y[k]);	//tile matrices into one
			}
			
			int[] dimY = {I,J,K};
			int innerIt = 1;	//only iterate once
			parafac(reshapeY, dimY, R, convCrit, factorIndicator, innerIt);	//call function to update factors
			
			fit = pf2fit();	//calculate fit
			
			tol = Math.abs(fit-oldfit)/oldfit;

			if(showOut)
			{
				System.out.println("FitErr: " + fit + " Iter: " + it + " Tol: " + tol + " Oldfit: " + oldfit);
			}

		}
		
		//scale slabs
		for(int k=0;k<K;++k)
		{
			if(wX[k]!=1)
			{
				for(int j=0;j<C.getColumnDimension();++j)
				{
					C.set(k, j, C.get(k, j)/wX[k]);
				}
			}
		}
		
		reconstruct(targetKey);	//reconstruct the matrix
		return predValues;
	}
	
	
	/**
	 * Fits the PARAFAC model Xk = A*Dk*B.' + E where Dk is a diagonal matrix holding the k'th row of C.
	 * @param X	Data
	 * @param dimX	Dimension of X
	 * @param fac	Number of factors
	 * @param crit	Convergence criterion (default 1e-6)
	 * @param factorIndicator	[a b c], if e.g. a=false => A need to be updated, a=true => A constant
	 * @param maxIt	Maximum iteration times
	 */
	public void parafac(Matrix X, int[] dimX, int fac, double crit, boolean[] factorIndicator, int maxIt)
	{
		//get the regularization of each factor
		Matrix regularA, regularB, regularC;
		regularA = Matrix.identity(fac, fac).times(lambda[0]);
		regularB = Matrix.identity(fac, fac).times(lambda[1]);
		regularC = Matrix.identity(fac, fac).times(lambda[2]);
		
		int I = dimX[0], J = dimX[1], K = dimX[2];
		
		boolean constA = factorIndicator[0], constB = factorIndicator[1], constC = factorIndicator[2];
		
		//calculate the original fit(the power sum of each elements)
		//====================remove it to save time====================
//		double sumSqX=0;
//		for(int i=0;i<X.getRowDimension();++i)
//		{
//			for(int j=0;j<X.getColumnDimension();++j)
//			{
//				sumSqX += Math.pow(X.get(i, j),2);
//			}
//		}
//		double fit = sumSqX;
//		double fitold = 2*fit;
		
//		while(Math.abs((fit-fitold)/fitold) > crit && --maxIt>=0 && fit > 10*eps)
		while(--maxIt>=0)
		{
//			fitold = fit;
			Matrix Xbc = new Matrix(I,J);
			
			//update A
			for(int k=0;k<K;++k)
			{
				int[] cols = new int [J];
				for(int i=0;i<J;++i)
				{
					cols[i] = i+k*J;
				}
				Matrix tmpDiag = new Matrix(C.getColumnDimension(),C.getColumnDimension());
				for(int i=0;i<tmpDiag.getRowDimension();++i)
				{
					tmpDiag.set(i, i, C.get(k, i));
				}
				//Xbc = Xbc + X(:,(k-1)*J+1:k*J)*conj(B*diag(sparse(C(k,:))))
				Xbc = Xbc.plus(JamaUtils.getcolumns(X, cols).times(B.times(tmpDiag)));
			}
			
			if(!constA)
			{
				//A = Xbc*pinv((B'*B).*(C'*C)+regularA).'
				A = Xbc.times(pinv(B.transpose().times(B).arrayTimes(C.transpose().times(C)).plus(regularA)).transpose());
			}
			
			//QR decomposition of A to update B
			QRDecomposition qr = new QRDecomposition(A);
			Matrix Qa = qr.getQ();
			Matrix Ra = qr.getR();
			Matrix x = Qa.transpose().times(X);
			
			Matrix Xac = new Matrix(J,J);
			
			//update B
			for(int k=0;k<K;++k)
			{
				int[] cols = new int [J];
				for(int i=0;i<J;++i)
				{
					cols[i] = i+k*J;
				}
				Matrix tmpDiag = new Matrix(C.getColumnDimension(),C.getColumnDimension());
				for(int i=0;i<tmpDiag.getRowDimension();++i)
				{
					tmpDiag.set(i, i, C.get(k, i));
				}
				//Xac = Xac + x(:,(k-1)*J+1:k*J).'*conj(Ra*diag(sparse(C(k,:))))
				Xac = Xac.plus(JamaUtils.getcolumns(x, cols).transpose().times(Ra.times(tmpDiag)));
			}	
			
			if(!constB)
			{
				//B = Xac*pinv((Ra'*Ra).*(C'*C)+regularB).'
				B = Xac.times(pinv(Ra.transpose().times(Ra).arrayTimes(C.transpose().times(C)).plus(regularB)).transpose());
			}
				
			//update C
			if(!constC)
			{
				//ab=pinv((Ra'*Ra).*(B'*B)+regularC)
				Matrix ab = pinv(Ra.transpose().times(Ra).arrayTimes(B.transpose().times(B)).plus(regularC));
				for(int k=0;k<K;++k)
				{
					int[] cols = new int [J];
					for(int i=0;i<J;++i)
					{
						cols[i] = i+k*J;
					}
					Matrix tmpDiag = new Matrix(C.getColumnDimension(),1);
					Matrix tmpM = Ra.transpose().times(JamaUtils.getcolumns(x, cols).times(B));
					for(int i=0;i<tmpDiag.getRowDimension();++i)
					{
						tmpDiag.set(i, 0, tmpM.get(i, i));
					}
					//C(k,:) = (ab*diag(Ra'* x(:,(k-1)*J+1:k*J)*conj(B))).'
					JamaUtils.setrow(C, k, ab.times(tmpDiag).transpose());
				}
				
			}
			
			//TODO
			//QR of jama cannot be used
//			QRDecomposition qrb = new QRDecomposition(B);
//			Matrix Qb = qrb.getQ();
//			Matrix Rb = qrb.getR();
//			QRDecomposition qrc = new QRDecomposition(C.transpose());
//			Matrix Z = qrc.getQ();
//			Matrix Rc = qrc.getR();
//			
//			Matrix tmpppp = ppp(Rb,Rc).transpose();
//			tmpppp = Ra.times(tmpppp);
//			
//			double tmpsum=0;
//			for(int i=0;i<tmpppp.getRowDimension();++i)
//			{
//				for(int j=0;j<tmpppp.getColumnDimension();++j)
//				{
//					tmpsum+=Math.pow(tmpppp.get(i, j),2);
//				}
//			}
//			
//			fit = sumSqX - tmpsum;	
			
			//=========calculate fit (remove to save time)====================
//			fit = 0;
//			for(int k=0;k<K;++k)
//			{
//				Matrix residual = new Matrix(X.getRowDimension(), J);
//				int[] cols = new int [J];
//				for(int i=0;i<J;++i)
//				{
//					cols[i] = i+k*J;
//				}
//				
//				Matrix tmpDiag = new Matrix(C.getColumnDimension(),C.getColumnDimension());
//				for(int i=0;i<tmpDiag.getRowDimension();++i)
//				{
//					tmpDiag.set(i, i, C.get(k, i));
//				}
//				//residual=X(:,(k-1)*J+1:k*J)-A*diag(C(k,:))*B.'
//				residual = JamaUtils.getcolumns(X, cols).minus(A.times(tmpDiag).times(B.transpose()));
//				
//				for(int i=0;i<residual.getRowDimension();++i)
//				{
//					for(int j=0;j<residual.getColumnDimension();++j)
//					{
//						//fit=fit+sum(sum((abs(residual).^2)))
//						fit += Math.pow(residual.get(i, j), 2);
//					}
//				}
//
//			}
			//=================================================================
		}
		
		//ORDER ACCORDING TO VARIANCE
		double[] tuck = new double[J];
		Matrix tmp = A.transpose().times(A).arrayTimes(B.transpose().times(B)).arrayTimes(C.transpose().times(C));
		for(int i=0;i<tmp.getRowDimension();++i)
		{
			//Tuck = diag((A'*A).*(B'*B).*(C'*C))
			tuck[i] = tmp.get(i, i);
		}
		double[] out = tuck.clone();
		int[] ID = new int[J];
		for(int i=0;i<ID.length;++i)
		{
			ID[i] = i;
		}
		bubbleSort(out,ID);
		A = JamaUtils.getcolumns(A, ID);
		B = JamaUtils.getcolumns(B, ID);
		C = JamaUtils.getcolumns(C, ID);
		
		//NORMALIZE A AND C (variance in B)
		
		Matrix normA = new Matrix(fac,1);
		Matrix normC = new Matrix(fac,1);
		
		for(int f=0;f<fac;++f)
		{
			normC.set(f, 0, JamaUtils.getcol(C, f).norm2());
			normA.set(f, 0, JamaUtils.getcol(A, f).norm2());
		}
		
		Matrix diagNormA = new Matrix(fac,fac);
		Matrix diagNormC = new Matrix(fac,fac);
		Matrix diagNorm_A = new Matrix(fac,fac);
		Matrix diagNorm_C = new Matrix(fac,fac);
		
		for(int i=0;i<fac;++i)
		{
			diagNormA.set(i, i, normA.get(i, 0));
			diagNormC.set(i, i, normC.get(i, 0));
			diagNorm_A.set(i, i, 1.0 / normA.get(i, 0));
			diagNorm_C.set(i, i, 1.0 / normC.get(i, 0));
		}
		
		B = B.times(diagNormC).times(diagNormA);
		A = A.times(diagNorm_A);
		C = C.times(diagNorm_C);
		
		//APPLY SIGN CONVENTION
		Matrix signA_tmp = new Matrix(A.getRowDimension(), A.getColumnDimension());
		
		int tmprow = signA_tmp.getRowDimension(), tmpcol = signA_tmp.getColumnDimension();
		for(int i=0;i<tmprow;++i)
		{
			for(int j=0;j<tmpcol;++j)
			{
				signA_tmp.set(i, j, Math.signum(A.get(i, j)));
			}
		}
		Matrix signA = new Matrix(1,tmpcol);
		
		for(int i=0;i<tmpcol;++i)
		{
			signA.set(0, i, Math.signum(JamaUtils.colsum(signA_tmp, i)+eps));
		}
		
		Matrix signC_tmp = new Matrix(C.getRowDimension(), C.getColumnDimension());
		
		tmprow = signC_tmp.getRowDimension(); tmpcol = signC_tmp.getColumnDimension();
		for(int i=0;i<tmprow;++i)
		{
			for(int j=0;j<tmpcol;++j)
			{
				signC_tmp.set(i, j, Math.signum(C.get(i, j)));
			}
		}
		Matrix signC = new Matrix(1,tmpcol);
		
		for(int i=0;i<tmpcol;++i)
		{
			signC.set(0, i, Math.signum(JamaUtils.colsum(signC_tmp, i)+eps));
		}
		
		Matrix spDiagA = new Matrix(signA.getColumnDimension(), signA.getColumnDimension());
		Matrix spDiagC = new Matrix(signC.getColumnDimension(), signC.getColumnDimension());
		
		for(int i=0;i<spDiagA.getColumnDimension();++i)
		{
			spDiagA.set(i, i, signA.get(0, i));
		}
		
		for(int i=0;i<spDiagC.getColumnDimension();++i)
		{
			spDiagC.set(i, i, signC.get(0, i));
		}
		
		//Update A, B and C
		A = A.times(spDiagA);
		C = C.times(spDiagC);
		B = B.times(spDiagA).times(spDiagC);
		
	}
	
	/**
	 * Sort the arrays as well as the sub by bubbleSort
	 * @param number	Array to be sort
	 * @param sub	Original sub of array
	 */
	public void bubbleSort(double[] number, int[] sub)
	{
		double temp;
		int tmpsub;
		for(int i=0;i<number.length;i++)
			for(int j=0;j<number.length-1-i;j++)
				if(number[j]>number[j+1]) {
					temp=number[j];
					tmpsub=sub[j];
					number[j]=number[j+1];
					sub[j]=sub[j+1];
					number[j+1]=temp;
					sub[j+1]=tmpsub;
				} 
	}
	
	
	/**
	 * Produces A^(-.5) even if rank-problems
	 * @param A	Input matrix
	 * @return	A^(-.5)
	 */
	public Matrix psqrt(Matrix A)
	{
		//SVD A
		SingularValueDecomposition svd = new SingularValueDecomposition(A);
		Matrix U = svd.getU();
		Matrix s = svd.getS();
		Matrix V = svd.getV();
		
		Matrix X;
		Matrix S;
		
		//set diagnal matrix into vector 
		S = new Matrix(s.getRowDimension(), 1);
		for(int i=0;i<s.getRowDimension();++i)
		{
			S.set(i, 0, s.get(i, i));
		}
		
		tol = Math.max(A.getRowDimension(), A.getColumnDimension()) * S.get(0, 0) * eps;
		
		int r=0;
		for(int i=0;i<S.getRowDimension();++i)
		{
			if(S.get(i, 0) > tol)
			{
				++r;
			}
		}
		
		if(r==0)
			X = new Matrix(A.getColumnDimension(),A.getRowDimension());
		else
		{
			Matrix tmpS = new Matrix(r,r);
			for(int i=0;i<r;++i)
			{
				tmpS.set(i, i, 1.0/Math.sqrt(S.get(i, 0)));	//diag(ones(r,1)./sqrt(S(1:r)))
			}
			
			int[] cols = new int[r];
			for(int i=0;i<r;++i)
			{
				cols[i]=i;
			}
			//X = V(:,1:r)*S*U(:,1:r)'
			X = JamaUtils.getcolumns(V, cols).times(tmpS).times(JamaUtils.getcolumns(U, cols).transpose());
		}
		return X;
	}
	
	
	/**
	 * Pseudoinverse
	 * X = PINV(A) produces a matrix X of the same dimensions as A' so that A*X*A = A, 
	 * X*A*X = X and A*X and X*A are Hermitian
	 * @param A	Input Matrix
	 * @return	Pseudoinverse of A
	 */
	public Matrix pinv(Matrix A)
	{	
		Matrix X;
		int m = A.getRowDimension(), n = A.getColumnDimension();
		
		if(n>m)	//if n>m, transpose A
			X = pinv(A.transpose());
		else
		{
			//SVD A
			SingularValueDecomposition svd = new SingularValueDecomposition(A);
			Matrix U = svd.getU();
			Matrix S = svd.getS();
			Matrix V = svd.getV();
			
			double[] s = new double[m];
			double maxs=Double.NEGATIVE_INFINITY;
			for(int i=0;i<m;++i)
			{
				s[i] = S.get(i, i);
				if(s[i]>maxs)
					maxs = s[i];
			}
			
			tol = Math.max(m, n) * eps * Math.floor(maxs);
			
			int r=0;
			for(int i=0;i<m;++i)
			{
				if(s[i]>tol)
					++r;
			}
			if(r==0)
				X = new Matrix(A.getColumnDimension(),A.getRowDimension());
			else
			{
				//tmps = diag(ones(r,1)./s(1:r))
				Matrix tmps = new Matrix(r,r);
				for(int i=0;i<r;++i)
				{
					tmps.set(i, i, 1.0/s[i]);
				}
				int[] cols = new int[r];
				for(int i=0;i<r;++i)
					cols[i] = i;
				
				//X = V(:,1:r)*s*U(:,1:r)'
				X = JamaUtils.getcolumns(V, cols).times(tmps).times(JamaUtils.getcolumns(U, cols).transpose());
			}
		}
		return X;
	}
	
	/**
	 * The parallel proportional profiles product - triple-P product
	 * @param A	Matrix A
	 * @param B	Matrix B
	 * @return	Triple-P product
	 */
	public Matrix ppp(Matrix A, Matrix B)
	{
		int I = A.getRowDimension(), F = A.getColumnDimension();
		int J = B.getRowDimension(), F1 = B.getColumnDimension();
		
		if(F!=F1)
		{
			System.out.println("Error! Matrices must have the same number of columns!");
		}
		
		Matrix AB = new Matrix(I*J, F);
		
		for(int f=0;f<F;++f)
		{
			//ab=A(:,f)*B(:,f).'
			Matrix ab = JamaUtils.getcol(A, f).times(JamaUtils.getcol(A, f).transpose());
			//AB(:,f)=ab(:)
			JamaUtils.setcol(AB, f, ab);
		}
		return AB;
	}
	
	/**
	 * Calculate fit and impute missing elements from model
	 * @return	Fit
	 */
	public double pf2fit()
	{
		double fit = 0;
		for(int k=0;k<K;++k)
		{
			Matrix M;
			int Ccol = C.getColumnDimension();
			Matrix tmpDiag = new Matrix(Ccol,Ccol);
			for(int i=0;i<tmpDiag.getRowDimension();++i)
			{
				tmpDiag.set(i, i, C.get(k, i));
			}
			
			//M = A*diag(sparse(C(k,:)))*(P{k}*H)'
			M = A.times(tmpDiag).times(P[k].times(B).transpose());
						
			if(EM)
			{
				int row = missingOnes[k].getRowDimension(), col = missingOnes[k].getColumnDimension();
				Matrix x = X[k];
				
				for(int i=0;i<row;++i)
				{
					for(int j=0;j<col;++j)
					{
						if(missingOnes[k].get(i, j)==1)
						{
							x.set(i, j, M.get(i, j));
						}
					}
				}
				X[k] = x;
			}
			
			Matrix tmp = X[k].minus(M);
			int row = X[k].getRowDimension(), col = X[k].getColumnDimension();
			for(int i=0;i<row;++i)
			{
				for(int j=0;j<col;++j)
				{
					//fit = fit + sum(sum(abs (X{k} - M ).^2))
					fit += Math.pow(tmp.get(i, j),2);
				}
			}
		}
		return fit;
	}
	
	/**
	 * Calculate reconstructive matrix
	 */
	public void reconstruct(int targetKey)
	{
		Matrix tmpDiag = new Matrix (C.getColumnDimension(), C.getColumnDimension());
		for(int i=0;i<tmpDiag.getColumnDimension();++i)
		{
			tmpDiag.set(i, i, C.get(targetKey, i));
		}
		
		//A*diag(sparse(C(1,:)))*(P{1}*H)'
		predValues = A.times(tmpDiag).times(P[0].times(B).transpose());
	}

}
