import java.io.*;
import java.sql.Time;
import java.util.Date;
import java.util.StringTokenizer;
import java.math.*;

public class  clustering{
////////////////////////////////////////////////////////////
	//ok New Akhare mordad
	final static double C = 1;
	final static double U = 0.75;
	final static double D = 0.5;
	final static double R = 0.25;

	// coefficients:
	final static double RT = 0.3;  // the condition for choosing step 2 or 3
	final static double Ts = 0.01; //for Pr - Pr2<Ts
	final static double Tr = 0;    //for Z-Z2<Tr
	
	/////////////////////////       Data Structures        /////////////////////////////////////
	
	public class Info   // a structure to store considered rows to sort and append by maximums
    {
		 public int ID;
		 public double Z;
		 public Info(){this.ID = 0; this.Z = 0;}
    };
    
   public class ServNode // Service Node : to store service information
    {
    	
    	public int Srow;
    	public int Scol;
    	public int Frow;
    	public int Fcol;
    	public int ID;
    	public ServNode Next;
  
    ServNode(ServNode node,int SR,int FR,int SC,int FC,int I)
        {
    	Srow =SR;  //Start  row
    	Scol =SC;  // Start  col
    	Frow =FR;  // Finish row
    	Fcol = FC; // finish col
    	ID = I;    // ID of service
    	Next = node; // point to next service
    	}
    };
  
	
	//
	////////////////////     Functions     //////////////////////////
	public static void QuickSort(double[] Acc,int[] ID,int left,int right)  // to sort considered rows
	{
	for(int i = 0;i<right;i++)
		for(int j=0;j<right;j++)
		 if(Acc[i]<Acc[j])
		 {
			 double t;
		 	 t= Acc[i];
			 Acc[i] = Acc[j];
			 Acc[j]= t;
			 int g;
			 g= ID[i];
			 ID[i] = ID[j];
			 ID[j] = g;
		 }
	}
	
	//////////////////////////////////////////////////////////////////
	public static double pow(double a , double b) // : a^b
	{
		double orgA=a;
		for(;b>1;b--)
			a *=orgA;
		return a;
	}
	////////////////////////////////////////////////////////////////
	public static void ColExch(double [][] matrix,int rows,int col1, int col2)  // exchange the column a( col1 ) with column b (col2)
	{
		double tmp;
		for(int i=0;i<=rows;i++)
		{
			tmp =matrix[i][col1];
			matrix[i][col1] = matrix[i][col2];
			matrix[i][col2] = tmp;
		}
	}
	////////////////////////////////////////////////////////////////
	public static void RowExch(double [][] matrix,int Cols,int row1, int row2) // exchange the row a( row1 ) with row b (row2)
	{
		double tmp;
		for(int i=0;i<=Cols+1;i++)
		{
			tmp =matrix[row1][i];
			matrix[row1][i] = matrix[row2][i];
			matrix[row2][i] = tmp;
		}
	}
	
		
	/////////////////////          MAIN           /////////////////////////////
	
	public static void main(String arg[]){
		double GT=0,ZT=0,AffT=0,coupT=0;
		int IDCNT = 1;
		ServNode  head =null;
		clustering X = new clustering();
		double [][] matrix;
		boolean condition = false;
		double orgG=0,orgcoup=0;
		int CrowCNT = 0;  // hold Number of rows that has C 
		String FileName = "in.txt";		// address of input matrix

		
		try {
		      
			 BufferedReader br = new BufferedReader(new FileReader(FileName));

		      try {
		    	  
		    	  //////////// READING MATRIX FROM FILE //////////////////
		    	  String line = null;
		    	  int row = 0,col=0;
		    	  line = br.readLine();
		    	  
		    	  
		    	  row = Integer.parseInt(line);
		    	  line = br.readLine();
		    	  col = Integer.parseInt(line);
		    	  matrix = new double [row+1][col+2];
		    	  System.out.println("matrix order is: " + row +": "+ col);
		    	  System.out.println();
		    	  System.in.read();
		    	  //Set ID for each row and column
		    	  for(int i=0;i<=row;i++)
		    		  matrix[i][col+1] = i;
		    	  for(int i=0;i<col;i++)
		    		  matrix[row][i] = i;
		    	  
		    	  System.out.println("   the matrix elements: ");
		    	  line = br.readLine();
		    	  int rowIndex = 0,colIndex = 0;
		    	  
		    	  while(line != null)
		    	  {
		    	  StringTokenizer theLine = new StringTokenizer(line);
		    	  colIndex = 0;
		    	  int HasC =0;
		    	  while(theLine.hasMoreTokens()){
		    		  matrix [rowIndex][colIndex] = Double.parseDouble(theLine.nextToken());//.trim();
		    		 
		    		  if(matrix [rowIndex][colIndex] ==1)
		    			  HasC =1; 
		    		  System.out.print(matrix [rowIndex][colIndex]+" ");
		    		  colIndex = colIndex + 1;
		    		  
		    		  }//end a row
		    	  CrowCNT += HasC;
		    	  matrix [rowIndex][col] = HasC;  // matrix[X][col] tell us that this row has C or not
		    	  System.out.println( "**" );
		    		  rowIndex = rowIndex + 1;
		    	  
		    		  
		    	  line = br.readLine();
		    	  }//end while
		    	  
		    	  //////////// END READING MATRIX FROM FILE //////////////////
		    	  
		    	  double Aff = 0,coup = 0,G=0,Z=0,Pr=0,tmp = 0;
		    	  System.in.read();
		    	  System.in.read();
		    	  ///////////////////// Starting Processing ///////////////////
		    	  rowIndex = colIndex = 0;
		    	  int Bcol=0,Brow=0;  // row and col Bound in this cluster
		    	 
		    	  int Sbrow = 0, Sbcol = 0;
		    	  
		    	  while(rowIndex <row-1)// end when there is no row to consider
		    	  {
		    		 	
		    		  Sbrow  = Brow;
		    		  Sbcol = Bcol;
		    		  for(int i = Brow;i<row;i++) // bring the Cs of a row together
		    		  {
		    			  if(matrix [i][col]!=C)
		    				  continue;
		    			 
		    			  for(int j =Bcol; j<col;j++)
		    				if(matrix[i][j] == C )
		    				{ 
		    					ColExch(matrix,row,j,Bcol);
		    					Bcol++;
		    				}
		    			  
		    			  System.out.println("first row in cluster ID ********  : " + (matrix[i][col+1]+1));
		    		  if(Bcol-Sbcol > 0)
		    		  {
		    			  RowExch(matrix, col, i, Brow);
  					  Brow++;
		    			  break;
		    		  }
		    		  }
		    		   rowIndex++;
		    		  
		    		  // Initial testing :
		    		  
		    		  
		    		Aff = 0;
		    		  for(int i= Sbcol;i<Bcol;i++) // calculate Aff
		    			  Aff += matrix[Sbrow][i];
		    		  
		    		  G = Aff/(Bcol-Sbcol);
		    		  tmp=0;
		    		  for(int i= Sbcol;i<col;i++)
		    			  tmp += matrix[Sbrow][i];
		    		 
		    		  coup = 1+tmp - Aff;
		    		  Z =( pow(G,3) * pow(Aff ,2))/pow(coup,2) ;
		    		
		    		  
		    		 double sum = 0;
		    		 int NonC = 0;
		    		 tmp=0;
		    			  	for(int i= Brow;i<row;i++)
		    			  	{
		    			  		int HasR;//**
		    			  		HasR = 0;//**
		    			  		for(int j= Sbcol;j<Bcol;j++)
		    			  		{	tmp += matrix[i][j];
		    			  			sum ++;
		    			  		if(matrix[i][j] == R)//**
		    			  			HasR=1;//**
		    			  		}
		    			  		if(matrix[i][col] ==0 && HasR ==1)//**
		    			  			NonC++;//**
		    			  	}
		    			  	
		    			  //** : for checking RT condition
		    			if(Brow==row)
		    				Pr = 0;
		    			else
		    	            Pr = 1 - tmp/sum;
		    	    
                     
		    	     // check RT condition.
		    	     if(NonC/sum <RT)
		    	     { 
////// Step 2  /////////////// consider C semantic R.///////////////////
		    	    	 System.out.println("Step 2 first with Z: " + Z);
		    	    	// System.out.println("Brow is :"+Brow + " Sbrow: "+ Sbrow);	
						    
		    	    	 {
		    	         double Aff2 = 0,coup2 = 0,G2=0,Z2=0,Pr2=0,tmp2 = 0;
			    	     int AccCnt=0;
			    	     int Bcol2 = Bcol, Brow2=Brow; //Hold original bound of cluster
			    	     int i;
			    	   
			    	     double[] Acc = new double[CrowCNT];
				    	   int[] ID = new int[CrowCNT]; 

				    	  // System.out.println("Brow is :"+Brow +" Sbrow: "+ Sbrow);
				    	 for( i = Sbrow+1; i<row;i++)
		    	     {
			    		// System.out.println("cheking : " + (matrix[i][col+1]+1));
				    	    
			    		 Brow = Brow2;
			    		 Bcol = Bcol2;
		    	    	 if(matrix[i][col] ==0)
		    	    		 continue;
		    	    	 
		    	    	
		    	    	
	    				for(int j =Bcol; j<col;j++) // bring all Cs of a row together
			    				if(matrix[i][j] == C )
			    				{ 
			    					ColExch(matrix,row,j,Bcol);
			    					Bcol++;
			    					
			    				} 
	    				
	    				RowExch(matrix, col, i, Brow);
	    					Brow++;
	    				
	    				//	System.out.println("Bcol" +i+" :"+Bcol +"Brow: " +Brow);
	    				//////////// calculating Metrics ///////
	    				Aff2=0;
	    				for( int l= Sbrow;l<Brow;l++)
	    				for( int k= Sbcol;k<Bcol;k++)
			    			  Aff2 += matrix[l][k];
	    				
			    		  G2 = Aff2/(Bcol-Sbcol);
			    		 tmp2=0;
			    		  for( int l= Sbrow;l<Brow;l++)
			    			  for( int k=Sbcol;k<col;k++)
			    			  tmp2 += matrix[l][k];
			    		 
			    		  coup2 = 1+tmp2 - Aff2;
			    		  Z2 =( pow(G2,3) * pow(Aff2 ,2))/pow(coup2,2) ;
			    		
			    		 double sum2 = 0;
			    		 tmp2=0;
			    		 for( int l= Brow;l<row;l++)
			    			for(int j= Sbcol;j<Bcol;j++)
			    			{
			    				tmp2 += matrix[l][j];
			    			    sum2 ++;
			    			  		
			    			}
			    			if(Brow==row)
			    				Pr2 = Pr +1;
			    			else
			    	           Pr2 = 1 - tmp2/sum2;
	    				//////////////////////////////////////
	 
				    	  
	
      //ACTIVE    	     
			    	// System.out.println( i+ ":" +"coup: " +coup2 +"* Aff= " + Aff2 + " * G= "+ G2 + " *Pr= " +Pr2 + " * Z="  +Z2  );	
			    	// System.in.read();
			    
	
			    	     if(Pr - Pr2<Ts && Z-Z2<Tr) //accepted 
			    	     {
			    	    	 //System.out.println("accepted : " + (matrix[Brow-1][col+1]+1));
			    	    	 ID[AccCnt] = (int)matrix[Brow-1][col+1]; 
			    	    	 Acc[AccCnt]= Z-Z2;
			    	    	 AccCnt++;
			    	      }
		    	   
		    	     }//end of checking all C semantic
			    	
			    	 
				    	 Brow = Brow2;
			    		 Bcol = Bcol2;
			    	  Pr2 = Pr;
			    	 Z2 = Z;
			    	 
			    	 if(AccCnt>0)
			    	 {
			    	 // Now sort all accepted rows;
			    	 QuickSort(Acc,ID,0,AccCnt);
			    	
			    	 
			            //Add to cluster  start with maximum Z//////////////////////
			    	 int r = 0;
			    	 
			    	
			    //	 System.out.println("Bcol is " +Bcol);
			    	 double orgZ=0,orgPr=0;
			    	 int colPlus=0;
			    	 
			    	 
			    	 do
			    	  {
			    		 colPlus=0;
			    		// for(int v=0;v<AccCnt;v++)
			    			 //System.out.println((ID[v]+1));
			    		 orgG = G2;
			    		 orgcoup = coup2;
		    			 rowIndex++;
		    			orgZ = Z;
		    			orgPr = Pr;
		    			 Z = Z2;
		    			 Pr = Pr2;
		    			 G = G2;
			            Aff2 = 0;coup2 = 0;G2=0;Z2=0;Pr2=0;tmp2 = 0;
			            for(int y=0;y<row;y++)
			            	if((int)matrix[y][col+1] == ID[r] )
			            		{
			            		//System.out.println("ID[r]: "+(ID[r]+1)+ "   y: " +y + "   Brow: "+Brow);
			            		RowExch(matrix, col, y, Brow);
			            		break;
			            		}
			            
	    					
		    				
			            for(int j =Bcol; j<col;j++) // bring all Cs of a row together
				    				if(matrix[Brow][j] == C )
				    				{ 
				    					ColExch(matrix,row,j,Bcol);
				    					Bcol++;
				    					colPlus++;
				    					
				    				} 
		    				
			            Brow++;
    					r++;
		    			//	System.out.println("Now R is : " + r + " Bcol :"+Bcol +"Brow: " +Brow);
		    				//////////// calculating Metrics ///////
		    				Aff2=0;
		    				for( int l= Sbrow;l<Brow;l++)
		    				for( int k= Sbcol;k<Bcol;k++)
				    			  Aff2 += matrix[l][k];
		    				
				    		  G2 = Aff2/(Bcol-Sbcol);
				    		 tmp2=0;
				    		  for( int l= Sbrow;l<Brow;l++)
				    			  for( int k=Sbcol;k<col;k++)
				    			  tmp2 += matrix[l][k];
				    		 
				    		  coup2 = 1+tmp2 - Aff2;
				    		  Z2 =( pow(G2,3) * pow(Aff2 ,2))/pow(coup2,2) ;
				    		
				    		 double sum2 = 0;
				    		 tmp2=0;
				    		 for( int l= Brow;l<row;l++)
				    			for(int j= Sbcol;j<Bcol;j++)
				    			{
				    				tmp2 += matrix[l][j];
				    			    sum2 ++;
				    			  		
				    			}
				    		 if(Brow==row)
				    				Pr2 = Pr +1;
				    			else  	
				    	     Pr2 = 1 - tmp2/sum2;
		    				//System.out.println("checking: "+(ID[r-1]+1) + "Bcol is :  " + (Bcol - Sbcol));
			    	
		    				 	 
		   //System.out.println( "Pr-Pr2: "+ (Pr-Pr2) + "Z-Z2 :" + (Z-Z2) + " Acc: "+AccCnt);
			    	if(r ==1 && r<AccCnt)
			    		condition = true;
			    	else
			    		condition = Pr - Pr2<Ts && Z-Z2<Tr && r<AccCnt;
			    	  }while( condition);
			    	
			    	
			    	 
			    	 if(!(Pr - Pr2<Ts && Z-Z2<Tr) )
			    	 {
			    	 Brow--;
		    	     Bcol =Bcol- colPlus;
		    	    rowIndex--;
		    	    
			    	 }
			    	 Z= orgZ;
		    	     Pr = orgPr;
		    	     G= orgG ;
		    		 coup = orgcoup;
			    	 }
			    	 //System.out.println( "**** coup: " +coup +"* Aff= " + Aff + " * G= "+ G + " *Pr= " +Pr + " * Z="  +Z  );	
			    //	 System.out.println("Brow is :"+Brow + " Sbrow: "+ Sbrow + " RowIndex: " + rowIndex);	
					      
			    	 
		    	     }//end of block 2

//////// Step 3 /////////// consider non-C semantic R. 
		    	     System.out.println("Step 3 Started , Bcol is :" + (Bcol - Sbcol));
		    	    	
			    	     double Aff2 = 0,coup2 = 0,G2=0,Z2=0,Pr2=0,tmp2 = 0;
			    	     int AccCnt=0;
			    	     int Bcol2 = Bcol, Brow2=Brow; //Hold original bound of cluster
			    	     int i;
			    	   
			    	   double[] Acc = new double[row-CrowCNT];
			    	   int[] ID = new int[row-CrowCNT];
			    	//   System.out.println("Brow is :"+Brow + " Sbrow: "+ Sbrow);	
				    	
			    	//  System.out.println("Brow is :"+Brow);
			    	 for( i = Sbrow+1; i<row;i++)
		    	     {
			    		 Brow = Brow2;
			    		 Bcol = Bcol2;
		    	    	 if(matrix[i][col] ==1)//has c in row
		    	    		 continue;
		    	    	
	    				RowExch(matrix, col, i, Brow);
	    					Brow++;
	    			//	System.out.println("Bcol" +i+" :"+Bcol +"Brow: " +Brow);
	    				//////////// calculating Metrics ///////
	    				Aff2=0;
	    				for( int l= Sbrow;l<Brow;l++)
	    				for( int k= Sbcol;k<Bcol;k++)
			    			  Aff2 += matrix[l][k];
	    				
			    		  G2 = Aff2/(Bcol-Sbcol);
			    		 tmp2=0;
			    		  for( int l= Sbrow;l<Brow;l++)
			    			  for( int k=Sbcol;k<col;k++)
			    			  tmp2 += matrix[l][k];
			    		 
			    		  coup2 = 1+tmp2 - Aff2;
			    		  Z2 =( pow(G2,3) * pow(Aff2 ,2))/pow(coup2,2) ;
			    		
			    		 double sum2 = 0;
			    		 tmp2=0;
			    		 for( int l= Brow;l<row;l++)
			    			for(int j= Sbcol;j<Bcol;j++)
			    			{
			    				tmp2 += matrix[l][j];
			    			    sum2 ++;
			    			  		
			    			}
			    			  	
			    	     Pr2 = 1 - tmp2/sum2;
	    				//////////////////////////////////////
			    	    	
			    	    
			    	     if(Pr - Pr2<Ts && Z-Z2<Tr) //accepted 
			    	     {
			    	    	 ID[AccCnt] = (int)matrix[Brow-1][col+1]; 
			    	    	 Acc[AccCnt]= Z-Z2;
			    	    	 AccCnt++;
			    	      }
		    	   
		    	     }//end of checking all nonC semantic
			    	// System.out.println("ACCCnt : " + AccCnt);
			    	
			    	 Brow = Brow2;
		    		 Bcol = Bcol2;
			    	Pr2 = Pr;
			    	Z2 = Z;
			    	
			    	 
			    	if(AccCnt>0)
			    	{
			    	 // Now sort all accepted rows;
			    	 QuickSort(Acc,ID,0,AccCnt);
			    	
			    	
			    	 
			    	  //Add to cluster ://////////////////////
			    	 int r = 0;
			    	double orgZ = 0, orgPr = 0;
			    	do
			    	  {
			    	
			    		orgG = G2;
			    		 orgcoup = coup2;
		    			orgZ = Z;
		    			orgPr = Pr;
		    			
		    			 rowIndex ++;
		    			 Z = Z2;
		    			 Pr = Pr2;
		    			 Aff = Aff2;
			            Aff2 = 0;coup2 = 0;G2=0;Z2=0;Pr2=0;tmp2 = 0;
			            for(int y=0;y<row;y++)
			            	if((int)matrix[y][col+1] == ID[r] )
			            		{
			            	//	System.out.println("ID[r]: "+(ID[r]+1)+ "   y: " +y + "   Brow: "+Brow);
				            	
			            		RowExch(matrix, col, y, Brow);
			            		break;
			            		}
			            
	    					Brow++;
	    					r++;
		    				
			           
		    				
		    		//	System.out.println("Bcol" +i+" :"+Bcol +"Brow: " +Brow);
		    				//////////// calculating Metrics ///////
		    				Aff2=0;
		    				for( int l= Sbrow;l<Brow;l++)
		    				for( int k= Sbcol;k<Bcol;k++)
				    			  Aff2 += matrix[l][k];
		    				
				    		  G2 = Aff2/(Bcol-Sbcol);
				    		 tmp2=0;
				    		  for( int l= Sbrow;l<Brow;l++)
				    		  for( int k=Sbcol;k<col;k++)
				    			  tmp2 += matrix[l][k];
				    		 
				    		  coup2 = 1+tmp2 - Aff2;
				    		  Z2 =( pow(G2,3) * pow(Aff2 ,2))/pow(coup2,2) ;
				    		
				    		 double sum2 = 0;
				    		 tmp2=0;
				    		 for( int l= Brow;l<row;l++)
				    			for(int j= Sbcol;j<Bcol;j++)
				    			{
				    				tmp2 += matrix[l][j];
				    			    sum2 ++;
				    			  		
				    			}
				    			if(Brow ==row)
				    				Pr2 = Pr+1;
				    			else
				    	     Pr2 = 1 - tmp2/sum2;
				    	 
				 // 	 System.out.println( "coup: " +coup2 +"* Aff= " + Aff2 + " * G= "+ G2 + " *Pr= " +Pr2 + " * Z="  +Z2  );	
				  //  	System.out.println( "Pr-Pr2: "+ (Pr-Pr2) + "Z-Z2 :" + (Z-Z2) + " Acc: "+AccCnt);
				  //	 System.in.read();
				    	  // System.out.println("checking : " + (ID[r-1]+1));
				    	 if(r ==1 && r<AccCnt)
					    		condition = true;
					    	else
					    		condition = Pr - Pr2<Ts && Z-Z2<Tr && r<AccCnt;
					    	
				    	// System.out.println("rowIndex is : " + rowIndex);
					    	
			    	  }while( condition);
			    	
			    	    
			    	
			    	 if(!(Pr - Pr2<Ts && Z-Z2<Tr) )
			    	{
			    	 Brow--;
			    	rowIndex--;
			    	 Z= orgZ;
		    	     Pr = orgPr;
		    	     G=  orgG ;
		    		 coup = orgcoup;
		    		 }
			    	
			    	
			    	}
			    	    
			    	}//end if
		    	        
		    	     
		    	    else//first 3 and then 2
			    	{
//////// Step 3 /////////// consider non-C semantic R. 
		    	    	{  	   
		    	    		System.out.println("first Step 3 Started   **** Z is : "+ Z);
		    	    		//System.out.println("rowIndex is : " + rowIndex);
					    	
		    				    	     double Aff2 = 0,coup2 = 0,G2=0,Z2=0,Pr2=0,tmp2 = 0;
		    				    	     int AccCnt=0;
		    				    	     int Bcol2 = Bcol, Brow2=Brow; //Hold original bound of cluster
		    				    	     int i;
		    				    	   
		    				    	   double[] Acc = new double[row-CrowCNT];
		    				    	   int[] ID = new int[row-CrowCNT];

		    				    	//   System.out.println("Brow is :"+Brow);
		    				    	 for( i = Sbrow+1; i<row;i++)
		    			    	     {
		    				    			//System.out.println("Checking: "+(matrix[i][col+1]+1));
		    				            	
		    				    		 Brow = Brow2;
		    				    		 Bcol = Bcol2;
		    			    	    	 if(matrix[i][col] ==1)//has c in row
		    			    	    		 continue;
		    			    	    	
		    		    				RowExch(matrix, col, i, Brow);
		    		    					Brow++;
		    		    			//	System.out.println("Bcol" +i+" :"+Bcol +"Brow: " +Brow);
		    		    				//////////// calculating Metrics ///////
		    		    				Aff2=0;
		    		    				for( int l= Sbrow;l<Brow;l++)
		    		    				for( int k= Sbcol;k<Bcol;k++)
		    				    			  Aff2 += matrix[l][k];
		    		    				
		    				    		  G2 = Aff2/(Bcol-Sbcol);
		    				    		 tmp2=0;
		    				    		  for( int l= Sbrow;l<Brow;l++)
		    				    		   for( int k=Sbcol;k<col;k++)
		    				    			  	  tmp2 += matrix[l][k];
		    				    			 
		    				    		  
		    				     //System.out.println("tmp2 is : " +tmp2 + " ID is" + (matrix[i][col+1]+1) );
		    				    		  coup2 = 1+tmp2 - Aff2;
		    				    		  Z2 =( pow(G2,3) * pow(Aff2 ,2))/pow(coup2,2) ;
		    				    		
		    				    		 double sum2 = 0;
		    				    		 tmp2=0;
	    				    	
		    				    		// System.out.println("Bcol :"+Bcol +"Brow: " +Brow + "Sbrow: "+ Sbrow+ "Sbcol : " + Sbcol);
				    		    			
		    				    		
		    				    		 if(Brow==row)Brow--;
		    				    		 for( int l= Brow;l<row;l++)
		    				    			for(int j= Sbcol;j<Bcol;j++)
		    				    			{
		    				    				tmp2 += matrix[l][j];
		    				    			    sum2 ++;
		    				    			   
		    				    			}
		    				    			  Brow++;	
		    				    	     Pr2 = 1 - tmp2/sum2;
		    		    				//////////////////////////////////////
		    				    	//     System.out.println("sum2 : " + tmp2);
		    				     //    System.out.println( "coup: " +coup2 +"* Aff= " + Aff2 + " * G= "+ G2 + " *Pr= " +Pr2 + " * Z="  +Z2  );	
		    					 //  	 System.out.println( "Pr-Pr2: "+ (Pr-Pr2) + "Z-Z2 :" + (Z-Z2) + " Acc: "+AccCnt);
		    					 //  	 System.in.read();
		    						     
		    				    	    
		    				    	     if(Pr - Pr2<Ts && Z-Z2<Tr) //accepted 
		    				    	     {
		    				    	    	 ID[AccCnt] = (int)matrix[Brow-1][col+1]; 
		    				    	    	 Acc[AccCnt]= Z-Z2;
		    				    	    	 AccCnt++;
		    				    	      }
		    			    	   
		    			    	     }//end of checking all nonC semantic
		    				    	 Brow = Brow2;
		    			    		 Bcol = Bcol2;
		    			    		 Pr2 = Pr;
		    			    		 Z2 = Z;
		    				    	if(AccCnt>0)
		    				    	{
		    				    	 // Now sort all accepted rows;
		    				    	 QuickSort(Acc,ID,0,AccCnt);
		    				    	
		    				    
		    				    	// for(int v=0;v<AccCnt;v++)
		    			    			 //System.out.println((ID[v]+1));
		    				    	  //Add to cluster ://////////////////////
		    				    	 int r = 0;
		    				    	
		    			    		double orgZ = 0, orgPr = 0;
		    				    	do
		    				    	  {
		    				    		orgG = G2;
		    				    		 orgcoup = coup2;
		    			    			orgZ = Z;
		    			    			orgPr = Pr;
		    			    			 rowIndex ++;
		    			    			 Z = Z2;
		    			    			 Pr = Pr2;
		    			    			 Aff = Aff2;
		    				            Aff2 = 0;coup2 = 0;G2=0;Z2=0;Pr2=0;tmp2 = 0;
		    				            for(int y=0;y<row;y++)
		    				            	if((int)matrix[y][col+1] == ID[r] )
		    				            		{
		    				            		//System.out.println("ID[r]: "+(ID[r]+1)+ "   y: " +y + "   Brow: "+Brow);
		    					            	
		    				            		RowExch(matrix, col, y, Brow);
		    				            		break;
		    				            		}
		    				            
		    		    					Brow++;
		    		    					r++;
		    			    				
		    				           
		    			    				
		    			    				//System.out.println("Bcol" +i+" :"+Bcol +"Brow: " +Brow);
		    			    				//////////// calculating Metrics ///////
		    			    				Aff2=0;
		    			    				for( int l= Sbrow;l<Brow;l++)
		    			    				for( int k= Sbcol;k<Bcol;k++)
		    					    			  Aff2 += matrix[l][k];
		    			    				
		    					    		  G2 = Aff2/(Bcol-Sbcol);
		    					    		 tmp2=0;
		    					    		  for( int l= Sbrow;l<Brow;l++)
		    					    		  for( int k=Sbcol;k<col;k++)
		    					    			  tmp2 += matrix[l][k];
		    					    		 
		    					    		  coup2 = 1+tmp2 - Aff2;
		    					    		  Z2 =( pow(G2,3) * pow(Aff2 ,2))/pow(coup2,2) ;
		    					    		
		    					    		 double sum2 = 0;
		    					    		 tmp2=0;
		    					    		 for( int l= Brow;l<row;l++)
		    					    			for(int j= Sbcol;j<Bcol;j++)
		    					    			{
		    					    				tmp2 += matrix[l][j];
		    					    			    sum2 ++;
		    					    			  		
		    					    			}
		    					    			if(Brow ==row)
		    					    				Pr2 = Pr+1;
		    					    			else 	
		    					    	     Pr2 = 1 - tmp2/sum2;
		    					    	     //System.out.println("on adding a row : "+ (ID[(r-1)]+1));
		    					    	// System.out.println( "coup: " +coup2 +"* Aff= " + Aff2 + " * G= "+ G2 + " *Pr= " +Pr2 + " * Z="  +Z2  );	
		    					    	 //System.out.println( "Pr-Pr2: "+ (Pr-Pr2) + "Z-Z2 :" + (Z-Z2) + " Acc: "+AccCnt);
		    					    	 //System.in.read();
		    						     
		    					    	 if(r ==1 && r<AccCnt)
		    						    		condition = true;
		    						    	else
		    						    		condition = Pr - Pr2<Ts && Z-Z2<Tr && r<AccCnt;
		    						    	 
		    					        	
		    				    	  }while( condition);
		    				    	
		    				    	
		    				    	if(!(Pr - Pr2<Ts && Z-Z2<Tr) )
		    				    	{
		    				    	 Brow--;
		    				    	 rowIndex--;
		    				    	Z= orgZ;
		    			    	     Pr = orgPr;
		    			    	     G=  orgG ;
		    			    		 coup = orgcoup;
		    			    		 }
		    				    	 
		    				    	}
			    	}
////// Step 2  /////////////// consider C semantic R.///////////////////
		    	    	System.out.println("Step 2 started");
		    	    	
		    	    	double Aff2 = 0,coup2 = 0,G2=0,Z2=0,Pr2=0,tmp2 = 0;
			    	     int AccCnt=0;
			    	     int Bcol2 = Bcol, Brow2=Brow; //Hold original bound of cluster
			    	     int i;
			    	   
			    	     double[] Acc = new double[CrowCNT];
				    	   int[] ID = new int[CrowCNT]; 
				    	   //System.out.println("Brow is  : "+ Sbrow);
			    	 for( i = Sbrow+1; i<row;i++)
		    	     {
			    		 Brow = Brow2;
			    		 Bcol = Bcol2;
			    		 //System.out.println("cheking : " + (matrix[i][col+1]+1));
			    	    	
			    		 if(matrix[i][col] ==0)
		    	    		 continue;
		    	    	 
		    	    	
		    	    	
	    				for(int j =Bcol; j<col;j++) // bring all Cs of a row together
			    				if(matrix[i][j] == C )
			    				{ 
			    					ColExch(matrix,row,j,Bcol);
			    					Bcol++;
			    					
			    				} 
	    				
	    				RowExch(matrix, col, i, Brow);
	    					Brow++;
	    				//System.out.println("Bcol" +i+" :"+Bcol +"Brow: " +Brow);
	    				//////////// calculating Metrics ///////
	    				Aff2=0;
	    				for( int l= Sbrow;l<Brow;l++)
	    				for( int k= Sbcol;k<Bcol;k++)
			    			  Aff2 += matrix[l][k];
	    				
			    		  G2 = Aff2/(Bcol-Sbcol);
			    		 tmp2=0;
			    		  for( int l= Sbrow;l<Brow;l++)
			    		  for( int k= Sbcol;k<col;k++)
			    			  tmp2 += matrix[l][k];
			    		 
			    		  coup2 = 1+tmp2 - Aff2;
			    	//	  System.out.println("G3 : " + pow(G2,3) + "  pow(Aff2,2): " +pow(Aff2,2));
		                
			    		  Z2 =( pow(G2,3) * pow(Aff2 ,2))/pow(coup2,2) ;
			    		
			    		 double sum2 = 0;
			    		 tmp2=0;
			    		 for( int l= Brow;l<row;l++)
			    			for(int j= Sbcol;j<Bcol;j++)
			    			{
			    				tmp2 += matrix[l][j];
			    			    sum2 ++;
			    			  		
			    			}
			    			  	
			    	     Pr2 = 1 - tmp2/sum2;
	    				//////////////////////////////////////
	 
				    	  
	
     //ACTIVE    	     
			    	// System.out.println( "coup: " +coup2 +"* Aff= " + Aff2 + " * G= "+ G2 + " *Pr= " +Pr2 + " * Z="  +Z2  );	
			    //	 System.out.println( "Pr-Pr2: "+ (Pr-Pr2) + "Z-Z2 :" + (Z-Z2) + " Acc: "+AccCnt);
			    //	 System.in.read();
				     
			    	
			    
	
			    	     if(Pr - Pr2<Ts && Z-Z2<Tr) //accepted 
			    	     {
			    	    	 //System.out.println("accepted : " + (matrix[Brow-1][col+1]+1));
				    	    	
			    	    	 ID[AccCnt] = (int)matrix[Brow-1][col+1]; 
			    	    	 Acc[AccCnt]= Z-Z2;
			    	    	 AccCnt++;
			    	      }
		    	   
		    	     }//end of checking all C semantic
			    	
			    	 
			    	 
			    	 
			    	 Brow = Brow2;
		    		 Bcol = Bcol2;
			    	Pr2 = Pr;
			    	 Z2 = Z;
			    	 if(AccCnt>0)
			    	 {
			    	 // Now sort all accepted rows;
			    	 QuickSort(Acc,ID,0,AccCnt);
			    	
			            //Add to cluster  start with maximum Z//////////////////////
			    	 int r = 0;
			    	  
			    	// System.out.println("Bcol is " +Bcol);
			    	 
			    	 int colplus=0;
			    	double orgZ= 0,orgPr = 0;
			    	
			    	do
			    	  {
			    		colplus=0;
			    	//	for(int v=0;v<AccCnt;v++)
			    			 //System.out.println((ID[v]+1));
			    		orgG = G2;
			    		orgcoup = coup2;
			    		 orgZ = Z;
			    		 orgPr = Pr;
			    		 
		    			 rowIndex++;
		    			 Z = Z2;
		    			 Pr = Pr2;
		    			 G = G2;
			            Aff2 = 0;coup2 = 0;G2=0;Z2=0;Pr2=0;tmp2 = 0;
			            for(int y=0;y<row;y++)
			            	if((int)matrix[y][col+1] == ID[r] )
			            		{
			            		//System.out.println("ID[r]: "+(ID[r]+1)+ "   y: " +y + "   Brow: "+Brow);
				            	
			            		RowExch(matrix, col, y, Brow);
			            		break;
			            		}
			            
	    					
		    				
			            for(int j =Bcol; j<col;j++) // bring all Cs of a row together
				    				if(matrix[Brow][j] == C )
				    				{ 
				    					ColExch(matrix,row,j,Bcol);
				    					Bcol++;
				    					colplus++;
				    					
				    				} 
		    				
			            Brow++;
   					r++;
		    		//		System.out.println("Now R is : " + r + " Bcol :"+Bcol +"Brow: " +Brow);
		    				//////////// calculating Metrics ///////
		    				Aff2=0;
		    				for( int l= Sbrow;l<Brow;l++)
		    				for( int k= Sbcol;k<Bcol;k++)
				    			  Aff2 += matrix[l][k];
		    				
				    		  G2 = Aff2/(Bcol-Sbcol);
				    		 tmp2=0;
				    		  for( int l= Sbrow;l<Brow;l++)
				    		  for( int k= Sbcol;k<col;k++)
				    			  tmp2 += matrix[l][k];
				    		 
				    		  coup2 = 1+tmp2 - Aff2;
				    		 
				    		  Z2 =( pow(G2,3) * pow(Aff2 ,2))/pow(coup2,2) ;
				    		
				    		 double sum2 = 0;
				    		 tmp2=0;
				    		 for( int l= Brow;l<row;l++)
				    			for(int j= Sbcol;j<Bcol;j++)
				    			{
				    				tmp2 += matrix[l][j];
				    			    sum2 ++;
				    			  		
				    			}
				    			if(Brow ==row)
				    				Pr2 = Pr+1;
				    			else
				    	     Pr2 = 1 - tmp2/sum2;
				    	     //System.out.println("checking : " + (ID[r-1]+1));
			    	       
			    	 if(r ==1 && r<AccCnt)
				    		condition = true;
				    	else
				    		condition = Pr - Pr2<Ts && Z-Z2<Tr && r<AccCnt;
				    	
			    	  }while( condition);
			    			    	
			    	
			    	if(!(Pr - Pr2<Ts && Z-Z2<Tr) )
			    	{
			    	Brow--;
			    	 rowIndex--;
			    	 Bcol =Bcol- colplus;
		    	     
			    	 
			    	 G= orgG;
		    	     coup = orgcoup;
		    	     Z= orgZ;
		    	     Pr = orgPr;
			    	}
			    	 }
			    	}//end else
			    	
		    	   
			    	
		    	     
		    	     
			    	
			    	
			    	// Add this service and it's bound to Export List
		    		head= X.new ServNode(head,Sbrow,Brow,Sbcol,Bcol,IDCNT);
		    		Aff =  (G*(Bcol-Sbcol));
		    		GT += G;
		    		AffT +=Aff;
		    		coupT += coup;
		    		
		    		System.out.println();
		   	System.out.println("Aff: " + Aff + "  G: "+ G + "     coup:" + coup);
		    		double ZT1 = (pow(G,3) * pow(Aff,2))/pow(coup,2);
		    		
		    		for(int o = Sbrow ;o<Brow;o++)
		    			System.out.println("ID "+(int)(matrix[o][col+1]+1)+" is in this service");
		    		
		    		System.out.println("  Bounded Rows is : " + (Brow-Sbrow)+ "     Bounded Cols is: " + (Bcol-Sbcol));
		    		System.out.println("********** Service "+ IDCNT + "***** ZT: " + ZT1+"  **************");
		    	    IDCNT++;
		    	    System.in.read();
		    	    System.in.read();
		    	  }//end main while
		
		    	 ZT = (pow(GT,3) * pow(AffT,2))/pow(coupT,2);
		    	 
		    	 
		    	 
		    	 System.out.println(" **********************  ToTal Z :" + ZT + " ****************");
		    
		    	 System.in.read();
		    	 System.in.read();
		    	
		    	 
		    	 
		    	 //// Now We can write the Services list to .XML file , HERE:////
		    	  //
		    	  //
		    	  //
		    	  ////////////////////////////////////////////////////////////////
		    	  
		      } //end TRY
		     finally {
		        br.close();
		      }
				}
		    catch (IOException ex){
		      ex.printStackTrace();
		    }
	
 }//main
}//clustering
