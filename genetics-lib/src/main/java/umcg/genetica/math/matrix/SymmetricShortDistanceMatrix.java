/*
 * SymmetricShortDistancaMatrix.java
 *
 * Created on 04 June 2004, 16:03
 */

package umcg.genetica.math.matrix;

/**
 *
 * Symmetric Short Distance Matrix: A memory efficient matrix capable of performing calculations on all genes x all genes:
 *
 * @author  Lude Franke
 */
public class SymmetricShortDistanceMatrix {

    private int size;
    private short[] matrix;
    private final static double eLog10 = java.lang.Math.log(10);
    private final static int MAX_ALL_PAIRS_STEPS = 50;
    public final static int MAX_VALUE = Short.MAX_VALUE + 32768;
    private long[] elementIndex;
    
    /** Creates a new instance of SymmetricShortDistancaMatrix
     *
     * Defines a new Symmetric matrix, with a predefined size
     *
     * Initially all values will be Short.MAX_Value (65535)
     *
     * All data is stored in memory efficient one dimensional array, which costs size * (size + 1) bytes
     *
     */
    public SymmetricShortDistanceMatrix(int size) {
    
        this.size = size;
        long arraySize = ((long) size * (long) (size + 1)) / 2l;
        matrix = new short[(int) arraySize];
        elementIndex = new long[size]; for (int x=0; x<size; x++) elementIndex[x] = (long) x * (long) size - ((long) x * (long)(x+1)) / 2l;
        setMaxDistance();
        
    }
    
    public SymmetricShortDistanceMatrix(int size, boolean setMaxDistance) {
        
        this.size = size;
        long arraySize = ((long) size * (long) (size + 1)) / 2l;
        matrix = new short[(int) arraySize];
        elementIndex = new long[size]; for (int x=0; x<size; x++) elementIndex[x] = (long) x * (long) size - ((long) x * (long)(x+1)) / 2l;
        if (setMaxDistance) setMaxDistance();
    }

    public void setAllElements(int value) {
        for (int x=0; x<size; x++) {
            for (int y = x; y<size; y++) {
                matrix[(int) getElement(x,y)] = (short) (value - 32768);
            }
        }
    }
    
    private void setMaxDistance() {

        for (int x=0; x<size; x++) {
            for (int y = x; y<size; y++) {
                matrix[(int) getElement(x,y)] = Short.MAX_VALUE;
            }
        }
        
    }
    
    /* Get the array element location for point (x,y):
     */
    private long getElement(int x, int y) {
        if (x>y) {
            return elementIndex[y] + (long) x;
        } else {
            return elementIndex[x] + (long) y;
        }
    }
    
    /* Set the value for point (x,y):
     */
    public void set(int x, int y, int value) {
        matrix[(int) getElement(x,y)] = (short) (value - 32768);
    }
    
    /* Get the value for point (x,y):
     */
    public int get(int x, int y) {
        return matrix[(int) getElement(x,y)] + 32768;
    }
    
    /* Get size of matrix:
     */
    public int size() {
        return size;
    }
    
    public int maxValue() {
        return Short.MAX_VALUE + 32768;
    }
    
       
    /* Get the shortest path for all the pairs using the Floyd-Warshall algorithm: 
     *
     * Please be aware that this method replaces the values inside the current matrix with the shortest path values!
     */
    public void getAllPairsShortestPath() {

        //This method uses the Floyd-Warshall all pairs shortest path algorithm:
        //As we deal here with a symmetric matrix, the number of required calculations is 50% less.

        //Copy matrix:
        //System.out.println("Copying matrix to temporary array:");
        short[][] matrixCopy = new short[size()][size()];
        for (int x=0; x<size(); x++) {
            for (int y=x+1; y<size(); y++) {
                int value = get(x,y);
                if (value>32767) value = 32767;
                matrixCopy[x][y] = (short) (value);
                matrixCopy[y][x] = (short) (value);
            }
        }
        for (int v=0; v<size(); v++) {
            matrixCopy[v][v] = 0;
        }
        matrix = null;
        System.gc();

        /*
        int previousPercentage = 0;
        for (int i=0; i<size()&&i<MAX_ALL_PAIRS_STEPS; i++) {
            for (int v=0; v<size(); v++) {
                for (int w=v; w<size(); w++) {
                    int vi = get(v,i);
                    int iw = get(i,w);
                    if (vi!=maxValue()&&iw!=maxValue()) {
                        set(v, w, (short) java.lang.Math.min(get(v, w), vi + iw));
                    }
                }
            }
            int percentage = (int) (double) (100d * (double) i / (double) size());
            if (percentage!=previousPercentage) System.out.print(percentage + "% ");
            previousPercentage = percentage;
            System.out.print(".");
            
        }
        System.out.println("");
        */
        
        /*
        // Warm up the virtual machine, and time it
        short[][] matrixForHotspot = new short[500][500];
        long before = System.currentTimeMillis();
        performLoops(matrixForHotspot, (short) 500);
        long loopTime = System.currentTimeMillis() - before;
        //System.out.println("Loop time: " + Long.toString(loopTime) + " milliseconds");
         */
        
        //Perform actual all-pairs shortest path analysis:
        performLoops(matrixCopy, (short) size());

        //Copy the results back to the original matrix:
        matrix = new short[(size  * (size + 1)) / 2];
        for (int x=0; x<size(); x++) {
            for (int y=x+1; y<size(); y++) {
                set(x,y, matrixCopy[x][y]);
            }
        }
        
        //Garbage collect:
        matrixCopy = null;
        System.gc();
        
        //Set the total distance to genes itself to 0:
        //for (int v=0; v<size(); v++) {
            //set(v,v,0);
        //}
                
    }
    
    private void performLoops(short[][] matrixCopy, short sz) {
        //System.out.println("Performing actual all pairs shortest path calculation:");
        int previousPercentage = 0;
        int[] genes = new int[4];
        genes[0] = 925;
        genes[1] = 9077;
        genes[2] = 10207;
        genes[3] = 17440;
        /*
        genes[0] = 10;
        genes[1] = 200;
        genes[2] = 30;
        genes[3] = 40;
        matrixCopy[10][200] = 1;
        matrixCopy[200][10] = 1;
        matrixCopy[10][30] = 1;
        matrixCopy[30][10] = 1;
        matrixCopy[30][40] = 1;
        matrixCopy[40][30] = 1;
        matrixCopy[200][40] = 96;
        matrixCopy[40][200] = 96;
        matrixCopy[10][10] = 32000;
        matrixCopy[40][10] = 32000;
        matrixCopy[200][30] = 100;
        matrixCopy[30][200] = 100;
        */
        for (short k=0; k<sz; k++) {
            for (short i=0; i<sz; i++) {
                short matrix_ik = matrixCopy[i][k];
                for (short j=i; j<sz; j++) {
                    int value = (int) matrix_ik + (int) matrixCopy[k][j];
                    if (value < (int) matrixCopy[i][j]) {
                        short valueShort = (short) value;
                        matrixCopy[i][j] = valueShort;
                        matrixCopy[j][i] = valueShort;
                        matrix_ik = matrixCopy[i][k];
                    }
                }
            }
            //int percentage = (int) (double) (100d * (double) k / (double) size());
            //if (percentage!=previousPercentage) System.out.print(percentage + "% ");
            //previousPercentage = percentage;
            if (sz>1000 && k%10==9) System.out.println(k);
            /*
            if (k%10==9&&sz>20000) {
                System.out.print(k + " ");
            }
            if (k%100==99&&sz>20000) {
                System.out.println("");
                for (int x=0; x<4; x++) {
                    System.out.print("\t" + genes[x]);
                }
                System.out.println("");
                for (int x=0; x<4; x++) {
                    System.out.print(genes[x]);
                    for (int y=0; y<4; y++) {
                        System.out.print("\t" + matrixCopy[genes[x]][genes[y]]);
                    }
                    System.out.println("");
                }
            } 
             */
        }
        System.out.println("");
    }
    
    
    /* Get the shortest path from node X to node Y.
     *
     * Returns a Vector containing the nodes as Integers that in the correct order make up the path. Empty if no path exists.
     *
     * This algorithm uses the Dijkstra's shortest path algorithm.
     * 
     * Please be aware that it does not use a Fibonacci heap. Implementing this improves the performance considerably.
     */
    public java.util.Vector getShortestPath(int x, int y) {
        
        //This method uses the Dijkstra's shortest path algorithm:
        
        //Create path which eventually will contain all the nodes in the shortest path from x to y:
        java.util.Vector path =  new java.util.Vector();
        
        //Initialise S: Array indicating whether node is already included and traversed
        boolean S[] = new boolean[size()];
        
        //Initialise T: Two dimensional array, column 0: distance, column 1: nodeID predecessor:
        int T[][] = new int[2][size()];   
        for (int i=0; i<size(); i++) {
            S[i] = false;
            T[0][i] = maxValue();
            T[1][i] = -1;
        }
        
        //Starting from nodeID: nodeStart: set traversed true:
        S[x] = true;
        
        //Distance: distance from nodeStart to nodeStart is 0 as it is the beginning of path:
        T[0][x] = 0;
        
        //Set nodeStart, no predecessor, nodeStart = nodeStart:
        T[1][x] = x;
        
        //Fill T array with all directly connected nodes with node nodeStart:
        for (int i=0; i<size(); i++) {
            if (i!=x&&get(i, x)<maxValue()) {
                
                //Directly connected node, set distance:
                T[0][i] = get(i, x);
                
                //Set direct predecessor of this node:
                T[1][i] = x;
            }
        }
        
        //Check whether paths diverge from node nodeStart:
        boolean pathPossible = false;
        for (int i=0; i<size(); i++) {
            if (get(i, x)>0) {
                pathPossible = true;
            }
        }
        
        //Hypothesis: a path does exist between nodeStart and nodeEnd:
        boolean pathExists = true;

        //Traverse path, do not stop:
        boolean traversePath = true;
        
        //Perform Dijkstra algorithm loop, as long as a path is possible, exists and no stop request has been raised:
        while (pathPossible && pathExists && traversePath) {
            
            //Set minimal distance T:
            int tMin = maxValue();
            
            //Find node v that is nearest to nodeStart:
            int v = -1;
            
            //Traverse all nodes:
            for (int i=0; i<size(); i++) {
                
                //Find node whose shortest path neighbour has not been found and is directly connected:
                if (!S[i]&&T[0][i]<tMin) {
                    
                    //Node, directly connected and nearer has been observed:
                    tMin = T[0][i];
                    v = i;
                    
                }
                
            }

            //Check whether any connected node is found:
            if (v==-1) {
                pathExists = false;
                break;
            }
            
            //The nearest not finalized node v has been found, process this node:
            S[v] = true;
        
            //Traverse all nodes:
            for (int i=0; i<size(); i++) {
            
                //Find node whose final shortest path has not been established and which is directly connected:
                if (!S[i]&&get(v,i)>0) {
                    
                    //Check whether the distance of the edge between the processed node and node v is the lowest:
                    if (T[0][v] + get(v, i) < T[0][i]) {
                        
                        //Change the distance and nodeIndex of the i:
                        T[0][i] = T[0][v] + (int) get(v,i);
                        T[1][i] = v;
                        
                    }
                }
            }
            
            //Check whether destination has been reached:
            if (S[y]) { 
                traversePath = false;
            } else {
                //Check whether the path still between nodeStart and nodeEnd is still likely to exist:
                pathExists = false;
                for (int i=0; i<size(); i++) {
                    if (!S[i]&&T[0][i] < maxValue()) {
                        pathExists = true;
                        break;
                    }
                }
            }
        }
        
        //Check whether the path exists and is possible:
        if (pathExists&&pathPossible) {

            //Select nodes that have been involved, start with nodeEnd:
            path.add(new Integer(y));
            
            //Select predecessor nodes:
            int totalDistance = T[0][y];
            int predecessor = T[1][y];
            path.add(new Integer(predecessor));
            
            while (predecessor!=x) {
                totalDistance += T[0][predecessor];
                predecessor = T[1][predecessor];
                path.add(new Integer(predecessor));
            }

        }

        return path;
        
    }
    

    public int getDistance(double pValue) {
        return (int) (java.lang.Math.max(java.lang.Math.log(pValue) / eLog10, -7.9998779296875) * 8192 + 65535);
    }

    public double getPValue(int distance) {
        return (double) java.lang.Math.pow(10, ((double) distance - 65535) / 8192);
    }

    /** The distance measure in short ranges from 0 to 65.000. The distance in correlation (Pearson, Spearman, whatever)
     *  ranges from -1 to 1. This method translates correlationdistance into short distance. */
    public int getDistanceFromCorrelation(double correlation){
        int distance;
        if (correlation > 0){
            distance = (int)(-(double)Short.MIN_VALUE + ((correlation) * ((double)Short.MAX_VALUE)));  
        } else {
            distance = (int)(-(double)Short.MIN_VALUE + ((correlation) * (double)Short.MIN_VALUE));    
        }
        return distance;
    }
    
    public double getCorrelationFromDistance(int distance){
        double correlation;
        if (distance > -Short.MIN_VALUE){
            correlation = (((double)distance + (double) Short.MIN_VALUE)/(double)Short.MAX_VALUE);
        }else{
            correlation = (((double)distance + (double) Short.MIN_VALUE)/(-(double)Short.MIN_VALUE));
        }
        return correlation;
    }

    public void save(java.io.File fileName) {
        
        try {
            
            //Initialize output stream:
            java.io.OutputStream out = new java.io.BufferedOutputStream(new java.io.FileOutputStream(fileName));

            //Define the buffer size: 
            int bufferLength = 2500 * 1024;
            byte[] buffer = new byte[bufferLength * 2];

            //Define length of read bytes:
            int len = 0;

            //Initialize location where data is and should go:
            int loc = 0;

            while(1==1){

                // Copy data
                if (loc + bufferLength <= matrix.length) {
                    len = bufferLength;
                    int bufferLoc = 0;
                    for (int x=loc; x<len + loc; x++) {
                        buffer[bufferLoc] = (byte)(matrix[x] >> 8);
                        buffer[bufferLoc+1] = (byte)(matrix[x] & 0xff);
                        bufferLoc+=2;
                    }
                    loc+= len;
                    System.out.print(".");
                    out.write(buffer);
                } else {
                    len = matrix.length - loc;
                    buffer = new byte[len * 2];
                    int bufferLoc = 0;
                    for (int x=loc; x<len + loc; x++) {
                        buffer[bufferLoc] = (byte)(matrix[x] >> 8);
                        buffer[bufferLoc+1] = (byte)(matrix[x] & 0xff);
                        bufferLoc+=2;
                    }
                    out.write(buffer);
                    System.out.print(".");
                    break;
                }
            }

            System.out.println("");

            //Dereference buffer:
            buffer = null;

            //Close file:
            out.close();

        } catch (Exception e) {
            System.out.println("Cannot write to file! (" + e.getMessage() + ")");
        }
    }

    public void load(java.io.File fileName) {
        
        try {
            
            //Get random access to file:
            //java.io.RandomAccessFile file = new java.io.RandomAccessFile(fileName,"rw");
            java.io.InputStream in = new java.io.BufferedInputStream(new java.io.FileInputStream(fileName));

            //Get the size of the file
            long length = fileName.length();

            //Read contents, check that size is same as defined matrix size:
            long matrixLength = 2l * (long) matrix.length;
            if (length==matrixLength) {

                //Define the buffer size: 
                int bufferLength = 5000 * 1024;
                byte[] buffer = new byte[bufferLength];

                //Define length of read bytes:
                int len = 0;
                
                //Initialize location where data is and should go:
                int loc = 0;
                
                while((len = in.read(buffer, 0, bufferLength)) != -1){

                    // copy bytes to short:
                    for (int x=0; x<len; x+=2) {
                        
                        matrix[loc] = (short) (buffer[x]<<8 | (buffer[x+1]&0xff));
                        loc++;
                        
                    }

//                    System.out.print(".");
                    //if (progressDialog!=null) progressDialog.setProgress((int) ((double) loc * 100d / (double) matrix.length));
                }
                
//                System.out.println("");
                //if (progressDialog!=null) progressDialog.setProgress(100);
                
                //Dereference buffer:
                buffer = null;
                
            } else {
                System.out.println("File does not adhere to the matrix size! Actual file length:\t" + length + "\t, should be:\t" + matrix.length*2);
            }

            //Close file:
            in.close();

        } catch (java.io.IOException e) {
            System.out.println("Cannot read from file! (" + e.getMessage() + ")");
        }
        
    }
    
}
