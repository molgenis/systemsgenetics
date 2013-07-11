/*
 * SymmetricByteDistancaMatrix.java
 *
 * Created on 04 June 2004, 15:53
 */

package umcg.genetica.math.matrix;

/**
 *
 * Symmetric Byte Distance Matrix: A memory efficient matrix capable of performing calculations on all genes x all genes:
 *
 * @author  Lude Franke
 */
public class SymmetricByteDistanceMatrix {

    private int size;
    private byte[] matrix;
    private final static double eLog10 = java.lang.Math.log(10);
    private final static int MAX_ALL_PAIRS_STEPS = 100;
    public final static int MAX_VALUE = Byte.MAX_VALUE + 128;
    private long[] elementIndex;
    public boolean loadSuccess = false;
    
    /** Creates a new instance of SymmetricByteDistancaMatrix
     *
     * Defines a new Symmetric matrix, with a predefined size
     *
     * Initially all values will be Byte.MAX_Value (127)
     *
     * All data is stored in memory efficient one dimensional array, which costs size * (size + 1) / 2 bytes
     *
     */

    public SymmetricByteDistanceMatrix(int size) {
    
        this.size = size;
        long arraySize = ((long) size * (long) (size + 1)) / 2l;
        matrix = new byte[(int) arraySize];
        elementIndex = new long[size]; for (int x=0; x<size; x++) elementIndex[x] = (long) x * (long) size - ((long) x * (long)(x+1)) / 2l;
        setMaxDistance();
        
    }
    
    public SymmetricByteDistanceMatrix(int size, boolean setMaxDistance) {
        
        this.size = size;
        long arraySize = ((long) size * (long) (size + 1)) / 2l;
        matrix = new byte[(int) arraySize];
        elementIndex = new long[size]; for (int x=0; x<size; x++) elementIndex[x] = (long) x * (long) size - ((long) x * (long)(x+1)) / 2l;
        if (setMaxDistance) setMaxDistance();
    }

    public void dispose() {
        size = 0;
        matrix = null;
        elementIndex = null;
    }

    public void setAllElements(int value) {
        for (int m=0; m<matrix.length; m++) matrix[m] = (byte) (value - 128);;
    }
    
    private void setMaxDistance() {
        for (int m=0; m<matrix.length; m++) matrix[m] = Byte.MAX_VALUE;
    }
    
    /* Get the element locaiton for point (x,y):
     */
    private long getElement(int x, int y) {
        if (x>y) {
            return elementIndex[y] + (long) x;
        } else {
            return elementIndex[x] + (long) y;
        }
    }
    
    /* Set the value for point (x,y), enter value between 0 and 255:
     */
    public void set(int x, int y, int value) {
        matrix[(int) getElement(x,y)] = (byte) (value - 128);
    }
    
    /* Get the value for point (x,y), returns value between 0 and 255:
     */
    public int get(int x, int y) {
        return matrix[(int) getElement(x,y)] + 128;
    }
    
    /* Get size of matrix:
     */
    public int size() {
        return size;
    }
     
    public int maxValue() {
        return Byte.MAX_VALUE + 128;
    }
    
    /* Get the shortest path for all the pairs using the Floyd-Warshall algorithm: 
     *
     * Please be aware that this method replaces the values inside the current matrix with the shortest path values!
     */
    public void getAllPairsShortestPath() {

        //This method uses the Floyd-Warshall all pairs shortest path algorithm:
        //As we deal here with a symmetric matrix, the number of required calculations is 50% less.

        //Traverse data:
        int previousPercentage = 0;
        for (int i=0; i<size()&&i<MAX_ALL_PAIRS_STEPS; i++) {
            for (int v=0; v<size(); v++) {
                for (int w=v; w<size(); w++) {
                    if (get(v,i)!=maxValue()&&get(i,w)!=maxValue()) {
                        set(v, w, (short) java.lang.Math.min((int) get(v, w), (int) get(v, i) + (int) get(i, w)));
                    }
                }
            }
            int percentage = (int) (double) (100d * (double) i / (double) size());
            if (percentage!=previousPercentage) System.out.print(percentage + "% ");
            previousPercentage = percentage;
            System.out.print(".");
            
        }
        System.out.println("");
        
        //Set the total distance to genes itself to 0:
        for (int v=0; v<size(); v++) {
            set(v,v,0);
        }
        
        
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
        boolean traversed[] = new boolean[size()];
        
        //Initialise T: Two dimensional array, column 0: distance, column 1: nodeID predecessor:
        int T[][] = new int[2][size()];   
        for (int i=0; i<size(); i++) {
            traversed[i] = false;
            T[0][i] = maxValue();
            T[1][i] = -1;
        }
        
        //Starting from nodeID: nodeStart: set traversed true:
        traversed[x] = true;
        
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
        
        System.out.println("pathPossible: " + pathPossible);
        
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
                if (!traversed[i]&&T[0][i]<tMin) {
                    
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
            
            System.out.println("pathExists: " + pathExists);
                
            //The nearest not finalized node v has been found, process this node:
            traversed[v] = true;
        
            //Traverse all nodes:
            for (int i=0; i<size(); i++) {
            
                //Find node whose final shortest path has not been established and which is directly connected:
                if (!traversed[i]&&get(v,i)>0) {
                    
                    //Check whether the distance of the edge between the processed node and node v is the lowest:
                    if (T[0][v] + get(v, i) < T[0][i]) {
                        
                        //Change the distance and nodeIndex of the i:
                        T[0][i] = T[0][v] + (int) get(v,i);
                        T[1][i] = v;
                        
                    }
                }
            }
            
            //Check whether destination has been reached:
            if (traversed[y]) { 
                traversePath = false;
            } else {
                //Check whether the path still between nodeStart and nodeEnd is still likely to exist:
                pathExists = false;
                for (int i=0; i<size(); i++) {
                    if (!traversed[i]&&T[0][i] < maxValue()) {
                        pathExists = true;
                        break;
                    }
                }
            }
        }
        
        System.out.println("pathExists: " + pathExists + "\tpathPossible: " + pathPossible + "\tx: " + x + "\ty: " + y);
        
        //Check whether the path exists and is possible:
        if (pathExists&&pathPossible) {

            //Select nodes that have been involved, start with nodeEnd:
            //System.out.println(y);
            path.add(new Integer(y));
            
            //Select predecessor nodes:
            int totalDistance = T[0][y];
            int predecessor = T[1][y];
            
            //System.out.println(predecessor);
            path.add(new Integer(predecessor));
            
            while (predecessor!=x) {
                totalDistance += T[0][predecessor];
                predecessor = T[1][predecessor];
                System.out.println(predecessor);
                path.add(new Integer(predecessor));
            }

        }

        
        
        return path;
        
    }
    
    public int getDistance(double pValue) {
        return (int) (java.lang.Math.max(java.lang.Math.log(pValue) / eLog10, -3.984375) * 64 + 255);
    }

    public double getPValue(int distance) {
        return (double) java.lang.Math.pow(10, ((double) distance - 255) / 64);
    }

    public void save(java.io.File fileName) {
        
        try {
            
            //Initialize output stream:
            java.io.OutputStream out = new java.io.BufferedOutputStream(new java.io.FileOutputStream(fileName));

            //Define the buffer size: 
            int bufferLength = 5000 * 1024;
            byte[] buffer = new byte[bufferLength];

            //Define length of read bytes:
            int len = 0;

            //Initialize location where data is and should go:
            int loc = 0;

            //(len = in.read(buffer, 0, bufferLength)) != -1
            while(1==1){

                // Copy data
                if (loc + bufferLength <= matrix.length) {
                    len = bufferLength;
                    System.arraycopy(matrix, loc, buffer, 0, len);
                    loc+= len;
                    out.write(buffer);
                    System.out.print(".");
                } else {
                    len = matrix.length - loc;
                    buffer = new byte[len];
                    System.arraycopy(matrix, loc, buffer, 0, len);
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
            if (length==matrix.length) {

                //Define the buffer size: 
                int bufferLength = 5000 * 1024;
                byte[] buffer = new byte[bufferLength];

                //Define length of read bytes:
                int len = 0;
                
                //Initialize location where data is and should go:
                int loc = 0;
                
                while((len = in.read(buffer, 0, bufferLength)) != -1){

                    // copy data
                    System.arraycopy(buffer, 0, matrix, loc, len);
                    loc+= len;

                    System.out.print(".");
                }
                
                System.out.println("");
                
                //Dereference buffer:
                buffer = null;
                loadSuccess = true;
                
            } else {
                System.out.println("File does not adhere to the matrix size!");
                loadSuccess = false;
            }

            //Close file:
            in.close();

        } catch (Exception e) {
            System.out.println("Cannot read from file! (" + e.getMessage() + ")");
        }
        
    }
    
}
