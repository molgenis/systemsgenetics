/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package eqtlmappingpipeline.graphics.curve;

/**
 *
 * @author harm-jan
 */
public class ROCCurve extends Curve {

    public ROCCurve(int width, int height, boolean outputPDF, String outputLoc){
        super(width, height, outputPDF, outputLoc);
    }

    // ROC curve: x-axis: 1-specificity (FPR), y-axis sensitivity (TPR)
    // x-axis and y-axis are always [0,1]
    public void plot(double[] TP, double[] TN){
        setMargins(50);
        calcScaling(1.0, 1.0);
        this.drawGrid(10);

        // draw x = y
        setColor(255,0,0,128);
        setStroke(2);
        
        int xCoord1 = getXCoord( 0.0 );
        int xCoord2 = getXCoord( 1.0 );
        int yCoord1 = getYCoord( 0.0 );
        int yCoord2 = getYCoord( 1.0 );
        drawLine(xCoord1, yCoord1, xCoord2, yCoord2);

        double xInterval = 1.0 / TP.length;

        setColor(0,0,0,128);
        setStroke(5);
        drawCurve(TP, TN, 0.0125);

        setColor(255,0,0,128);
        setStroke(5);
        drawCurve(TP, TN, 0.025);

        setColor(0,255,0,128);
        setStroke(5);
        drawCurve(TP, TN, 0.05);

        setColor(0,0,255,128);
        setStroke(5);
        drawCurve(TP, TN, 0.10);

        setColor(255,255,0,128);
        setStroke(5);
        drawCurve(TP, TN, 0.25);

        setColor(0,255,255,128);
        setStroke(5);
        drawCurve(TP, TN, 0.50);
    }

    public void drawCurve(double[] TP, double[] TN, double apriori){
        
        for(int i=0; i<TP.length-1; i++){
            int xCoord1 = getXCoord(TN[i]);
            int xCoord2 = getXCoord(TN[i+1]);
            int yCoord1 = getYCoord(TP[i]);
            int yCoord2 = getYCoord(TP[i+1]);

            drawLine(xCoord1, yCoord1, xCoord2, yCoord2);

        }
    }

    public void drawGrid(int numLines){
        java.text.DecimalFormat df = new java.text.DecimalFormat("0.00");

        // start by drawing a grid...
        // vertical lines
        for(int i=0; i<10; i++){

            if(i==0){
                setColor(0,0,0,255);
                setStroke(2);
            } else {
                setColor(0,0,0,128);
                setStroke(1);
            }

            double linelbl   = (double) i / numLines;

             // vertical lines
            int xCoord1 = getXCoord( linelbl );
            int xCoord2 = xCoord1;
            int yCoord1 = getYCoord( 1.0 );
            int yCoord2 = getYCoord( 0.0 );

            drawLine(xCoord1, yCoord1, xCoord2, yCoord2);

            String lineLabel = df.format(linelbl);
            int labelHeight  = getStringHeight(lineLabel);
            int labelWidth   = getStringWidth(lineLabel);
            
            if(i>0){    
                drawText(xCoord1, yCoord2 + labelWidth + 5, -90, lineLabel);
            }

            // horizontal lines
            xCoord1 = getXCoord( 0.0 );
            xCoord2 = getXCoord( 1.0 );
            yCoord1 = getYCoord( linelbl );
            yCoord2 = yCoord1;

            drawLine(xCoord1, yCoord1, xCoord2, yCoord2);

            if(i>0){
                drawText(xCoord1 - labelWidth - 5, yCoord2, lineLabel);
            }
        }

        setColor(0,0,0,128);
        setStroke(1);

        String lineLabel = df.format(0.0);
        int labelHeight  = getStringHeight(lineLabel);
        int labelWidth   = getStringWidth(lineLabel);

        int xCoord1 = getXCoord( 0.0 );
        int xCoord2 = xCoord1;
        int yCoord1 = getYCoord( 0.0 );
        int yCoord2 = yCoord1;
        drawText(xCoord1 - labelWidth, yCoord1 + labelHeight + 10, -45, lineLabel);

        xCoord1 = getXCoord( 0.0 );
        yCoord1 = getYCoord( 1.0 );

        lineLabel       = df.format(1.0);
        labelHeight     = getStringHeight(lineLabel);
        labelWidth      = getStringWidth(lineLabel);
        drawText(xCoord1 - labelWidth - 5, yCoord1 + 5, lineLabel);

        xCoord1 = getXCoord( 1.0 );
        yCoord1 = getYCoord( 0.0 );
        drawText(xCoord1, yCoord1 + labelWidth + 5, -90, lineLabel);

        setColor(0,0,0,255);
        setStroke(2);
        // draw top line and right line
        // vertical lines
        xCoord1 = getXCoord( 1.0 );
        xCoord2 = xCoord1;
        yCoord1 = getYCoord( 1.0 );
        yCoord2 = getYCoord( 0.0 );

        drawLine(xCoord1, yCoord1, xCoord2, yCoord2);

        // horizontal lines
        xCoord1 = getXCoord( 0.0 );
        xCoord2 = getXCoord( 1.0 );
        yCoord1 = getYCoord( 1.0 );
        yCoord2 = yCoord1;

        drawLine(xCoord1, yCoord1, xCoord2, yCoord2);

        setStroke(1);

        
    }

    public void plot(double[][] TPR, double[][] FPR, double[] AUC, double apriori) {
	
	setMargins(50);
        calcScaling(1.0, 1.0);
        this.drawGrid(10);

	setTitle("Reciever Operator Curve");
	
	setColor(255,255,255,255);
	drawRect(graphWidth - 40 - 200, graphHeight- 40  - 200, 150, 150, true);

        // draw x = y
        setColor(255,0,0,128);
        setStroke(2);

        int xCoord1 = getXCoord( 0.0 );
        int xCoord2 = getXCoord( 1.0 );
        int yCoord1 = getYCoord( 0.0 );
        int yCoord2 = getYCoord( 1.0 );
        drawLine(xCoord1, yCoord1, xCoord2, yCoord2);

	java.text.DecimalFormat df = new java.text.DecimalFormat("0.000");
	java.text.DecimalFormat df2 = new java.text.DecimalFormat("0.0");

	drawText(graphWidth - 20 - 200 + 20, graphHeight - 20 - 200 + 10, "Apriori");
	drawText(graphWidth - 20 - 200 + 20, graphHeight - 20 - 200 + 20, "P-value");

	drawText(graphWidth - 20 - 200 + 70, graphHeight - 20 - 200 + 20, "AUC");

        for(int i=TPR.length-1; i!=1; i--){
            

	    if(i == 5){
		setColor(0,0,0,255);
	    }
	    
	    if(i == 4){
		setColor(0,0,0,225);
	    }
	    if(i == 3){
		setColor(0,0,0,200);
	    }

	    if(i == 2){
		setColor(0,0,0,175);
	    }

	    if(i == 1){
		setColor(0,0,0,150);
	    }
            setStroke(5);
            drawCurve(TPR[i], FPR[i], apriori);

	    drawText(graphWidth - 20 - 200 + 20, graphHeight - 20 - 200 + 30 + ((5-i)*10), df2.format( 1-(i*apriori) ));
	    drawText(graphWidth - 20 - 200 + 70, graphHeight - 20 - 200 + 30 + ((5-i)*10), df.format(AUC[i]));
        }


	


        
    }
}
