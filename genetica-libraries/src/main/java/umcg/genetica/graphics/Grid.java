package umcg.genetica.graphics;

import com.itextpdf.text.DocumentException;
import umcg.genetica.graphics.panels.AssociationPanel;
import umcg.genetica.graphics.panels.Panel;

import java.io.IOException;

/**
 * Created by hwestra on 7/15/15.
 */
public class Grid {

    int marginBetweenPanels;
    int marginX;
    int marginY;
    int figureHeight;
    int figureWidth;
    int nrRows;
    int nrCols;
    int panelHeight;
    int panelWidth;

    Panel[][] grid;
    boolean[][] panelOccupied;
    boolean isInitialized = false;

    int currentPanel = 0;

    public enum SIZE {
        LETTER(2250, 3900),
        A4(2480, 3508);

        private int width;
        private int height;

        SIZE(int width, int height) {
            this.width = width;
            this.height = height;
        }
    }

    public Panel getPanel(int x, int y) {
        return grid[x][y];
    }

    public enum EXPAND {
        ROW,
        COL,
        BOTH
    }

    public Grid(SIZE size, int panelWidth, int panelHeight, int nrRows, int nrCols, int margin, int marginBetweenPanels) {
        this.marginX = margin;
        this.marginY = margin;
        this.marginBetweenPanels = marginBetweenPanels;
        this.nrCols = nrCols;
        this.nrRows = nrRows;
        this.panelWidth = panelWidth;
        this.panelHeight = panelHeight;
        initFixedSize(size, nrRows, nrCols);
    }

    public Grid(int panelWidth, int panelHeight, int nrRows, int nrCols, int margin, int marginBetweenPanels) {
        this.marginX = margin;
        this.marginY = margin;
        this.marginBetweenPanels = marginBetweenPanels;
        this.nrCols = nrCols;
        this.nrRows = nrRows;
        this.panelWidth = panelWidth;
        this.panelHeight = panelHeight;
        init(nrRows, nrCols);
        System.out.println("Grid: " + nrRows + "x" + nrCols);
    }

    private void initFixedSize(SIZE size, int nrRowsTmp, int nrColsTmp) {
        nrCols = nrColsTmp;
        nrRows = nrRowsTmp;
        System.out.println("Init: " + nrRowsTmp + " x " + nrColsTmp);

        this.figureHeight = size.height;
        this.figureWidth = size.width;

        if (isInitialized) {
            // copy contents
            Panel[][] tmpgrid = new Panel[nrRows][nrCols];
            boolean[][] tmpgridOccupied = new boolean[nrRows][nrCols];
            for (int i = 0; i < grid.length; i++) {
                for (int j = 0; j < grid[0].length; j++) {
                    tmpgrid[i][j] = grid[i][j];
                    tmpgridOccupied[i][j] = panelOccupied[i][j];
                }
            }
            grid = tmpgrid;
            panelOccupied = tmpgridOccupied;
        } else {
            grid = new Panel[nrRows][nrCols];
            panelOccupied = new boolean[nrRows][nrCols];
            isInitialized = true;
        }
    }

    private void init(int nrRowsTmp, int nrColsTmp) {
        nrCols = nrColsTmp;
        nrRows = nrRowsTmp;
        System.out.println("Init: " + nrRowsTmp + " x " + nrColsTmp);

        this.figureWidth = (nrCols * panelWidth) + ((nrCols - 1) * marginBetweenPanels) + (2 * marginX);
        if (this.figureWidth > 14400) {
            System.out.println("WARNING: figureWidth > 14400: " + this.figureWidth);
            this.figureWidth = 14000;
        }
        this.figureHeight = (nrRows * panelHeight) + ((nrRows - 1) * marginBetweenPanels) + (2 * marginY);
        if (this.figureHeight > 14400) {
            System.out.println("WARNING: figureheight > 14400: " + this.figureHeight);
            this.figureHeight = 14000;

        }

        if (isInitialized) {
            // copy contents
            Panel[][] tmpgrid = new Panel[nrRows][nrCols];
            boolean[][] tmpgridOccupied = new boolean[nrRows][nrCols];
            for (int i = 0; i < grid.length; i++) {
                for (int j = 0; j < grid[0].length; j++) {
                    tmpgrid[i][j] = grid[i][j];
                    tmpgridOccupied[i][j] = panelOccupied[i][j];
                }
            }
            grid = tmpgrid;
            panelOccupied = tmpgridOccupied;
        } else {
            grid = new Panel[nrRows][nrCols];
            panelOccupied = new boolean[nrRows][nrCols];
            isInitialized = true;
        }
    }

    // default is to fill from left to right..
    // perform magical auto-layout
    public void addPanel(Panel p) {
        int nrRowsPanel = p.getNrRows();
        int nrColsPanel = p.getNrCols();

        boolean panelPlaced = false;
        while (!panelPlaced) {
            for (int i = 0; i < grid.length; i++) {
                if (panelPlaced) {
                    break;
                } else {
                    for (int j = 0; j < grid[i].length; j++) {
                        if (!panelOccupied[i][j] && !panelPlaced) {
                            // check whether there is space near
                            if (nrRowsPanel == 1 && nrColsPanel == 1) {
                                grid[i][j] = p;
                                panelOccupied[i][j] = true;
                                panelPlaced = true;
                                break;
                            } else {
                                int nrCellsRequired = nrRowsPanel * nrColsPanel;
                                int availableCells = 0;
                                for (int row = i; row < i + nrRowsPanel && row < grid.length; row++) {
                                    for (int col = j; col < j + nrColsPanel && col < grid[i].length; col++) {
                                        if (!panelOccupied[row][col]) {
                                            availableCells++;
                                        }
                                    }
                                }

                                if (availableCells == nrCellsRequired) {
                                    grid[i][j] = p;
                                    for (int row = i; row < i + nrRowsPanel && row < grid.length; row++) {
                                        for (int col = j; col < j + nrColsPanel && col < grid[i].length; col++) {
                                            panelOccupied[row][col] = true;
                                        }
                                    }
                                    panelPlaced = true;
                                    break;
                                }
                            }
                        }
                    }
                }

            }
            if (!panelPlaced) {
                init(nrRows + 1, nrCols + 1);
            }
        }
    }

    // allow for user defined layouts
    public void addPanel(Panel p, int row, int col) {
        int nrRowsPanel = p.getNrRows();
        int nrColsPanel = p.getNrCols();
        grid[row][col] = p;
        for (int q = row; q < row + nrRowsPanel && q < grid.length; q++) {
            for (int r = col; r < col + nrColsPanel && r < grid[0].length; r++) {
                panelOccupied[q][r] = true;
            }
        }
    }

    public void draw(String filename) throws IOException, DocumentException {
        if (filename.endsWith(".png")) {
            draw(filename, DefaultGraphics.Output.PNG);
        } else if (filename.endsWith(".pdf")) {
            draw(filename, DefaultGraphics.Output.PDF);
        }
    }


    public void draw(String filename, DefaultGraphics.Output outputtype) throws IOException, DocumentException {

        // remove unneeded rows and columns
        trim();

        // set up the canvas
        DefaultGraphics g = new DefaultGraphics(filename, figureWidth, figureHeight);
        for (int row = 0; row < grid.length; row++) {
            for (int col = 0; col < grid[0].length; col++) {
                Panel p = grid[row][col];
                if (p != null) {
                    int starty = marginY + (row * panelHeight) + (row * marginBetweenPanels);
                    int startx = marginX + (col * panelWidth) + (col * marginBetweenPanels);


                    p.setX0(startx);
                    p.setY0(starty);

                    // System.out.println("Plotting!! " + startx + "\t" + starty);
                    int panelRows = p.getNrRows();
                    int panelCols = p.getNrCols();

                    // set pixel dimensions, for easy calculation when drawing
                    if (panelCols == 1 && panelRows == 1) {
                        p.setDimensions(panelWidth, panelHeight);
                    } else {
                        int panelW = (panelCols * panelWidth) + ((panelCols - 1) * marginBetweenPanels);
                        int panelH = (panelRows * panelHeight) + ((panelRows - 1) * marginBetweenPanels);
                        System.out.println(panelW + " x " + panelH + " panel size");
                        p.setDimensions(panelW, panelH);
                    }
                    if (p instanceof AssociationPanel) {
                        System.out.println("Pre: " + ((AssociationPanel) p).getTitle() + "\t" + ((AssociationPanel) p).getMaxP());
                        p.draw(g);
                        System.out.println("Post: " + ((AssociationPanel) p).getTitle() + "\t" + ((AssociationPanel) p).getMaxP());
//						System.exit(-1);
                    } else {
                        p.draw(g);
                    }

                }
            }
        }

        g.close();
    }

    private void trim() {

        boolean[] removeRow = new boolean[grid.length];
        int nrRowsToRemove = 0;
        boolean[] removeCol = new boolean[grid[0].length];
        int nrColsToRemove = 0;

        for (int row = 0; row < grid.length; row++) {
            int colsfilled = 0;
            for (int col = 0; col < grid[row].length; col++) {
                if (panelOccupied[row][col]) {
                    colsfilled++;
                }
            }
            if (colsfilled == 0) {
                removeRow[row] = true;
                nrRowsToRemove++;
            }
        }

        for (int col = 0; col < grid[0].length; col++) {
            int rowsFilled = 0;
            for (int row = 0; row < grid.length; row++) {
                if (panelOccupied[row][col]) {
                    rowsFilled++;
                }
            }
            if (rowsFilled == 0) {
                removeCol[col] = true;
                nrColsToRemove++;
            }
        }
        if (nrColsToRemove > 0 || nrRowsToRemove > 0) {

            nrRows = grid.length - nrRowsToRemove;
            nrCols = grid[0].length - nrColsToRemove;

            this.figureWidth = (nrCols * panelWidth) + ((nrCols - 1) * marginBetweenPanels) + (2 * marginX);
            this.figureHeight = (nrRows * panelHeight) + ((nrRows - 1) * marginBetweenPanels) + (2 * marginY);

            Panel[][] tmpGrid = new Panel[nrRows][nrCols];
            System.out.println(nrColsToRemove + " cols to remove " + nrRowsToRemove + " rows to remove");
            if (nrRows == 0 || nrCols == 0) {
                throw new IllegalStateException("Error: nothing to plot");
            }
            System.out.println("new grid: " + tmpGrid.length + " x " + tmpGrid[0].length);
            boolean[][] tmpgridOccupied = new boolean[nrRows][nrCols];
            for (int row = 0; row < nrRows; row++) {

                for (int col = 0; col < nrCols; col++) {
                    tmpGrid[row][col] = grid[row][col];
                    tmpgridOccupied[row][col] = panelOccupied[row][col];
                }
            }

            panelOccupied = tmpgridOccupied;
            grid = tmpGrid;


        }
    }
}
