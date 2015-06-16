/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.chrContacts;

/**
 *
 * @author MaKKie_Admin
 */
public class DesiredChrContact extends ChrContact{
    public final String snpName;
    public final String probeName;

  
    public DesiredChrContact(int chrLocSmal, int chrLocLarge, double contactVal, String snpName, String probeName) {
        super(chrLocSmal, chrLocLarge, contactVal);
        this.snpName = snpName;
        this.probeName = probeName;
    }

    public String getSnpName() {
        return snpName;
    }


    public String getProbeName() {
        return probeName;
    }    
    
}
