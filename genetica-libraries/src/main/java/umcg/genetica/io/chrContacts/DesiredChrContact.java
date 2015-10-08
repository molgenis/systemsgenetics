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
public class DesiredChrContact extends ChrContact {

    private final String snpName;
    private final String probeName;
    private boolean contact = false;
    private double contactValue = 0.0d;
    private double normalizedContactValue = 0.0d;

    public boolean isContact() {
        return contact;
    }

    public void setContact(boolean contact) {
        this.contact = contact;
    }

    public double getNormalizedContactValue() {
        return normalizedContactValue;
    }

    public void setNormalizedContactValue(double normalizedContactValue) {
        this.normalizedContactValue = normalizedContactValue;
    }
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

    public boolean hasContact() {
        return contact;
    }

    public void setContact() {
        this.contact = true;
    }

    public DesiredChrContact clone() {
        DesiredChrContact d = new DesiredChrContact(this.chrLocationSmaller, this.chrLocationLarger, this.contactValue, this.snpName, this.probeName);
        return d;
    }

}
