package betaqtl.junk;

public class SplitTest {

    public static void main(String[] args) {

        String data = "-";
        for (int i = 0; i < 100; i++) {
            if (i == 50) {
                data += "\t0.1";
            } else {
                data += "\t";
            }
        }
        String[] elems = data.split("\t");
        System.out.println(elems.length);
        for(int i=0;i<elems.length;i++){
            System.out.println(i+"\t"+elems[i]);
        }
    }
}
