package nl.systemsgenetics.datg;

import static org.testng.Assert.*;

public class DatgConverterTest {

    @org.testng.annotations.BeforeMethod
    public void setUp() {



    }

    @org.testng.annotations.AfterMethod
    public void tearDown() {
    }

    @org.testng.annotations.Test
    public void testMain() {

        try {
            DatgConverter.main(new String[]{""});
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }

    }
}