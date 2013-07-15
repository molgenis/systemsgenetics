/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.util;

import java.util.ArrayList;

/**
 * As taken from jsc.jar
 *
 * @author MarcJan
 */
public class Sort {

    public static ArrayList<String> getLabels(String[] paramArrayOfString) {
        int i = paramArrayOfString.length;

        String[] arrayOfString = new String[i];
        System.arraycopy(paramArrayOfString, 0, arrayOfString, 0, i);

        sortStrings(arrayOfString, null, 0, i - 1, true);
        ArrayList<String> localVector = new ArrayList<String>();
        for (int j = 0; j < i; j++) {
            if (!localVector.contains(arrayOfString[j])) {
                localVector.add(arrayOfString[j]);
            }
        }
        return localVector;
    }

    public static void sort(double[] paramArrayOfDouble, int[] paramArrayOfInt, int paramInt1, int paramInt2, boolean paramBoolean) {
        if ((paramArrayOfDouble == null) || (paramArrayOfDouble.length < 2)) {
            return;
        }

        int i = paramInt1;
        int j = paramInt2;
        Double localDouble = new Double(paramArrayOfDouble[((paramInt1 + paramInt2) / 2)]);
        do {
            if (paramBoolean) {
                while (true) {
                    i++;
                    if (i < paramInt2) {
                        if (localDouble.compareTo(new Double(paramArrayOfDouble[i])) <= 0) {
                            break;
                        }
                    }
                }
                do {
                    j--;
                    if (j <= paramInt1) {
                        break;
                    }
                } while (localDouble.compareTo(new Double(paramArrayOfDouble[j])) < 0);
            } else {
                do {
                    i++;
                    if (i >= paramInt2) {
                        break;
                    }
                } while (localDouble.compareTo(new Double(paramArrayOfDouble[i])) < 0);
                while ((j > paramInt1) && (localDouble.compareTo(new Double(paramArrayOfDouble[j])) > 0)) {
                    j--;
                }
            }
            if (i < j) {
                double d = paramArrayOfDouble[i];
                paramArrayOfDouble[i] = paramArrayOfDouble[j];
                paramArrayOfDouble[j] = d;
                if (paramArrayOfInt != null) {
                    int k = paramArrayOfInt[i];
                    paramArrayOfInt[i] = paramArrayOfInt[j];
                    paramArrayOfInt[j] = k;
                }
            }
            if (i <= j) {
                i++;
                j--;
            }
        } while (i <= j);
        if (paramInt1 < j) {
            sort(paramArrayOfDouble, paramArrayOfInt, paramInt1, j, paramBoolean);
        }
        if (i < paramInt2) {
            sort(paramArrayOfDouble, paramArrayOfInt, i, paramInt2, paramBoolean);
        }
    }

    public static void sort(Double[] paramArrayOfDouble, int[] paramArrayOfInt, int paramInt1, int paramInt2, boolean paramBoolean) {
        if ((paramArrayOfDouble == null) || (paramArrayOfDouble.length < 2)) {
            return;
        }

        int i = paramInt1;
        int j = paramInt2;
        Double localDouble1 = paramArrayOfDouble[((paramInt1 + paramInt2) / 2)];
        do {
            if (paramBoolean) {
                while (true) {
                    i++;
                    if (i < paramInt2) {
                        if (localDouble1.compareTo(paramArrayOfDouble[i]) <= 0) {
                            break;
                        }
                    }
                }
                do {
                    j--;
                    if (j <= paramInt1) {
                        break;
                    }
                } while (localDouble1.compareTo(paramArrayOfDouble[j]) < 0);
            } else {
                do {
                    i++;
                    if (i >= paramInt2) {
                        break;
                    }
                } while (localDouble1.compareTo(paramArrayOfDouble[i]) < 0);
                while ((j > paramInt1) && (localDouble1.compareTo(paramArrayOfDouble[j]) > 0)) {
                    j--;
                }
            }
            if (i < j) {
                Double localDouble2 = paramArrayOfDouble[i];
                paramArrayOfDouble[i] = paramArrayOfDouble[j];
                paramArrayOfDouble[j] = localDouble2;
                if (paramArrayOfInt != null) {
                    int k = paramArrayOfInt[i];
                    paramArrayOfInt[i] = paramArrayOfInt[j];
                    paramArrayOfInt[j] = k;
                }
            }
            if (i <= j) {
                i++;
                j--;
            }
        } while (i <= j);
        if (paramInt1 < j) {
            sort(paramArrayOfDouble, paramArrayOfInt, paramInt1, j, paramBoolean);
        }
        if (i < paramInt2) {
            sort(paramArrayOfDouble, paramArrayOfInt, i, paramInt2, paramBoolean);
        }
    }

    public static void sortDoubles(double[] paramArrayOfDouble1, double[] paramArrayOfDouble2, int paramInt1, int paramInt2, boolean paramBoolean) {
        if ((paramArrayOfDouble1 == null) || (paramArrayOfDouble1.length < 2)) {
            return;
        }

        int i = paramInt1;
        int j = paramInt2;
        Double localDouble = new Double(paramArrayOfDouble1[((paramInt1 + paramInt2) / 2)]);
        do {
            if (paramBoolean) {
                while (true) {
                    i++;
                    if (i < paramInt2) {
                        if (localDouble.compareTo(new Double(paramArrayOfDouble1[i])) <= 0) {
                            break;
                        }
                    }
                }
                do {
                    j--;
                    if (j <= paramInt1) {
                        break;
                    }
                } while (localDouble.compareTo(new Double(paramArrayOfDouble1[j])) < 0);
            } else {
                do {
                    i++;
                    if (i >= paramInt2) {
                        break;
                    }
                } while (localDouble.compareTo(new Double(paramArrayOfDouble1[i])) < 0);
                while ((j > paramInt1) && (localDouble.compareTo(new Double(paramArrayOfDouble1[j])) > 0)) {
                    j--;
                }
            }
            if (i < j) {
                double d = paramArrayOfDouble1[i];
                paramArrayOfDouble1[i] = paramArrayOfDouble1[j];
                paramArrayOfDouble1[j] = d;
                if (paramArrayOfDouble2 != null) {
                    d = paramArrayOfDouble2[i];
                    paramArrayOfDouble2[i] = paramArrayOfDouble2[j];
                    paramArrayOfDouble2[j] = d;
                }
            }
            if (i <= j) {
                i++;
                j--;
            }
        } while (i <= j);
        if (paramInt1 < j) {
            sortDoubles(paramArrayOfDouble1, paramArrayOfDouble2, paramInt1, j, paramBoolean);
        }
        if (i < paramInt2) {
            sortDoubles(paramArrayOfDouble1, paramArrayOfDouble2, i, paramInt2, paramBoolean);
        }
    }

    public static void sortInts(int[] paramArrayOfInt1, int[] paramArrayOfInt2, int paramInt1, int paramInt2, boolean paramBoolean) {
        if ((paramArrayOfInt1 == null) || (paramArrayOfInt1.length < 2)) {
            return;
        }

        int i = paramInt1;
        int j = paramInt2;
        Integer localInteger = new Integer(paramArrayOfInt1[((paramInt1 + paramInt2) / 2)]);
        do {
            if (paramBoolean) {
                while (true) {
                    i++;
                    if (i < paramInt2) {
                        if (localInteger.compareTo(new Integer(paramArrayOfInt1[i])) <= 0) {
                            break;
                        }
                    }
                }
                do {
                    j--;
                    if (j <= paramInt1) {
                        break;
                    }
                } while (localInteger.compareTo(new Integer(paramArrayOfInt1[j])) < 0);
            } else {
                do {
                    i++;
                    if (i >= paramInt2) {
                        break;
                    }
                } while (localInteger.compareTo(new Integer(paramArrayOfInt1[i])) < 0);
                while ((j > paramInt1) && (localInteger.compareTo(new Integer(paramArrayOfInt1[j])) > 0)) {
                    j--;
                }
            }
            if (i < j) {
                int k = paramArrayOfInt1[i];
                paramArrayOfInt1[i] = paramArrayOfInt1[j];
                paramArrayOfInt1[j] = k;
                if (paramArrayOfInt2 != null) {
                    int m = paramArrayOfInt2[i];
                    paramArrayOfInt2[i] = paramArrayOfInt2[j];
                    paramArrayOfInt2[j] = m;
                }
            }
            if (i <= j) {
                i++;
                j--;
            }
        } while (i <= j);
        if (paramInt1 < j) {
            sortInts(paramArrayOfInt1, paramArrayOfInt2, paramInt1, j, paramBoolean);
        }
        if (i < paramInt2) {
            sortInts(paramArrayOfInt1, paramArrayOfInt2, i, paramInt2, paramBoolean);
        }
    }

    public static void sortStrings(String[] paramArrayOfString1, String[] paramArrayOfString2, int paramInt1, int paramInt2, boolean paramBoolean) {
        if ((paramArrayOfString1 == null) || (paramArrayOfString1.length < 2)) {
            return;
        }

        int i = paramInt1;
        int j = paramInt2;
        String str1 = paramArrayOfString1[((paramInt1 + paramInt2) / 2)];
        do {
            if (paramBoolean) {
                while (true) {
                    i++;
                    if (i < paramInt2) {
                        if (compareStrings(str1, paramArrayOfString1[i]) <= 0) {
                            break;
                        }
                    }
                }
                do {
                    j--;
                    if (j <= paramInt1) {
                        break;
                    }
                } while (compareStrings(str1, paramArrayOfString1[j]) < 0);
            } else {
                do {
                    i++;
                    if (i >= paramInt2) {
                        break;
                    }
                } while (compareStrings(str1, paramArrayOfString1[i]) < 0);
                while ((j > paramInt1) && (compareStrings(str1, paramArrayOfString1[j]) > 0)) {
                    j--;
                }
            }
            if (i < j) {
                String str2 = paramArrayOfString1[i];
                paramArrayOfString1[i] = paramArrayOfString1[j];
                paramArrayOfString1[j] = str2;
                if (paramArrayOfString2 != null) {
                    str2 = paramArrayOfString2[i];
                    paramArrayOfString2[i] = paramArrayOfString2[j];
                    paramArrayOfString2[j] = str2;
                }
            }
            if (i <= j) {
                i++;
                j--;
            }
        } while (i <= j);
        if (paramInt1 < j) {
            sortStrings(paramArrayOfString1, paramArrayOfString2, paramInt1, j, paramBoolean);
        }
        if (i < paramInt2) {
            sortStrings(paramArrayOfString1, paramArrayOfString2, i, paramInt2, paramBoolean);
        }
    }

    public static int compare(Object paramObject1, Object paramObject2, boolean paramBoolean) {
        if ((paramObject1 != null) && (paramObject2 == null)) {
            return paramBoolean ? -1 : 1;
        }
        if ((paramObject1 == null) && (paramObject2 == null)) {
            return 0;
        }
        if ((paramObject1 == null) && (paramObject2 != null)) {
            return paramBoolean ? 1 : -1;
        }
        Object localObject1;
        Object localObject2;
        if (((paramObject1 instanceof Double)) && ((paramObject2 instanceof Double))) {
            localObject1 = (Double) paramObject1;
            localObject2 = (Double) paramObject2;
            return ((Double) localObject1).compareTo((Double) localObject2);
        }
        if (((paramObject1 instanceof Integer)) && ((paramObject2 instanceof Integer))) {
            localObject1 = (Integer) paramObject1;
            localObject2 = (Integer) paramObject2;
            return ((Integer) localObject1).compareTo((Integer) localObject2);
        }
        return compareStrings(paramObject1.toString(), paramObject2.toString());
    }

    public static int compareStrings(String paramString1, String paramString2) {
        try {
            Double localDouble1 = new Double(Double.parseDouble(paramString1));
            Double localDouble2 = new Double(Double.parseDouble(paramString2));
            return localDouble1.compareTo(localDouble2);
        } catch (NumberFormatException localNumberFormatException) {
        }
        return paramString1.compareToIgnoreCase(paramString2);
    }

    static class Test {

        public static void main(String[] paramArrayOfString) {
            double[] arrayOfDouble = {7.0D, 3.0D, -4.0D, 3.0D, 1.0D, 8.0D, 2.0D, 1.0D, 6.0D, -5.0D, -2.0D, 0.0D, 7.0D};
            Sort.sortDoubles(arrayOfDouble, null, 0, arrayOfDouble.length - 1, true);
            for (int i = 0; i < arrayOfDouble.length; i++) {
                System.out.println(arrayOfDouble[i] + ", ");
            }
            String[] arrayOfString = {"X-ray", "1211", "0", "Golf", "99", "x-ray", "111", "Hotel", "Alpha", "Zero", "25", "Juliette", "£", "Foxtrot", "Tango", "Romeo", "1111", "Bravo", "November", "Charlie", "*", "-10", "Victor", "10", "romeo", "12", "11", "121", "1", "-100", "-1", "-0"};

            Sort.sortStrings(arrayOfString, null, 0, arrayOfString.length - 1, true);
            for (int i = 0; i < arrayOfString.length; i++) {
                System.out.println(arrayOfString[i] + ", ");
            }
        }
    }
}
