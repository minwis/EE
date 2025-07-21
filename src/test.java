package uk.ac.ebi.uniprot.dataservice.client.examples;//package uk.ac.ebi.uniprot.dataservice.client.examples;

import jdk.incubator.vector.*;

//in Terminal, type: java --add-modules jdk.incubator.vector test.java

import java.util.Arrays;

public class test {

    static final VectorSpecies<Integer> SPECIES = IntVector.SPECIES_256;
    static final VectorSpecies<Float> SPECIES2 = FloatVector.SPECIES_64;

    public static void main(String[] args) {

        int[] a = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

        int[] b = new int[] {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
        int[] c = new int[a.length];

        for( int i = 0; i < a.length; i+=SPECIES.length() ) {
            var mask = SPECIES.indexInRange(i, a.length);
            var va = IntVector.fromArray(SPECIES, a, i, mask);
            var vb = IntVector.fromArray(SPECIES, b, i, mask);
            var vc = va.add(vb, mask);
            vc.intoArray(c, i, mask);

        }

        IntVector vc = IntVector.fromArray(SPECIES, c, 0);

        System.out.println(vc);
        System.out.println(vc.lanewise(VectorOperators.LSHL, 2));


        String s = "hello world";
        int n = 'h';
        System.out.println(s.charAt(0) == n);
    }

}
