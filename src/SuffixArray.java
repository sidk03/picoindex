import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

/*
 * This file creates a suffix array it uses a overridden compareTo method and sorts it using QuickSort
 */

public class SuffixArray {
    String genome;
    String outputName;
    int[] suffixArray;

    class Suffix implements Comparable<Suffix>{
        int offset;
        int length;

        private Suffix(int offset, int length){
            this.offset = offset;
            this.length = length;
        }

        @Override
        public int compareTo(Suffix other){
            int len = Math.min(this.length, other.length);
            for(int i = 0; i < len; i++){
                char c1 = genome.charAt(this.offset + i);
                char c2 = genome.charAt(other.offset + i);
                if (c1 != c2){
                    return c1 - c2;
                }
            }
            return this.length - other.length;
        }
    }

    public SuffixArray(String genome, String outputName){
        this.genome = genome + "$";
        this.outputName = outputName;
        this.suffixArray = this.suffixArrayBuilder();
    }

    private int[] suffixArrayBuilder(){
        int n = this.genome.length();
        Suffix[] suffixes = new Suffix[n];
        for(int i = 0; i < n; i++){
            suffixes[i] = new Suffix(i, n-i);
        }
        Arrays.sort(suffixes);
        int [] sufArr  = new int[n];
        for(int i = 0; i < n; i++){
            sufArr[i] = suffixes[i].offset;
        }
        return sufArr;
    }

    /*
     * Serializing the created suffix array into binary
     */

    public void serializeSA(){
        try {
            DataOutputStream dataStream = new DataOutputStream(new FileOutputStream(this.outputName));
            dataStream.writeInt(genome.length());
            dataStream.writeBytes(genome);
            dataStream.writeInt(suffixArray.length);
            for(int i = 0; i < suffixArray.length; i++){
                dataStream.writeInt(suffixArray[i]);
            }
            dataStream.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}


