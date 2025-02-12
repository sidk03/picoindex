/*
 * The pico index runs the code to parse the genome, create teh suffix array from it 
 * and then serialize the suffix array
 */


public class Picoindex {
    public static void main(String[] args) {
        String filePath = args[0];
        String outputName = args[1];
        FastaFile genomeFile = new FastaFile(filePath);
        SuffixArray suffix = new SuffixArray(genomeFile.records.get(0).seq, outputName);
        suffix.serializeSA();  
    }
}