import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;

public class Picomap {

    /*
     * 
     */
    private class Locations{
        boolean exactMatch;
        int numLocations;
        ArrayList<Integer> matchLocations;
        ArrayList<String> cigarStrings;
        int score;


        private Locations(boolean exactMatch, int numLocations, ArrayList<Integer> matchLocations){
            this.exactMatch = exactMatch;
            this.numLocations = numLocations;
            this.matchLocations = matchLocations;
        }
    }

    FastaFile readsFile;
    String genome;
    int[] suffixArray;
    int mismatch_penalty;
    int gap_penalty;
    String outputFile;
    
    public Picomap(String index_file, String read_file, int mismatch_penalty, int gap_penalty, String output_file){
        this.deserializeSA(index_file);
        this.readsFile = new FastaFile(read_file);
        this.mismatch_penalty = mismatch_penalty;
        this.gap_penalty = gap_penalty;
        this.outputFile = output_file;
    }

    /*
     * De serializes the suffix array binary file
     */
    private void deserializeSA(String filePath){
        try {
            DataInputStream dataStream = new DataInputStream(new FileInputStream(filePath));
            // Read in genome
            int geneLength = dataStream.readInt();
            byte[] bytes = new byte[geneLength];
            dataStream.readFully(bytes);
            this.genome = new String(bytes, StandardCharsets.UTF_8);
            // Read in SuffixArray
            int saLength = dataStream.readInt();
            int[] sufArr = new int[saLength];
            for(int i = 0; i < saLength; i++){
                sufArr[i] = dataStream.readInt();
            }
            this.suffixArray = sufArr;
            dataStream.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /*
     * Does an accelrated search of the suffix array using the longest common prefix
     */
    private int[] accelBs(String query){
        int l = 0;
        int r = this.suffixArray.length;
        int comp = 0;
        int lcpL = 0;
        int lcpR = 0;

        while(l < r){
            int c = (l+r)/2;
            int[] res = this.compareSuffixAccel(query,this.suffixArray[c],lcpL,lcpR);
            comp += res[1];
            if(res[0] < 0){
                if(c == l+1){
                    return new int[]{c,comp};
                }
                r = c;
                lcpR = res[2];
            }else{
                if(c == r-1){
                    return new int[]{r,comp};
                }
                l = c;
                lcpL = res[2];
            }
        }
        return new int[]{0,0};
    }

    /*
     * Compares the given query to an entry in the suffix array
     * skips comparisons by using the longest common prefix
     */
    private int[] compareSuffixAccel(String q, int s, int lcpL, int lcpR){
        int comp = 0;
        int n = genome.length();
        int minLen = Math.min(q.length(), n-s);
        int i;
        for(i = Math.min(lcpL,lcpR); i < minLen; i++){
            char cQuery = q.charAt(i);
            char cSuffix = this.genome.charAt(s+i);
            comp += 1;
            if (cQuery != cSuffix){
                return new int[]{cQuery - cSuffix,comp,i};
            }
        }
        return new int[]{q.length() - (n-s),comp,i};
    }

    /*
     * Divides a large query into smaller substrings of length 25
     */
    private ArrayList<String> divideQuery(String query){
        ArrayList<String> substrings = new ArrayList<>();
        for (int i = 0; i < query.length(); i += 25) {
            substrings.add(query.substring(i, Math.min(i + 25, query.length())));
        }
        return substrings;
    }

    private int[] exactMatch(String query){
        int l = this.accelBs(query+"#")[0];
        int r = this.accelBs(query+"}")[0];
        if(l == r){ 
            return new int[] {0};
        } else{
            return new int[] {r-l, l, r};
        }
    }

    /*
     * Heuristic to find potential places in out large genome where the query may match
     */
    private Locations potentialCandidates(String query){

        // Runs an exact search with the entire query to see if there are any exact matches
        int[] exactCheck = this.exactMatch(query);
        
        // If there are exact matches we add that to our locations class
        // We change exactMatch boolean to true
        if(exactCheck[0] != 0){
            int l = exactCheck[1];
            int r = exactCheck[2];
            ArrayList<Integer> exactLocations = new ArrayList<>();
            for(int i = l; i < r; i++){
                exactLocations.add(this.suffixArray[i]);
            }
            return new Locations(true, r-l, exactLocations);
        }

        // No exact matches were found time to use our heuristic to find potential matches
        else{
            // We divide our query
            ArrayList<String> seeds = this.divideQuery(query);
            ArrayList<ArrayList<Integer>> seedMap = new ArrayList<>();

            // For each seed (mini query) we search out genome for potential matches 
            for(String s : seeds){
                int l = this.accelBs(s+"#")[0];
                int r = this.accelBs(s+"}")[0];

                // Add potential locations to our seedMap
                ArrayList<Integer> seedLocations = new ArrayList<>();
                for(int i = l; i < r; i++){
                    seedLocations.add(this.suffixArray[i]);
                }
                seedMap.add(seedLocations);
            }
            ArrayList<Integer> seedLocations = this.filterCandidates(seedMap);
            ArrayList<ProblemOutput> output = this.alignCandidates(query, seedLocations);
            ArrayList<Integer> matchLocations = new ArrayList<>();
            ArrayList<String> cigarStrings = new ArrayList<>();
            for(ProblemOutput out : output){
                matchLocations.add(out.y_start);
                cigarStrings.add(out.cigar);
            }
            Locations loc = new Locations(false, output.size(), matchLocations);
            loc.cigarStrings = cigarStrings;
            loc.score = output.isEmpty() ? 0 :output.get(0).score;
            return loc;
        }
    }

    /*
     * Filters the potential matches we found using the mini seeds by checking if multiple seeds
     * align to similar parts of the genome (max distance here is 15 but can be changed depending on compute)
     */
    private ArrayList<Integer> filterCandidates(ArrayList<ArrayList<Integer>> potentialMatches){
       HashMap<Integer,Integer> candidateCounts = new HashMap<>();
       for(int i = 0; i < potentialMatches.size(); i++){
        ArrayList<Integer> matches = potentialMatches.get(i);
        int seedOffset = i*25;

        for(int match : matches){
            int idx = match - seedOffset;
            candidateCounts.put(idx, candidateCounts.getOrDefault(idx, 0) + 1);
            for(int j = 1; j <= 15; j++){
                if(candidateCounts.containsKey(idx-j)){
                    candidateCounts.put(idx-j, candidateCounts.getOrDefault(idx-j, 0) + 1);
                }
            }
        }
       }
       int maxSupport = 0;
        for (int support : candidateCounts.values()) {
            maxSupport = Math.max(maxSupport, support);
        }
       ArrayList<Integer> filteredCandidates = new ArrayList<>();
       for(var entry: candidateCounts.entrySet()){
        if(entry.getValue() == maxSupport){
            filteredCandidates.add(entry.getKey());
        }
       }
       return filteredCandidates;
    }

    private class ProblemOutput{
        int score;
        int y_start;
        int y_end;
        String cigar;

        private ProblemOutput(int score, int y_start, int y_end, String cigar){
            this.score = score;
            this.y_start = y_start;
            this.y_end = y_end;
            this.cigar = cigar;
        }
    }

    /*
     * Aligns queries with potential matches using 2D DP algorithm and builds a CIGAR string 
     */
    private ArrayList<ProblemOutput> alignCandidates(String query, ArrayList<Integer> candidates){

        int min_score = Integer.MAX_VALUE;
        ArrayList<ProblemOutput> outputList = new ArrayList<>();
        for(int i : candidates){
            int start = Math.max(0,i-15);
            int end = Math.min(this.genome.length(), i+query.length()+15);
            String reference = this.genome.substring(start,end);
            ProblemOutput output = this.fittingAlignment(query, reference, mismatch_penalty, gap_penalty,start,end);
            if(output.score < min_score){
                min_score = output.score;
                outputList.clear();
                outputList.add(output);
            }
            else if(output.score == min_score){
                outputList.add(output);
            }
        }
        return outputList;
    }
   
    private ProblemOutput fittingAlignment(String query_x, String reference_y, int mismatch_penalty, int gap_penalty, int start_genome, int end_genome){
        int cols = query_x.length()+1;
        int rows = reference_y.length()+1;
        int[][] dp = new int[rows][cols];

        for(int i = 0; i < cols; i++){dp[rows-1][i] = i*gap_penalty;}
        
        for(int r = rows-2; r >= 0; r--){
            for(int c = 1; c < cols; c++){
                char x = query_x.charAt(c-1);
                char y = reference_y.charAt((rows-2)-r);
                int cost = (x == y) ? 0 : mismatch_penalty;
                int val = Integer.max(gap_penalty+dp[r+1][c],gap_penalty+dp[r][c-1]);
                val = Integer.max(val,cost+dp[r+1][c-1]);
                dp[r][c] = val;
            }
        }

        int score = dp[0][cols-1];
        int start_row = 0; int end_row;
        for(int i = 0; i < rows; i++){
            if(dp[i][cols-1] >= score){
                score = dp[i][cols-1];
                start_row = i;
            }
        }

        StringBuilder cigar = new StringBuilder();
        int r = start_row; int c = cols-1;
        while(r < rows-1 && c > 0){
            char x = query_x.charAt(c-1);
            char y = reference_y.charAt(rows-2-r);
            int cost = (x==y)?0:mismatch_penalty;
            if(dp[r][c] == cost+dp[r+1][c-1]){
                cigar.append((x == y) ? '=' : 'X');
                r++; c--;
            }
            else if(dp[r][c] == dp[r+1][c]+gap_penalty){
                cigar.append('D');
                r++; 
            }
            else{
                cigar.append('I');
                c--;
            }
        }
        while(c > 0){
            cigar.append('I');
            c--;
        }
        end_row = r;
        String extended = cigar.reverse().toString();
        
        return new ProblemOutput(score,start_genome+reference_y.length()-end_row,start_genome+reference_y.length()-start_row,this.compressString(extended));
    }

    private String compressString(String str){
        StringBuilder compressed = new StringBuilder();
        int count = 1;

        for(int i = 1; i < str.length(); i++){
            if(str.charAt(i) == str.charAt(i-1)){
                count += 1;
            }
            else{
                compressed.append(count).append(str.charAt(i - 1));
                count = 1;
            }
        }
        compressed.append(count).append(str.charAt(str.length() - 1));
        return compressed.toString();
    }

    /*
     * Runs the program for every query and writes the output out to a new file 
     */
    private void runQueries(){
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(this.outputFile));
            for(FastaRecord queryRecord : this.readsFile.records){
                Locations loc = this.potentialCandidates(queryRecord.seq);
                if(loc.exactMatch){
                    writer.write(queryRecord.header + "\t" + loc.numLocations);
                    writer.newLine();
                    for(int i : loc.matchLocations){
                        writer.write(i + "\t" + 0 + "\t" + queryRecord.seq.length() + "=");
                        writer.newLine();
                    }
                }else{
                    writer.write(queryRecord.header + "\t" + loc.numLocations);
                    writer.newLine();
                    for(int i = 0; i < loc.numLocations; i++){
                        writer.write(loc.matchLocations.get(i) + "\t" + loc.score + "\t" + loc.cigarStrings.get(i));
                        writer.newLine();
                    }
                }
            }
            writer.close();
        } catch (IOException e) {
        }
    }

    
    public static void main(String[] args) {
        String index_file = args[0];
        String read_file = args[1];
        int mismatch_penalty = Integer.parseInt(args[2])*-1;
        int gap_penalty = Integer.parseInt(args[3])*-1;
        String output_file = args[4];

        Picomap picomap = new Picomap(index_file, read_file, mismatch_penalty, gap_penalty, output_file);
        picomap.runQueries();
    }
}
