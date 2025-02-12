import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/*
 * FASTA file parser used to parse the reference genome provided
 * Custom class object that is used to later create either a Suffix Array or an FM Index
 */
public class FastaFile {
    String filePath;
    ArrayList<FastaRecord> records;
    int no_records;

    public FastaFile(String path) {
        this.filePath = path;
        this.records = new ArrayList<>();
        this.no_records = 0;
        this.parseFile();
      }

    private void parseFile(){
        try {
            File file = new File(this.filePath);
            BufferedReader reader = new BufferedReader(new FileReader(file));

            String line;
            StringBuilder curr_str = null;
            String header = null;

            while ((line = reader.readLine()) != null){
                if (line.startsWith(">")){
                    if(curr_str != null){
                        FastaRecord curr_record = new FastaRecord(header,curr_str.toString());
                        this.records.add(curr_record);
                        this.no_records += 1;
                    }
                    curr_str = new StringBuilder();
                    header = line.substring(1);
                } else{
                    curr_str.append(line.trim());
                }
            }

            if(curr_str != null){
                FastaRecord curr_record = new FastaRecord(header,curr_str.toString());
                this.records.add(curr_record);
                this.no_records += 1;
            }

            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

/*
 * Each FASTA file can have multiple records which are stored individually as FASTARecord objects
 */

class FastaRecord{
    String header;
    String seq;
    int seq_len;

    public FastaRecord(String header, String seq){
        this.header = header;
        this.seq = seq;
        this.seq_len = seq.length();
    }
}