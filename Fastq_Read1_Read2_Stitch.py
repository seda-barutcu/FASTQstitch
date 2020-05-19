# 1: open the text file that has the list of all read1 file names (do not include read2 file names!). Open each read1 file one by one.

with open('filelistR1.txt') as f:
    for line in f:
        read1file = line.strip()
        
        # 2: Create ReadIdentifier-Sequence dictionary and ReadIdentifier-QualityScore dictionary for each read1-fastq file.
        
        ident_seq_dict = dict()
        ident_qc_dict = dict()
        it=2
        with open(read1file) as f1:
            for line in f1:
                info = line.strip()
                it = it+1
                if it%4==3:
                    infol = info.split()
                    ident = infol[0]
                elif it%4==0:
                    seq = info
                    ident_seq_dict[ident] = seq
                elif it%4==2:
                    qc = info
                    ident_qc_dict[ident] = qc
                    
        # 3: For each read1 file, open paired read2 file.
        # First line below replaces the last 12 character of the read1-filename with the last 12 characters of read2-filename.
        # If necessary, this line should be modified according to the file name format.
       
        read2file = read1file[:-12] + "R2_001.fastq"
        it=2
        with open(read2file) as f2:
            for line in f2:
                info = line.strip()
                it = it+1
                if it%4==3:
                    infol = info.split()
                    ident = infol[0]
                elif it%4==0:
                    seq = info
                    
                    # 4: convert each read2 sequence that matches a read1-identifier to its reverse-complement.
                    
                    if ident in ident_seq_dict :
                        letterdic = dict()
                        letn = 0
                        for letter in seq :
                            letn = letn -1
                            letterdic[letn] = seq[letn]
                        dicn = 0
                        for dicletter in letterdic :
                            dicn = dicn - 1
                            if letterdic[dicn] == "A" : letterdic[dicn] = "T"
                            elif letterdic[dicn] == "T" : letterdic[dicn] = "A"
                            elif letterdic[dicn] == "G" : letterdic[dicn] = "C"
                            elif letterdic[dicn] == "C" : letterdic[dicn] = "G"
                            else : letterdic[dicn] = letterdic[dicn]
                        nseq = ""
                        nseqn = -1
                        for nuc in letterdic :
                            nseq = nseq + letterdic[nseqn]
                            nseqn = nseqn-1
                            
                        # 5: If the last 12nt of read1-sequence is found in rev-comp of read2-sequence:
                        # add the read2-sequence to read1-sequence, excluding redundent parts.
                        
                        seq1 = ident_seq_dict[ident]
                        if seq1[-12:] in nseq :
                            key = seq1[-12:]
                            num = nseq.find(key)
                            pos = num+12
                            ident_seq_dict[ident] = ident_seq_dict[ident] + nseq[pos:]
                        
                        #if can not stitch, write read1/2 sequences separated by ">"
                        
                        else : ident_seq_dict[ident] = ident_seq_dict[ident] + ">" + nseq
                
                # 6: stitch reverse of read2 quality-score line to read1-qc-line in accordance with sequence line stitching:
                
                elif it%4==2 :
                    qc = info
                    qcdic = dict()
                    qcn = 0
                    for letter in qc :
                        qcn = qcn -1
                        qcdic[qcn] = qc[qcn]
                    nqc = ""
                    nqcn = -1
                    for let in qcdic :
                        nqc = nqc + qcdic[nqcn]
                        nseqn = nseqn-1
                    ident_qc_dict[ident] = ident_qc_dict[ident] + nqc[pos:]
        
        # 7: write only stitched reads to output file in fastq format.
        
        outputfilename = read1file[8:18] + "stitsched.fastq"
        fileout = open(outputfilename, "w")
        for ids in ident_seq_dict:
            if ">" not in ident_seq_dict[ids] :
                fileout.write(ids + "\n" + ident_seq_dict[ids] + "\n"+ "+" + "\n" + ident_qc_dict[ids] + "\n")
        
        # Optional: un-stitched reads can be written to a separate file.
