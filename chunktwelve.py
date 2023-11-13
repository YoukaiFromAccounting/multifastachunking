from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os

#Take in user input for algorithm parameters
if len(sys.argv) != 4:
    print("Incorrect script call, ensure call is formatted as: python3 chunketc.py <input> <output directory> <chunksize>")
    sys.exit(1)
  
#Read in name of the file and the desired size of chunks in nucleotides
input_file = sys.argv[1]
output_directory = sys.argv[2]
chunk_size = int(sys.argv[3])

#Use SeqIO.parse instead of read for multifasta file IMPORTANT
records = list(SeqIO.parse(input_file, "fasta"))

def create_batch(records, chunk_size):
    record_it = iter(records)

    record = next(record_it)
    current_base = 0

    batch = []
    batch_size = 0

    # While there are new records, keep creating new batches.
    while record:
        # Loop over records untill the batch is full. (or no new records)
        while batch_size != chunk_size and record:

            end = current_base + chunk_size - batch_size
            seq = record[current_base:end]

            end_of_slice = current_base + len(seq) - 1
            #Use record.id for pure name, record.description for contig details
            fasta_header = record.description + ":{}-{}".format(current_base, end_of_slice)
            
            seq.id = seq.name = fasta_header
            seq.description = ''
            batch.append(seq)

            current_base += len(seq)
            batch_size += len(seq)

            # Current record is exhausted, get a new one.
            if current_base >= len(record):
                record = next(record_it, None)
                current_base = 0

        # We have a batch with the correct size
        yield batch
        batch = []
        batch_size = 0

#Ensure output directory exists
os.makedirs(output_directory, exist_ok = True)

#Range over batches and export each chunk as a fasta file
for i, batch in enumerate(create_batch(records, chunk_size)):
    #Add input_file to prevent file overwrite from chunk loading
    filename = os.path.join(output_directory,(input_file+"chunk{}.fasta").format(i))
    SeqIO.write(batch, filename, "fasta")