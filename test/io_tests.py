from biodive import io

# Test reading in table
infile = "./test_table_dist_1_ksize_30.txt.gz"
up,down = io.read_table(infile)
