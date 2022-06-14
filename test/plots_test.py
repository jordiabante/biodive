from biodive import plots

infile = "./test_table_dist_1_ksize_30.txt.gz"
outfile1 = "./test_table_dist_1_ksize_30_hist.png"
outfile2 = "./test_table_dist_1_ksize_30_hist_upvsdown.png"
plots.plt_histo(infile,outfile1)
plots.plt_histo_up_vs_down(infile,outfile2)
