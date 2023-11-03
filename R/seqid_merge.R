# by Sisi Liu: sisi.liu@awi.de; sisi.liu.research@gmail.com
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS 14.1
#==== merge seqids ====
seqid=list.files("seqid", pattern = "*_seqid.txt", full.names = T)
input_strings = c("1152:Pseudanabaena:genus", "1158:Oscillatoria:genus", "1177:Nostoc:genus", "3798:Saxifraga:genus", "4479:Poaceae:family",
                  "5748:Nannochloropsis:genus", "7953:Cyprinidae:family", "8015:Salmonidae:family", "9789:Equus:genus", "9859:Cervus:genus",
                  "9903:Bos:genus", "9977:Ochotona:genus", "13228:Potamogeton:genus", "13398:Carex:genus", "21880:Salvia:genus",
                  "23204:Potentilla:genus", "24952:Myriophyllum:genus", "40685:Salix:genus", "43174:Pedicularis:genus",
                  "47251:Leptolyngbya:genus", "54304:Planktothrix:genus", "102804:Asteroideae:subfamily", "167375:Cyanobium:genus",
                  "202994:Rhodiola:genus", "217161:Chamaesiphon:genus", "337677:Cricetidae:family")
input_strings=gsub(":", "_", input_strings)

for (i in input_strings) {
  file_list <- list.files(path = "seqid", pattern = paste0("^", i, "_.*_seqid.txt$"), full.names = TRUE)
  print(file_list)
  #assuming the same header/columns for all files
  datafr = do.call("rbind", lapply(file_list, function(x)read.table(x, header=F)) ) 
  write.table(datafr, paste0(i, ".txt"), quote = F, col.names = F, row.names = F)
}
#==== END ====