args <- commandArgs(trailingOnly = TRUE)
#load data
Input <- read.table(args[1], sep = "\t", header = T)


#Log scale transformation of data
output_result <- log(Input+1)


# Write data in a file
write.table(output_result,file = args[2], sep = "\t",  row.names=F, quote = F)

# Command to run the code:
# Rscript script.R input_data.txt output.txt