import sys

chr_dict = {"chrX":"NC_003284.9", "chrM":"NC_001328.1", "chrI":"NC_003279.8", "chrII":"NC_003280.10" , "chrIII":"NC_003281.10", "chrIV":"NC_003282.8", "chrV":"NC_003283.11"} #dictionary to rename the .bed12 file

downloaded_bed = "/media/rik/Blast/Course6/C.Elegans2.bed12"
new_bed = "/media/rik/Blast/Course6/test"

new_bed = open(new_bed, "w+")

with open(downloaded_bed, "r") as f:
	for line in f:
		split_line = line.strip().split("\t")
		print(split_line)
		if split_line:
			chr_name = split_line[0] #which column to chose

			new_bed.write(line.replace(chr_name, chr_dict.get(chr_name))) #replace that with the new name
new_bed.close()
