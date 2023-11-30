# Detailed Instructions for Converting Datasets Used in Artifact Evaluation

## TW Dataset from Konect

1. Download the dataset:
	```bash
	wget http://konect.cc/files/download.tsv.twitter.tar.bz2
	```
2. Unzip the dataset
	```bash
	tar -xvjf download.tsv.twitter.tar.bz2
	```
3. Remove any header/comments from the file (comments can be either % or #)
	```bash
	sed '/^[#%]/d' out.twitter > out.twitter.snap
	```
4. Follow the instructions from [Section 2.3.1 of the SVSIG_README](SVSIG_README.md#231-preparing-streaming-datasets-for-any-graph)


## TT Dataset from Konect

1. Download the dataset:
	```bash
	wget http://konect.cc/files/download.tsv.twitter_mpi.tar.bz2
	```
2. Unzip the dataset
	```bash
	tar -xvjf download.tsv.twitter_mpi.tar.bz2
	```
3. Remove any header/comments from the file (comments can be either % or #)
	```bash
	sed '/^[#%]/d' out.twitter_mpi > out.twitter_mpi.snap
	```
4. Follow the instructions from [Section 2.3.1 of the SVSIG_README](SVSIG_README.md#231-preparing-streaming-datasets-for-any-graph)



## FT Dataset from Konect

1. Download the dataset:
	```bash
	wget http://konect.cc/files/download.tsv.friendster.tar.bz2
	```
2. Unzip the dataset
	```bash
	tar -xvjf download.tsv.friendster.tar.bz2
	```
3. Remove any header/comments from the file (comments can be either % or #)
	```bash
	sed '/^[#%]/d' out.friendster > out.friendster.snap
	```
4. Follow the instructions from [Section 2.3.1 of the SVSIG_README](SVSIG_README.md#231-preparing-streaming-datasets-for-any-graph)

## UK Dataset from Konect

1. Download the dataset:
	```bash
	wget http://konect.cc/files/download.tsv.dimacs10-uk-2007-05.tar.bz2
	```
2. Unzip the dataset
	```bash
	tar -xvjf download.tsv.dimacs10-uk-2007-05.tar.bz2
	```
3. Remove any header/comments from the file (comments can be either % or #)
	```bash
	sed '/^[#%]/d' out.friendster > out.friendster.snap
	```
4. Follow the instructions from [Section 2.3.1 of the SVSIG_README](SVSIG_README.md#231-preparing-streaming-datasets-for-any-graph)