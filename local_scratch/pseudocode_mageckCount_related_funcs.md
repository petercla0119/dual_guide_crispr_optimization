Pseudocode breakdown of mageckCount related functions
=====================================================

Common variables
	dict0 - holds sgRNA count data

# File: mageckCount.py
## mageckcount_parseargs
## mageckcount_checkargs
## mageckcount_mergedict
### Pseudocode:
	Function mageckcount_mergedict(dict0, dict1):
		Initialize nsample to 0
		
		If dict0 is not empty:
			Set nsample to the length of the first list in dict0
	
		For each key-value pair (k, v) in dict0:
			If key k exists in dict1:
				Append the value from dict1[k] to the list v in dict0
			Else:
				Append 0 to the list v in dict0
		
		For each key-value pair (k, v) in dict1:
			If key k does not exist in dict0:
				If nsample > 0:
					Create a new entry in dict0 with key k and a list of zeros of length nsample
				Else:
					Create a new entry in dict0 with key k and an empty list
				Append the value v from dict1 to this new list in dict0
	
		End Function (no return, modifies dict0 in place)

### Line-by-line:
## mageckcount_printdict
## mageck_printdict
## mageckcount_checklists
## mageckcount_processfastq
### Pseudocode:
	Function mageckcount_processfastq(args, genedict, sgdict):

    Initialize 'paired' to False
    
    If paired-end fastq files are provided (args.fastq_2 is not None):
        Set 'paired' to True
        Split 'args.fastq_2' into a list of paired fastq files ('pairedfq')

    Split 'args.fastq' into a list of fastq files ('listfq')
    Set 'nsample' to the number of samples (length of 'listfq')

    Initialize an empty dictionary 'datastat' for QC statistics

    Get all sample labels from 'args.sample_label'
    If 'args.sample_label' is empty:
        Create default sample labels ('sample1', 'sample2', ...)
    Else:
        Use the provided sample labels

    For each sample (indexed by 'i') in 'listfq':
        For each fastq file 'fi' in the sample:
            Initialize 'datastat[fi]' as an empty dictionary
            Set 'datastat[fi]["label"]' to the corresponding sample label

    Initialize an empty dictionary 'alldict' to store all merged sgRNA counts

    Set 'adjust' to True
    If paired-end reads should not be adjusted (args.count_pair is 'TRUE'):
        Set 'adjust' to False

    Initialize 'i' to 0

    For each list of filenames ('filenamelist') in 'listfq':
        Initialize an empty dictionary 'dict0' to store sgRNA counts for the current sample

        Initialize 'j' to 0

        For each fastq file 'filename' in 'filenamelist':
            If paired-end reads are used:
                Set 'pairedfile' to the corresponding paired-end file

            Initialize an empty dictionary 'dict00' to store counts from the current file

            If 'filename' ends with 'BAM':
                Process the BAM file with mageckcount_processonefile_bam
            Else If 'filename' ends with 'SAM':
                Process the SAM file with mageckcount_processonefile_sam
            Else:
                Process the fastq file with mageckcount_processonefile, considering 'pairedfile' and 'adjust'

            For each key-value pair (k, v) in 'dict00':
                If 'k' is not in 'dict0':
                    Initialize 'dict0[k]' to 0
                Add 'v' to 'dict0[k]'

            Increment 'j' by 1

        Increment 'i' by 1

        Merge 'dict0' into 'alldict' using mageckcount_mergedict

    Open the output file 'args.output_prefix + .count.txt' for writing ('ofilel')

    If 'unmapped_to_file' option is specified in args:
        Open the unmapped file 'args.output_prefix + .unmapped.txt' for writing ('ounmappedfilel')
    Else:
        Set 'ounmappedfilel' to None

    Print the contents of 'alldict' to 'ofilel' using mageckcount_printdict
    Close 'ofilel'

    If 'ounmappedfilel' was opened:
        Close 'ounmappedfilel'

    If 'sgdict' is not empty:
        Filter 'alldict' to only include known sgRNAs, storing the result in 'allmappeddict'
    Else:
        Set 'allmappeddict' to 'alldict'

    Return 'allmappeddict' and 'datastat'
### Line-by-line:

## getcounttablefromfile
## mageckcount_checkcontrolsgrna
## mageckcount_processcounttable
## mageckcount_main
# File: mageckCountNorm.py
## mageckcount_gettotalnormfactor
## mageckcount_getzerocounts
### Pseudocode:
### Line-by-line

## mageckcount_getmediannormfactor
### Pseudocode:
    # Step 1: Initialize variables
        n = length of the list corresponding to the first sgRNA in ctable  # Number of samples
        m = number of sgRNAs in ctable  # Number of sgRNAs
    
    # Step 2: Calculate the geometric mean of read counts for each sgRNA
    meanval = empty dictionary
        For each sgRNA (key k) and its associated read counts (list v) in ctable:
            If the total read count across all samples is greater than 0:
                Calculate the geometric mean of the read counts in v
                Store this mean in meanval with the key k
    
    # Step 3: Ensure geometric mean values are valid (no zeros)
        For each sgRNA (key k) and its geometric mean (value v) in meanval:
            If the geometric mean is 0 or negative, replace it with 1
    
    # Step 4: Initialize the list to store median normalization factors
        medianfactor = list of zeros with length n  # One factor for each sample
    
    # Step 5: Calculate the median normalization factor for each sample
        For each sample index ni from 0 to n-1:
            meanfactor = empty list
            For each sgRNA (key k) and its associated read counts (list v) in ctable:
                If k is in meanval:
                    Calculate the ratio of the read count at index ni to the geometric mean of k
                    Append this ratio to meanfactor
            
            # Calculate the median of the ratios
            xfactor = median of the sorted meanfactor list
            
            # If the median ratio is greater than 0, calculate the normalization factor
            If xfactor > 0:
                medianfactor[ni] = 1.0 / xfactor
        
    # Step 6: Return the list of median normalization factors
        Return medianfactor
### Line-by-line
#### def mageckcount_getmediannormfactor(ctable):
	# ctable: expected dictionary with sgRNA IDs as keys and lists of read counts as values
#### n = len(ctable[list(ctable.keys())[0]])  # samples
	# n = len(ctable[list(ctable.keys())[0]])': determines the number of samples (n). Assume each entry in ctable (each sgRNA) has a list of read counts with length == number of samples.
	# 'list(ctable.keys())[0]' retrieves the first key (sgRNA) in ctable, and len(ctable[first_key]) gives the number of samples
#### m = len(ctable)  # sgRNAs
	# 'm = len(ctable)':
		# Determines number of sgRNAs ('m') by calculating length of 'ctable' (# of keys in dict
#### meanval = {k: math.exp((sum([math.log(v2+1.0) for v2 in v])*1.0/n)) for (k, v) in ctable.items() if sum(v) > 0}  # geometric mean
	# meanval = {...}:
		# Calculate geometric mean of read counts for each sgRNA across all samples
		# k = sgRNA ID; v = list of read counts for that sgRNA
        # sum([math.log(v2+1.0) for v2 in v]):
			# calculate sum ln of the read counts plus 1 (to avoid taking the log of 0).
	# mathexp(...):
		# Calculate the sum of the logs then use exponential function (math.exp) to compute the geometric mean
	# if sum(v) > 0
		# Condition to ensure to only include sgRNAs with non-zero total counts
#### meanval = {k: (lambda x: x if x > 0 else 1)(v) for (k, v) in meanval.items()}  # delete those with all 0 read counts
	# meanval = {...}:
		# Adjust previously calculated geometric mean values
	# (lambda x: x if x>0 else 1)(v):
		# Avoid negative or zero mean values by replacing them with 1
    ## samplefactor=[0]*n
#### medianfactor = [0.0]*n
	# Initialize list of 0's w/ len = 'n' which stores median normalization factors for each sample
#### for ni in range(n): meanfactor = [v[ni]/meanval[k] for (k, v) in ctable.items() if k in meanval]
	# for ni in range(n):
		# Iterate over each sample index 'ni'
	# meanfactor = [...]:
		# for each sample (ni), calculate the ratio of the read count of each sgRNA (v[ni]) to its geometric mean (meanval[k])
	## print(str(sorted(meanfactor)))
#### xfactor = sorted(meanfactor)[len(meanfactor)//2]  # corrected
	# Calculate the median of 'meanfactor' list by sorting it then select the middle value
	# len(meanfactor)//2: give index of median in sorted list
#### if xfactor > 0.0: medianfactor[ni] = 1.0/xfactor
	# Check if median factor 'xfactor' is >0, 
		# if yes, then assign the inverse of 'xfactor' as the normalization factor for the current sample
	# medianfactor[ni] = 1.0/xfactor:
		# Assign the reciprocal of xfactor to medianfactor[ni] - used to normalize the read counts for the sample
            ## logging.debug('xfactor:'+str(xfactor))
#### return medianfactor
	# return list 'medianfactor' containing the median normalization factors for all samples

## normalizeCounts

# File: mageckMathFunc.py
## Pseudocode:
## Line-by-line


# mageckCountIO.py
## Pseudocode:
## Line-by-line


# mageckCountQC.py

/Users/Claire/Downloads/liulab-mageck-c491c3874dca/mageck/mageckCountNorm.py
