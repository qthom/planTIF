#maxim.ivanov@plen.ku.dk
#Modified by Quentin Thomas
#qt@plen.ku.dk
#latest update = 06/11/2017


chroms = ['1', '2', '3', '4', '5'] # have to match chromosome names in BAM file

# report only isoforms with at least this number of supporting reads

import pysam, collections, sys
from operator import itemgetter

bam_filename = sys.argv[1] # bam file is expected to be coordinate-sorted (not name-sorted!)
bg_filename = bam_filename[:bam_filename.rfind('.')] + '.bed'

offset = 1
offset2 = 20
min_support = 3

#function that checks if the bam file is sorted:

def checkDeque(pointer, clusters):
    while clusters and clusters[-1][2] < (pointer - offset): # this and further reads have no chances to associate with the oldest cluster(s)
        results.append(clusters.pop()) # move this cluster to results

def checkDeque2(pointer, clusters):
    while clusters and clusters[-1][1] < (pointer - offset2): # this and further reads have no chances to associate with the oldest cluster(s)
        Final_clusters.append(clusters.pop()) # move this cluster to results

#left = read start
#right = mate alignment start
#chrom = chromosome
#strand = sens of strand for the read
#the main function parse the previous clusters and if the reads analyzed is corresponding to previous clusters it is a "hit" then hi=true
#if the hit remain false then the read does not match any other cluster then it creates a cluster for that specific read.

def screenClusters(clusters, left, right, chrom, strand):
    hit = False
    for i in range(len(clusters)):
        cluster = clusters[i] #create a cluster with the first current analyzed read
        min_left, max_left, min_right, max_right, value, chrom, strand = cluster
        if left >= (min_left - offset) and left <= (max_left + offset) and right >= (min_right - offset) and right <= (max_right + offset):
            hit = True # read can be associated with given cluster
            #this below will uptdate the hit cluster with the minimum value of the previous value and the current value, ad the maximum of the maximun and the current value for both ends
            updated_cluster = [min(left, min_left), max(left, max_left), min(right, min_right), max(right, max_right), (value + 1), chrom, strand]            
            clusters[i] = updated_cluster # this update this cluster with the new values            
            break # more iterations through clusters are not required
    if hit == False: # if read does not match any cluster, or deque was empty
        #creates a new cluster with the current read coordinates       
        new_cluster = [left, left, right, right, 1, chrom, strand]      
        clusters.appendleft(new_cluster) #add the new cluster at the left of the cluster list to increase parsing rapidity

#the following function will merge previousy identified cluster which larger bondaries overlap with other clusters of a 9bp window

def Merge_Clusters(clusters, Start, End, Value, chrom, strand):
    hit = False  
    for i in range(len(clusters)):
        cluster = clusters[i] #create a cluster with the first current analyzed read
        min_left, max_right, value_new, chrom, strand = cluster
        if Start <= (min_left + offset2) and Start >= (min_left - offset2) and End >= (max_right - offset2) and End <= (max_right + offset2):
            hit = True # read can be associated with given cluster
            #this below will uptdate the hit cluster with the minimum value of the previous value and the current value, ad the maximum of the maximun and the current value for both end
            updated_cluster = [min(Start, min_left), max(End, max_right), (Value + value_new), chrom, strand]             
            clusters[i] = updated_cluster # this update this cluster with the new values               
            break # more iterations through clusters are not required
    if hit == False: # if read does not match any cluster, or deque was empty
        #creates a new cluster with the current read coordinates 
        new_cluster = [Start,End, Value, chrom, strand]        
        clusters.appendleft(new_cluster) #add the new cluster at the left of the cluster list to increase parsing rapidity
#....................................................................................................................................


results = []

with pysam.AlignmentFile(bam_filename, 'rb') as bamfile:
    print ('clustering:')
    #select the chromosome from 1 to 5 and run the functions over each chromosomal reads
    for chrom in chroms:
        print('\n', chrom, end='', flush = True)
        progress = 0
        fw_clusters, rev_clusters = collections.deque(), collections.deque()   #create list deque for the foward strand clusters and the reverse strand clusters
        fw_pointer, rev_pointer = 0, 0
        reads = bamfile.fetch(chrom) # iterate over all reads on only the given chrom
        for read in reads:       
            if read.is_paired and read.is_proper_pair and not read.is_reverse: # consider only forward reads because they are guaranteed to appear in increasing order (given that file is coordinate-sorted))
                progress += 1               
                if progress % 100000 == 0:          
                    print('.', end = '', flush = True)               
                left, right = read.reference_start, (read.reference_start + read.template_length) # leftmost and rightmost coordinates of the read pair                
                if read.is_read1: # if read is first-in-pair
                    strand = '+' # this is transcript orientation, not read orientation
                    fw_pointer = left
                    checkDeque(fw_pointer, fw_clusters) # check if any clusters can be moved to results
                    screenClusters(fw_clusters, left, right, chrom, strand) # check current read pair vs current clusters                
                elif read.is_read2:                  
                    strand = '-'
                    rev_pointer = left
                    checkDeque(rev_pointer, rev_clusters)
                    screenClusters(rev_clusters, left, right, chrom, strand)       
        # when reached chromosome end, dump current clusters to the results:
       
        results = results + list(fw_clusters) + list(rev_clusters) # convert deque to list for successful concatenation

#sort the results by chromosome and strand to fasten the iteration to merge clusters 

results = [[r[0], r[1], r[2], r[3], r[4], r[5], r[6]]  for r in results]
results = sorted(results, key = itemgetter(5, 0, 3, 6))


Final_clusters = []

for chrom in chroms:
    Cluster_of_Cluster_fw, Cluster_of_Cluster_rev = collections.deque(), collections.deque()
    for i in results:
        if i[5] == chrom:
            if i[6] == '+':
                Start, End, Value, chrom, strand = i[0] , i[3], i[4], i[5], i[6]
                fw_pointer = i[0]
                checkDeque2(rev_pointer, Cluster_of_Cluster_fw)
                Merge_Clusters(Cluster_of_Cluster_fw, Start, End, Value, chrom, strand)        
            else:
                Start, End, Value, chrom, strand = i[0] , i[3], i[4], i[5], i[6]
                rev_pointer = i[0]
                checkDeque2(rev_pointer, Cluster_of_Cluster_rev)           
                Merge_Clusters(Cluster_of_Cluster_rev, Start, End, Value, chrom, strand)
                
    Final_clusters = Final_clusters + list(Cluster_of_Cluster_fw) + list(Cluster_of_Cluster_rev)



# convert results to BED6 format (min_left and max_right coordinates are considered as cluster borders)



bed6_results = [[r[3], r[0], r[1], r[2], ".", r[4]]  for r in Final_clusters]
bed6_results = sorted(bed6_results, key = itemgetter(0, 1, 2))


#write a bedgraph file for the clusters:

with open(bg_filename, 'w') as outfile:
    countTIF = 0
    for result in bed6_results:
        if result[3] >= min_support:
            countTIF = countTIF + 1
            
            #un-hide the two following to print clusters name in the pseudo bed
            #TIF = 'TIF' +  str(countTIF)
            #result.insert(3, TIF)
            

            outfile.write('\t'.join([str(f) for f in result]) + '\n')

print ('Done')