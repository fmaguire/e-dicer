#!/usr/bin/env python

import sys

# get top hits

def parse_top_hits(db, collisions):
    """
    Sort collisions by expression level in db
    """
    expression = {}

    with open(db) as fh:
        next(fh)
        for line in fh:
            line = line.strip().split('\t')
            expression.update({line[0]: float(line[13])})

    collision_data = []
    with open(collisions) as fh:
        for line in fh:
            line = line.strip().split('\t')
            host = line[0].split()[0]
            seq = line[4]
            hit = line[2]
            collision_data.append([host, seq, hit])

    for collision in collision_data:
        collision.append(expression[collision[0]])


    collision_data.sort(key = lambda x: x[3], reverse=True)

    top100 = collision_data[:100]

    with open(collisions + "_top100.tsv", 'w') as fh:
        fh.write("Host_Transcript\tCollision_Sequence\tCollision_Accession\tHost_Expression_(Kallisto_TMM_Overall)\n")
        for hit in top100:
            hit = map(str, hit)
            fh.write("\t".join(hit) + '\n')


    already_seen = []
    unique_collisions = []
    for collision in collision_data:
        if collision[0] not in already_seen:
            unique_collisions.append(collision)
            already_seen.append(collision[0])

    with open(collisions + "_unique_top100.tsv", 'w') as fh:
        fh.write("Host_Transcript\tCollision_Sequence\tCollision_Accession\tHost_Expression_(Kallisto_TMM_Overall)\n")
        for hit in unique_collisions[:100]:
            hit = map(str, hit)
            fh.write("\t".join(hit) + '\n')

if __name__=='__main__':

    # first file is db with expression
    # second is bowtie file with collisions hits


    parse_top_hits(sys.argv[1], sys.argv[2])

