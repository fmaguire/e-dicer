#!/usr/bin/env python3

# Script to get sample

import sys
import random
import wget

def get_file(path):
    """
    Ensembl occasionally capitalise the species name at the end for reasons
    """

    try:
        retcode = wget.download(path)
    except IOError:
        try:
            print("Capitalisation error")
            new_fn = path.split('/')[-1]
            first_char = new_fn[0].upper()
            new_fn = first_char + new_fn[1:]
            new_path = "/".join(path.split("/")[:-1] + [new_fn])
            retcode = wget.download(new_path)
        except IOError:
            retcode = -1

    return retcode

def subsample_dataset(dataset):
    """
    Get 100 random samples from the dataset from different species i.e.
    don't sample 2 strains together
    """

    species_db = {}
    for entry in dataset:
        species_name = " ".join(entry[0].split(" ")[:2])

        if species_name in species_db:
            species_db[species_name].append(entry)
        else:
            species_db.update({species_name: [entry]})

    subsample = random.sample(species_db.values(), 150)

    subsample_dataset = []

    for entry in subsample:

        if len(entry) > 1:
            entry  = random.sample(entry, 1)

        value = entry[0]
        subsample_dataset.append(value)


    return subsample_dataset


if __name__=='__main__':

    with open(sys.argv[1]) as fh:
        dataset = [line.strip().split('\t') for line in fh]

    dataset = subsample_dataset(dataset)

    counter = 0
    for line in dataset:

        collection = "_".join(line[12].split('_')[:3])

        path = "ftp://ftp.ensemblgenomes.org/pub/bacteria/current/fasta/" + collection + '/' + line[1] + '/cdna/' + line[1] + '.' + line[4] + '.' + 'cdna.all.fa.gz'

        retcode = get_file(path)

        if retcode != -1:
            counter += 1

        if counter == 99:
            sys.exit(0)
