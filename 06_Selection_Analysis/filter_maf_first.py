#!/usr/bin/env python3
import sys

def filter_maf(maf_file, out_file):
    with open(maf_file) as fin, open(out_file, 'w') as fout:
        # first copy the header (until the first a)
        for line in fin:
            if line.startswith('a'):
                # start processing the first block, put this line into the block buffer
                block = [line]
                break
            else:
                fout.write(line)
        else:
            return  # no a line

        # process all blocks
        for line in fin:
            if line.startswith('a'):
                # process the current block
                process_block(block, fout)
                block = [line]
            else:
                block.append(line)
        # process the last block
        if block:
            process_block(block, fout)

def process_block(block, fout):
    seen = set()
    out_block = []
    for line in block:
        if line.startswith('s'):
            parts = line.split()
            if len(parts) < 2:
                out_block.append(line)
                continue
            full_name = parts[1]
            # extract species identifier (before the dot)
            species = full_name.split('.')[0]
            if species not in seen:
                seen.add(species)
                parts[1] = species
                out_block.append('\t'.join(parts) + '\n')
            # otherwise skip
        else:
            out_block.append(line)
    fout.write(''.join(out_block))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: filter_maf_first.py input.maf output.maf")
        sys.exit(1)
    filter_maf(sys.argv[1], sys.argv[2])
