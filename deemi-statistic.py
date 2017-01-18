from pysam import AlignmentFile

def pad_softclip(filename):
    namedict, posdict = {}, {}

    for read in AlignmentFile(filename).fetch():
        adjusted_start = read.reference_start - read.query_alignment_start
        adjusted_end   = adjusted_start + read.infer_query_length()
        name           = read.query_name
        reference_name = read.reference_name

        try_append(namedict, name, (reference_name, adjusted_start, adjusted_end))

    for k, v in namedict.items():
        if len(v) == 1:
            ref, start, end = v[0]
            try_add(posdict, (ref, start, end - start + 1), 0, 1)
        elif len(v) == 2:
            if v[0][0] != v[1][0]:
                print("what the fuck {}:{}".format(k, v))
                continue
            start  = min(map(lambda x: x[1], v))
            length = max(map(lambda x: x[2], v)) - start
            try_add(posdict, (v[0][0], start, length), 1, 0)
        else:
            print("what the fuck {}:{}".format(k, v))

    return posdict

def try_append(d, k, v):
    if k in d:
        d[k].append(v)
    else:
        d[k] = [v]

def try_add(d, k, a, b):
    x, y = d.get(k, (0, 0))
    d[k] = x+a, y+b

length_dict = {}

for ((ref, start, length), (pair, single)) in pad_softclip("S006_cfdna_rg.bam").items():
    if length in length_dict:
        length_dict[length] += pair
    else:
        length_dict[length] = pair
    if pair > 1:
        print("In {}:{}, there are {} pairs have length {}".format(ref, start, pair, length))
    if single > 1:
        print("In {}:{}, there are {} single reads have length {}".format(ref, start, single, length))

for length, count in length_dict.items():
    print("{} pairs total have length {}".format(count, length))
