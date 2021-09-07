def uniq(rlist):
    """
    """
    return list(set(rlist))

def cat(outfilename, *infilenames):
    """
    """
    with open(outfilename, 'w') as outfile:
        for infilename in infilenames:
            with open(infilename) as infile:
                for line in infile:
                    if line.strip():
                        outfile.write(line)

def cat_files(*infilenames):
    """
    """
    str_to_ret = ""
    for infilename in infilenames:
        with open(infilename) as infile:
            str_to_ret+=infile.read()
    return str_to_ret
                        
def cat_list(*infilenames):
    """
    """
    list_to_ret = []
    for infilename in infilenames:
        with open(infilename) as infile:
            for line in infile.read().split("\n"):
                list_to_ret.append(line)
    return list_to_ret