from metalift.ir import *

def filter_types(inp, inp_type):
    return list(filter(lambda x: parseTypeRef(x.type) == inp_type, inp))