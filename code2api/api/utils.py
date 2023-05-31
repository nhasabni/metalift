from metalift.ir import *

def filter_types(inp, inp_type):
    def match(x, x_type):
        return parseTypeRef(x.type) == x_type
    return list(filter(lambda x: match(x, inp_type), inp))