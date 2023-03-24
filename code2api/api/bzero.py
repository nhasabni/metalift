from metalift.analysis import CodeInfo
from metalift.ir import *

# bzero(a:list) =
#   if (len(a) > 0):
#     b[0] = 0
#     b[1:] = bzero(a[1:])
#   return b
#
# Intention is that the programmar will call it as:
#   out = bzero(in) # out is zeroed list of same size as in
#   in = out

def getGrammar(ci: CodeInfo):
    name = ci.name
       
    if name.startswith("inv"):  #inv
        if len(ci.readVars) != 1 or len(ci.modifiedVars) != 2:
            return BoolLit(False)
        
        a = ci.readVars[0]
        b = ci.modifiedVars[0] # b is out also
        i = ci.modifiedVars[1]

        if (parseTypeRef(a.type) != Type("MLList", Int()) or
            parseTypeRef(b.type) != Type("MLList", Int())):
            return BoolLit(False)
    
        i_in_bounds = And(Ge(i, IntLit(0)), Le(i, Call("list_length", Int(), a)))
        list_a_len = Ge(Call("list_length", Int(), a), IntLit(0))
        #list_len_eq = Eq(Call("list_length", Int(), a), Call("list_length", Int(), b))
        bzero_output = Eq(Call("list_take", ListT(Int()), b, i),
                        Call("bzero", ListT(Int()), Call("list_take", ListT(Int()), a, i)))
        api_choice = Choose(bzero_output)
        #inv = And(api_choice, And(i_in_bounds, And(list_a_len, list_len_eq)))
        inv = And(api_choice, And(i_in_bounds, list_a_len))
        return inv
    else:  #ps
        if len(ci.modifiedVars) == 0:
            return BoolLit(False)
        b = ci.modifiedVars[0] # b is out also
        if parseTypeRef(b.type) != Type("MLList", Int()):
            return BoolLit(False)
        possible_apis = Choose(Call("bzero", ListT(Int()), *ci.readVars))
        summary = Eq(b, possible_apis)
        return summary

def getTargetLang():
    a = Var("a", ListT(Int()))
    bzero = FnDeclRecursive(
        "bzero",
        ListT(Int()),
        Ite(
          Gt(Call("list_length", Int(), a), IntLit(0)),
          Call("list_concat", ListT(Int()), IntLit(0),
                Call("bzero", ListT(Int()), Call("list_tail", ListT(Int()), a, IntLit(1)))),
          Call("list_empty", ListT(Int()))
        ),
        a)
    return [bzero]