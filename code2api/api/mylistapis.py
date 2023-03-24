from metalift.analysis import CodeInfo
from metalift.ir import *

def getGrammar(ci: CodeInfo):
    name = ci.name 
    
    if name.startswith("inv"):
        if len(ci.readVars) != 1 or len(ci.modifiedVars) != 2:
             return BoolLit(False)

        a = ci.readVars[0]
        ret = ci.modifiedVars[0]
        i = ci.modifiedVars[1]

        if parseTypeRef(a.type) != Type("MLList", Int()) or parseTypeRef(ret.type) != Int():
            return BoolLit(False)
        
        i_in_bounds = And(Ge(i, IntLit(0)), Le(i, Call("list_length", Int(), a)))
        list_len = Ge(Call("list_length", Int(), a), IntLit(0))
        sum_output = Eq(ret, Call("mylistsum", Int(), 
                                  Call("list_take", ListT(Int()), a, i)))
        len_output = Eq(ret, Call("mylistlen", Int(),
                                  Call("list_take", ListT(Int()), a, i)))
        api_choice = Choose(sum_output, len_output)
        inv = And(api_choice, And(i_in_bounds, list_len))
        return inv
    else:
        if len(ci.modifiedVars) < 1:
            return BoolLit(False)
        
        ret = ci.modifiedVars[0]
        if parseTypeRef(ret.type) != Int():
            return BoolLit(False)

        possible_apis = Choose(
            Call("mylistlen", Int(), *ci.readVars),
            Call("mylistsum", Int(), *ci.readVars)
        )
        summary = Eq(ret, possible_apis)
        return summary
    
#
# Synthesize 2 APIs - 
#
# mylistlen(x : list) =
#   if len(x) > 0:
#       1 + mylistlen(x[1:])
#   else
#       0
#
# mylistsum(x: list) =
#   if len(x) > 0:
#       x[0] + mylistsum(x[1:])
#   else:
#       0
#
def getTargetLang():
    a = Var("a", ListT(Int()))
    mylistlen = FnDeclRecursive(
        "mylistlen",
        Int(),
        Ite(
            Gt(Call("list_length", Int(), a), IntLit(0)),
            Add(IntLit(1), Call("mylistlen", Int(),
                                Call("list_tail", ListT(Int()), a, IntLit(1)))),
            IntLit(0)
        ),
        a,
    )
    mylistsum = FnDeclRecursive(
        "mylistsum",
        Int(),
        Ite(
            Gt(Call("list_length", Int(), a), IntLit(0)),
            Add(Call("list_get", Int(), a, IntLit(0)), 
                Call("mylistsum", Int(),
                    Call("list_tail", ListT(Int()), a, IntLit(1)))),
            IntLit(0)
        ),
        a,
    )
    return [mylistlen, mylistsum]
