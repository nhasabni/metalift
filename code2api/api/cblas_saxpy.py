from metalift.analysis import CodeInfo
from metalift.ir import *

from .utils import filter_types

#
# saxpy (a, x : list, y : list) = 
#   z : list = a * x + y
#
def getGrammar(ci: CodeInfo):
    name = ci.name

    int_vars = filter_types(ci.readVars, Int())
    list_vars = filter_types(ci.readVars, ListT(Int()))
    x_y = Choose(*list_vars)

    if name.startswith("inv"):  #inv
      z = ci.modifiedVars[0] # z is out also
      i = ci.modifiedVars[1]

      intChoices = Choose(IntLit(0), IntLit(1))
    
      # i >= 0
      i_lower_bound = Choose(Ge(i, intChoices), Gt(i, intChoices),
                           Le(i, intChoices), Lt(i, intChoices),
                           Eq(i, intChoices))
      # i <= listlen
      i_upper_bound = Choose(Ge(i, Call("list_length", Int(), x_y)),
                           Gt(i, Call("list_length", Int(), x_y)),
                           Le(i, Call("list_length", Int(), x_y)),
                           Lt(i, Call("list_length", Int(), x_y)),
                           Eq(i, Call("list_length", Int(), x_y)))
      saxpy_output = Choose(Eq(z, Call("cblas_saxpy", ListT(Int()), Choose(*int_vars),
                                       Call("list_take", ListT(Int()), x_y, i),
                                       Call("list_take", ListT(Int()), x_y, i))))
      inv = Choose(Implies(Eq(Call("list_length", Int(), x_y), Call("list_length", Int(), x_y)),
                           And(saxpy_output, And(i_lower_bound, i_upper_bound))))
      return inv

    else:  #ps
      z = ci.modifiedVars[0]
      choices = Choose(Implies(Eq(Call("list_length", Int(), x_y), Call("list_length", Int(), x_y)),
                               Eq(z, Call("cblas_saxpy", ListT(Int()), Choose(*int_vars), Choose(*list_vars), Choose(*list_vars)))))
      return choices

#
# cblas_saxpy(a, x, y) =
#   return Concat(a * x[0] + y[0]), cblas_saxpy(a, x[1:], y[1:])
#
# TODO: Support case of incx, incy (look at stored patch)
def getTargetLang():
    a = Var("a", Int())
    x = Var("x", ListT(Int()))
    y = Var("y", ListT(Int()))
    cblas_saxpy = FnDeclRecursive(
        "cblas_saxpy",
        ListT(Int()),
        Ite(And(And(
                  Gt(Call("list_length", Int(), x), IntLit(0)),
                  Gt(Call("list_length", Int(), y), IntLit(0))),
                Eq(Call("list_length", Int(), x), Call("list_length", Int(), y))),
            Call("list_prepend", ListT(Int()),
              Add(Mul(a, Call("list_get", Int(), x, IntLit(0))), 
                  Call("list_get", Int(), y, IntLit(0))),
              Call("cblas_saxpy", ListT(Int()),
                    a, Call("list_tail", ListT(Int()), x, IntLit(1)),
                    Call("list_tail", ListT(Int()), y, IntLit(1)))),
            Call("list_empty", ListT(Int()))
        ),
        a, x, y
    )
    return [cblas_saxpy]