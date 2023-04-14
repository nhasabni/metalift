from metalift.analysis import CodeInfo
from metalift.ir import *

from .utils import filter_types

# 
# cblas_sgemv is implemented using smaller functions.
#
# sdot (x:n, z:n) -> ret:int =
#   sum(x[0] * z[0], sdot(x[1:], z[1:]))
#
# sgemv(a:mxn, x:n) -> ret:m = 
#   concat(sdot(a[0], x), sgemv(a[1:], x))
#
# sgemv_result = sgemv(a, x)
#
# cblas_sgemv(alpha, sgemv_result, beta, y) -> ret : m =
#   concat(alpha * sgemv_result[0] + beta * y[0],
#          cblas_sgemv(alpha, sgemv_result[1:], beta, y[1:]))
#
def getGrammar(ci: CodeInfo):
    name = ci.name

    int_vars = filter_types(ci.readVars, Int())
    list_vars = filter_types(ci.readVars, ListT(Int()))
    matrix_vars = filter_types(ci.readVars, ListT(ListT(Int())))

    x_y = Choose(*list_vars)
    a = Choose(*matrix_vars)

    i_j = Choose(filter_types(ci.modifiedVars, Int()))
    z = ci.modifiedVars[0] # z is out also

    if name.startswith("inv"):  #inv
      intChoices = Choose(IntLit(0), IntLit(1))
    
      # i >= 0
      i_lower_bound = Choose(Ge(i_j, intChoices), Gt(i_j, intChoices),
                             Le(i_j, intChoices), Lt(i_j, intChoices),
                             Eq(i_j, intChoices))
      # i <= listlen
      i_upper_bound = Choose(Ge(i_j, Call("list_length", Int(), x_y)),
                             Gt(i_j, Call("list_length", Int(), x_y)),
                             Le(i_j, Call("list_length", Int(), x_y)),
                             Lt(i_j, Call("list_length", Int(), x_y)),
                             Eq(i_j, Call("list_length", Int(), x_y)))
      # a:m_n * x:n_1 = sgemv_output : m_1
      sgemv_output = Call("sgemv", ListT(Int()), Call("list_get", ListT(Int()), a, i_j), x_y)
      # cblas_sgemv : m_1 = alpha * sgemv_output + beta * y
      cblas_sgemv_output = Choose(Eq(z, Call("cblas_sgemv", ListT(Int()), Choose(*int_vars),
                                       Call("list_take", ListT(Int()), sgemv_output, i_j),
                                       Choose(*int_vars),
                                       Call("list_take", ListT(Int()), x_y, i_j))))
      inv = Choose(Implies(And(
                                Eq(Call("list_length", Int(), x_y),
                                   Call("list_length", Int(), Call("list_get", ListT(Int()), a, IntLit(0)))),
                                Eq(Call("list_length", Int(), x_y), Call("list_length", Int(), a))
                              ),
                           And(cblas_sgemv_output, And(i_lower_bound, i_upper_bound))))
      return inv
    else:  #ps
      sgemv_output = Call("sgemv", ListT(Int()), Call("list_get", ListT(Int()), a, i_j), x_y)
      choices = Choose(Implies(And(
                                    Eq(Call("list_length", Int(), x_y),
                                       Call("list_length", Int(), Call("list_get", ListT(Int()), a, IntLit(0)))),
                                    Eq(Call("list_length", Int(), x_y), Call("list_length", Int(), a))
                                  ),
                               Eq(z, Call("cblas_sgemv", ListT(Int()),
                                          Choose(*int_vars), sgemv_output, Choose(*int_vars), x_y))))
      return choices

# 
# cblas_sgemv is implemented using smaller functions.
#
# sdot (x:n, z:n) -> ret:int =
#   sum(x[0] * z[0], sdot(x[1:], z[1:]))
#
# sgemv(a:mxn, x:n) -> ret:m = 
#   concat(sdot(a[0], x), sgemv(a[1:], x))
#
# sgemv_result = sgemv(a, x)
#
# cblas_sgemv(alpha, sgemv_result, beta, y) -> ret : m =
#   concat(alpha * sgemv_result[0] + beta * y[0],
#          cblas_sgemv(alpha, sgemv_result[1:], beta, y[1:]))
#
# TODO: Support case of incx, incy (look at stored patch)
def getTargetLang():
    a = Var("a", ListT(ListT(Int())))
    alpha = Var("alpha", Int())
    beta = Var("beta", Int())
    x = Var("x", ListT(Int()))
    y = Var("y", ListT(Int()))
    # performs sdot(x, y) - used to implement single output element of MatMul
    sdot = FnDeclRecursive(
       "sdot",
       Int(),
       Ite(
        # condition
        And(And(Gt(Call("list_length", Int(), x), IntLit(0)),
                Gt(Call("list_length", Int(), y), IntLit(0))),
            Eq(Call("list_length", Int(), x), Call("list_length", Int(), y))),
        # if
        Add(
            Mul(Call("list_get", Int(), x, IntLit(0)), Call("list_get", Int(), y, IntLit(0))),
            Call("sdot", Int(), Call("list_tail", ListT(Int()), x, 1), Call("list_tail", ListT(Int()), y, 1))
           ),
        # else
        IntLit(0)
       ),
      x, y
    )
    # performs a : m_n * x : n_1 matrix-vector multiplication to produce m_1 output.
    sgemv = FnDeclRecursive(
       "sgemv",
       ListT(Int()),
       Ite( # condition
            Eq(Call("list_length", Int(), x),
               Call("list_length", Int(), Call("list_get", Int(), a, IntLit(0)))),
            # if
            Call("list_prepend", ListT(Int()),
                  Call("sdot", Int(), Call("list_get", ListT(Int()), a, IntLit(0)), x),
                  Call("sgemv", ListT(Int()), Call("list_tail", ListT(Int()), a, IntLit(1)), x)),
            # else
            Call("list_empty", ListT(Int()))
       ),
       a, x
    )

    def cblas_sgemv_body (alpha, a, x, beta, y):
      # sgemv_result is of size m_1.
      sgemv_result = Call("sgemv", ListT(Int()), a, x)
      return Ite(# condition
              And(
                And(
                  Gt(Call("list_length", Int(), sgemv_result), IntLit(0)),
                  Gt(Call("list_length", Int(), y), IntLit(0))
                ),
                Eq(Call("list_length", Int(), sgemv_result), Call("list_length", Int(), y))
              ),
              # if
              Call("list_prepend", ListT(Int()),
                Add(Mul(alpha, Call("list_get", Int(), sgemv_result, IntLit(0))),
                    Mul(beta,  Call("list_get", Int(), y, IntLit(0)))),
                Call("cblas_sgemv", ListT(Int()),
                      alpha, Call("list_tail", ListT(Int()), sgemv_result, IntLit(1)),
                      beta,  Call("list_tail", ListT(Int()), y, IntLit(1)))),
              # else
              Call("list_empty", ListT(Int()))
      )
      
    # implements semantics of cblas_sgemv: output size will be m_1.
    cblas_sgemv = FnDeclRecursive(
        "cblas_sgemv",
        ListT(Int()),
        cblas_sgemv_body(alpha, a, x, beta, y),
        alpha, a, x, beta, y
    )
    return [cblas_sgemv]