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

    int_read_vars = filter_types(ci.readVars, Int())
    list_read_vars = filter_types(ci.readVars, ListT(Int()))
    matrix_vars = filter_types(ci.readVars, ListT(ListT(Int())))
    int_modified_vars = filter_types(ci.modifiedVars, Int())
    list_modified_vars = filter_types(ci.modifiedVars, ListT(Int()))

    alpha_beta = Choose(*int_read_vars)
    alpha = int_read_vars[0]
    beta = int_read_vars[1]
    #x_y = Choose(*list_read_vars)
    a = matrix_vars[0]
    x = list_read_vars[0]
    y = list_read_vars[1]

    if name.startswith("inv0"):  #inv for outer loop
      intChoices = Choose(IntLit(0), IntLit(1))
      z = list_modified_vars[0]
      #i_j_ret = Choose(*int_modified_vars)
      #print(*int_modified_vars)
      #i = int_modified_vars[0]
      #i = Choose(*int_modified_vars)
      i = int_modified_vars[2] # i13
      temp = int_modified_vars[0] # i14
      j = int_modified_vars[1]
      
      # i >= 0
      i_lower_bound = Choose(Ge(i, intChoices), Gt(i, intChoices),
                             Le(i, intChoices), Lt(i, intChoices),
                             Eq(i, intChoices))
      # i <= listlen
      i_upper_bound = Choose(Ge(i, Call("list_length", Int(), y)),
                             Gt(i, Call("list_length", Int(), y)),
                             Le(i, Call("list_length", Int(), y)),
                             Lt(i, Call("list_length", Int(), y)),
                             Eq(i, Call("list_length", Int(), y)))
      # single iteration of outer loop computes one element of output vector by performing:
      # z[i] = alpha * sdot(a[i], x) + beta * y[i]
      outer_loop_output = Choose(Eq(z, 
                                    Call("cblas_sgemv", ListT(Int()), alpha_beta,
                                         Call("list_list_take", ListT(Int()), a, i),
                                         x, alpha_beta, Call("list_take", ListT(Int()), y, i))))
                                        #Add(
                                        #    Mul(alpha_beta,
                                        #        Call("sdot", Int(), Call("list_list_get", ListT(Int()), a, i), x)),
                                        #    Mul(alpha_beta, Call("list_get", Int(), y, i))
                                        #)))) 
      inv = Choose(Implies(And(Eq(Call("list_length", Int(), x),
                                   Call("list_length", Int(),
                                        Call("list_list_get", ListT(Int()), a, IntLit(0)))),
                               Eq(Call("list_length", Int(), y), Call("list_list_length", Int(), a)),
                               Gt(Call("list_length", Int(), a), IntLit(1))
                           ),
                           And(outer_loop_output, And(i_lower_bound, i_upper_bound))
                           ))
      return inv
    
    elif name.startswith("inv1"):  #inv for inner loop
      intChoices = Choose(IntLit(0), IntLit(1))
      #i_j_ret = Choose(*int_modified_vars)
      z = list_modified_vars[0]
      i = int_modified_vars[2] # i13
      temp = int_modified_vars[0] # i14
      j = int_modified_vars[1] # i15
    
      # i >= 0
      i_lower_bound = Choose(Ge(i, intChoices), Gt(i, intChoices),
                             Le(i, intChoices), Lt(i, intChoices),
                             Eq(i, intChoices))
      # i <= listlen
      i_upper_bound = Choose(Ge(i, Call("list_length", Int(), a)),
                             Gt(i, Call("list_length", Int(), a)),
                             Le(i, Call("list_length", Int(), a)),
                             Lt(i, Call("list_length", Int(), a)),
                             Eq(i, Call("list_length", Int(), a)))
      # j >= 0
      j_lower_bound = Choose(Ge(j, intChoices), Gt(j, intChoices),
                             Le(j, intChoices), Lt(j, intChoices),
                             Eq(j, intChoices))
      # j <= listlen
      j_upper_bound = Choose(Ge(j, Call("list_length", Int(), x)),
                             Gt(j, Call("list_length", Int(), x)),
                             Le(j, Call("list_length", Int(), x)),
                             Lt(j, Call("list_length", Int(), x)),
                             Eq(j, Call("list_length", Int(), x)))
      # One iteration of inner loop is performing SDOT: for some i, j=0 to n: ret += a[i][j] * x[j]
      #inner_loop_output = Choose(Eq(temp, Add(temp,
      #                                        Mul(Call("list_get", Int(),
      #                                               Call("list_list_get", ListT(Int()), a, i), j),
      #                                            Call("list_get", Int(), x, j)))))
      inner_loop_output = Choose(Eq(temp,
                                    Call("sdot", Int(),
                                         Call("list_take", ListT(Int()),
                                              Call("list_list_get", ListT(Int()), a, i), j),
                                         Call("list_take", ListT(Int()), x, j)
                                    ))
                                 )
      # and using SDOT output to perform cblas_sgemv:
      # contents of z in inner_loop = alpha * sdot(a[i-1] * x) + beta * y[i-1]
      #outer_loop_output = Choose(Eq(z,
      #                              Call("list_append", ListT(Int()), z, 
      #                                  Add(
      #                                      Mul(alpha,
      #                                         Call("sdot", Int(), Call("list_list_get", ListT(Int()), a, i), x)),
      #                                      Mul(beta, Call("list_get", Int(), y, i))
      #                                  ))))
      outer_loop_output = Choose(Eq(z, 
                                    Call("cblas_sgemv", ListT(Int()), alpha,
                                         Call("list_list_take", ListT(Int()), a, i), x, beta, Call("list_take", ListT(Int()), y, i))))
      inv = Choose(Implies(And(
                                Eq(Call("list_length", Int(), x),
                                   Call("list_length", Int(), Call("list_list_get", ListT(Int()), a, IntLit(0)))),
                                Eq(Call("list_length", Int(), y), Call("list_list_length", Int(), a)),
                                Gt(Call("list_length", Int(), a), IntLit(1))
                              ),
                           And(
                              And(inner_loop_output,
                                  And(
                                    i_lower_bound, j_lower_bound,
                                    i_upper_bound, j_upper_bound
                                  ),
                              ),
                              outer_loop_output)))
      return inv

    else:  #ps
      z = Choose(*list_modified_vars)
      choices = Choose(Implies(And(
                                    Eq(Call("list_length", Int(), x),
                                       Call("list_length", Int(), Call("list_list_get", ListT(Int()), a, IntLit(0)))),
                                    Eq(Call("list_length", Int(), y), Call("list_list_length", Int(), a)),
                                    Gt(Call("list_length", Int(), a), IntLit(1))
                                  ),
                               Eq(z, Call("cblas_sgemv", ListT(Int()),
                                          alpha_beta, a, x, alpha_beta, y))))
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
            Call("sdot", Int(), Call("list_tail", ListT(Int()), x, IntLit(1)), Call("list_tail", ListT(Int()), y, IntLit(1)))
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
               Call("list_length", Int(), Call("list_list_get", Int(), a, IntLit(0)))),
            # if
            Call("list_prepend", ListT(Int()),
                  Call("sdot", Int(), Call("list_list_get", ListT(Int()), a, IntLit(0)), x),
                  Call("sgemv", ListT(Int()), Call("list_list_tail", ListT(Int()), a, IntLit(1)), x)),
            # else
            Call("list_empty", ListT(Int()))
       ),
       a, x
    )
      
    # implements semantics of cblas_sgemv: output size will be m_1.
    cblas_sgemv = FnDeclRecursive(
        "cblas_sgemv",
        ListT(Int()),
        Ite(# condition
              And(
               And(
                  Gt(Call("list_length", Int(), a), IntLit(0)),
                  Gt(Call("list_length", Int(), y), IntLit(0))
               ),
               Eq(Call("list_length", Int(), a), Call("list_length", Int(), y))
              ),
              # if
              Call("list_prepend", ListT(Int()),
                Add(Mul(alpha, Call("sdot", Int(), Call("list_list_get", ListT(Int()), a, IntLit(0)), x)),
                    Mul(beta,  Call("list_get", Int(), y, IntLit(0)))),
                Call("cblas_sgemv", ListT(Int()),
                      alpha, Call("list_list_tail", ListT(Int()), a, IntLit(1)), x,
                      beta,  Call("list_tail", ListT(Int()), y, IntLit(1)))),
              # else
              Call("list_empty", ListT(Int()))
        ),
        alpha, a, x, beta, y
    )
    return [sdot, sgemv, cblas_sgemv]