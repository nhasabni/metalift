from analysis import CodeInfo
import subprocess
import pyparsing as pp
import os

from ir import *

import typing
from typing import Any, Callable, Dict, Generator, Union


def flatten(L: typing.List[Any]) -> Generator[str, str, None]:
    for l in L:
        if isinstance(l, list):
            yield "("
            yield from flatten(l)
            yield ")"
        else:
            yield l
    return None


def generateAST(expr: str) -> typing.List[Any]:
    s_expr = pp.nestedExpr(opener="(", closer=")")
    parser = pp.ZeroOrMore(pp.Suppress("(exit)") | s_expr)
    ast = parser.parseString(expr, parseAll=True).asList()
    return ast  # type: ignore


def extractFuns(
    targetLang: typing.List[Expr],
) -> typing.Tuple[typing.List[str], typing.List[Type]]:
    funName, returnType = (
        [],
        [],
    )
    for t in targetLang:
        funName.append(t.args[0])
        returnType.append(t.type)
    return funName, returnType


def generateCandidates(
    invAndPs: typing.List[Expr],
    line: str,
    funName: typing.List[str],
    returnType: typing.List[Type],
) -> typing.Tuple[typing.List[Expr], Dict[str, Expr]]:
    candidates, candidatesExpr = [], {}
    ast = generateAST(line)
    for ce in invAndPs:
        name = ce.args[0]
        for a in ast[0]:
            if name in a:
                candidatesExpr[a[0]] = toExpr(a[1], funName, returnType)
                if isinstance(a[1], str):
                    candidates.append(FnDecl(ce.args[0], ce.type, a[1], *ce.args[2:]))
                else:
                    candidates.append(
                        FnDecl(
                            ce.args[0],
                            ce.type,
                            "(" + " ".join(list(flatten(a[1]))) + ")",
                            *ce.args[2:],
                        )
                    )
    return candidates, candidatesExpr


def toExpr(
    ast: typing.List[Any], funName: typing.List[str], returnType: typing.List[Type]
) -> Expr:
    expr_bi: Dict[str, Callable[..., Expr]] = {
        "=": Eq,
        "+": Add,
        "-": Sub,
        "*": Mul,
        "<": Lt,
        "<=": Le,
        ">": Gt,
        ">=": Ge,
        "and": And,
        "or": Or,
        "=>": Implies,
    }
    expr_uni = {"not": Not}
    if ast[0] in expr_bi.keys():
        return expr_bi[ast[0]](
            toExpr(ast[1], funName, returnType), toExpr(ast[2], funName, returnType)
        )
    elif ast[0] in expr_uni.keys():
        return expr_uni[ast[0]](toExpr(ast[1], funName, returnType))
    elif ast[0] in funName:
        index = funName.index(ast[0])
        return Call(ast[0], returnType[index], toExpr(ast[1], funName, returnType))
    else:
        return ast  # type: ignore


def synthesize(
    basename: str,
    targetLang: typing.List[Expr],
    vars: typing.List[Expr],
    invAndPs: typing.List[Expr],
    preds: Union[str, typing.List[Expr]],
    vc: Expr,
    loopAndPsInfo: typing.List[CodeInfo],
    cvcPath: str,
) -> typing.List[Expr]:
    synthDir = "./synthesisLogs/"
    if not os.path.exists(synthDir):
        os.mkdir(synthDir)
    sygusFile = synthDir + basename + ".sl"

    # Generate sygus file for synthesis
    toSMT(
        targetLang,
        "\n\n".join([str(i) for i in invAndPs]),
        vars,
        preds,
        vc,
        sygusFile,
        True,
    )

    logfile = synthDir + basename + "_logs.txt"
    synthLogs = open(logfile, "w")
    # Run synthesis subprocess
    proc = subprocess.Popen(
        [cvcPath, "--lang=sygus2", "--output=sygus", "--no-sygus-pbe", sygusFile],
        stdout=synthLogs,
    )  # ,stderr=subprocess.DEVNULL)

    funName, returnType = extractFuns(targetLang)
    logfileVerif = open(logfile, "r")
    while True:
        line = logfileVerif.readline()
        if "fail" in line:
            print("SyGuS failed")
            break
        elif "sygus-candidate" in line:
            print("Current PS and INV Guess ", line)
            candidates, _ = generateCandidates(invAndPs, line, funName, returnType)
            candDef = "\n\n".join(str(d) for d in candidates)
            smtFile = synthDir + basename + ".smt"
            toSMT(targetLang, candDef, vars, preds, vc, smtFile, False)

            # run external verification subprocess
            procVerify = subprocess.Popen(
                [cvcPath, "--lang=smt", "--fmf-fun", "--tlimit=10000", smtFile],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
            )
            output = procVerify.stdout.readline()  # type: ignore
            output = output.decode("utf-8")
            if output.count("unsat") == 1:
                print("UNSAT\n")
                print("Verified PS and INV Candidates ", candDef)

                return candidates
            else:
                print("CVC4 verification Result for Current Guess")
                print("SAT\n")
        else:
            continue

    raise Exception("SyGuS failed")


# print to a file
def toSMT(
    targetLang: typing.List[Expr],
    invAndPs: str,
    vars: typing.List[Expr],
    preds: Union[str, typing.List[Expr]],
    vc: Expr,
    outFile: str,
    isSynthesis: bool,
) -> None:
    # order of appearance: inv and ps grammars, vars, non inv and ps preds, vc
    with open(outFile, mode="w") as out:
        out.write("\n\n".join([str(t) for t in targetLang]))

        out.write("\n\n%s\n\n" % invAndPs)

        var_decl_command = "declare-var" if isSynthesis else "declare-const"

        declarations: typing.List[typing.Tuple[str, Type]] = []
        for v in vars:
            genVar(v, declarations)

        declsString = "\n".join(
            ["(%s %s %s)" % (var_decl_command, d[0], d[1]) for d in declarations]
        )  # name and type
        out.write("%s\n\n" % declsString)

        if isinstance(preds, str):
            out.write("%s\n\n" % preds)

        # a list of Exprs - print name, args, return type
        elif isinstance(preds, list):
            preds = "\n".join(
                [
                    "(define-fun %s (%s) (%s) )"
                    % (
                        p.args[0],
                        " ".join("(%s %s)" % (a.args[0], a.type) for a in p.args[1:]),
                        p.type,
                    )
                    for p in preds
                ]
            )
            out.write("%s\n\n" % preds)

        else:
            raise Exception("unknown type passed in for preds: %s" % preds)

        if isSynthesis:
            out.write("%s\n\n" % Constraint(vc))
            out.write("(check-synth)")
        else:
            out.write("%s\n\n" % MLInst_Assert(Not(vc)))
            out.write("(check-sat)\n(get-model)")