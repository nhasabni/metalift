import sys

from metalift.analysis import analyze, CodeInfo
#from metalift.analysis_new import analyze, VariableTracker
from metalift.synthesize_auto import synthesize
from metalift.ir import *

# Import grammar and target language for specific APIs
import api.mylistapis as mylistapis
import api.bzero as bzero
import api.cblas_saxpy as cblas_saxpy
import api.cblas_sdot as cblas_sdot
import api.cblas_sgemv as cblas_sgemv

def grammar(ci: CodeInfo):
    #g = Or(mylistapis.getGrammar(ci),
    #        bzero.getGrammar(ci))
    #g = cblas_saxpy.getGrammar(ci)
    g = cblas_sgemv.getGrammar(ci)
    #g = cblas_sdot.getGrammar(ci)
    #g = bzero.getGrammar(ci)
    #g = mylistapis.getGrammar(ci)
    return Synth(ci.name, g, *ci.modifiedVars, *ci.readVars)


def targetLang():
    return cblas_sgemv.getTargetLang()
    #return cblas_sdot.getTargetLang()
    #return cblas_saxpy.getTargetLang()
    #return bzero.getTargetLang()
    #return mylistapis.getTargetLang() + bzero.getTargetLang()
    #return mylistapis.getTargetLang()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Script require 2 inputs:<basename><function_name>")
        sys.exit(0)
    
    basename = sys.argv[1]
    fnName = sys.argv[2]

    filename = basename + ".ll"
    loopsFile = basename + ".loops"
    cvcPath = "cvc5"

    (vars, invAndPs, preds, vc, loopAndPsInfo) = analyze(filename, fnName, loopsFile)

    print("====== synthesis")
    invAndPs = [grammar(ci) for ci in loopAndPsInfo]

    lang = targetLang()

    candidates = synthesize(
        basename,
        lang,
        vars,
        invAndPs,
        preds,
        vc,
        loopAndPsInfo,
        cvcPath,
        noVerify=True
    )

    # test_analysis = analyze(filename, fnName, loopsFile)

    # variable_tracker = VariableTracker()
    # base = variable_tracker.variable("base", ListT(Int()))

    # synth_fun = grammar(fnName, [base], Var("ret", Int()))

    # vc = test_analysis.call(base)(variable_tracker, lambda ret: Call(
    #     fnName,
    #     Bool(),
    #     ret,
    #     base
    # ))

    # vars = variable_tracker.all()
    # loopAndPsInfo = [synth_fun]

    # print("====== grammar")
    # invAndPs = [synth_fun]

    # lang = targetLang()

    # # rosette synthesizer  + CVC verfication
    # print("====== synthesis")
    # candidates = synthesize(
    #     basename, lang, vars, invAndPs, [], vc, loopAndPsInfo, cvcPath
    # )

    print("====== verified candidates")
    for c in candidates:
        print(c, "\n")
