# https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("rtCol", "responseCol", "response_num",
                                "response_char", ".", "acc",
                                "accAdjust", "n", "a", "response",
                           "responseSim", "rt", "rt0", "rt0Sim", "rt1",
                           "rt1Sim", "rt1Var", "rtOverall", "rtOverallSim",
                           "rtVar", "t0_Ter", "v"))
}
