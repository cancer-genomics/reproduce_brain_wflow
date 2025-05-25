intended_use <- function(N, theta) tibble(cancer=rbinom(N, 1, theta))

test <- function(IU, se, sp){
    cancer <- IU$cancer
    calls <- rep(NA, nrow(IU))
    Sp <- rbeta(1, sp$shape1, sp$shape2)
    Se <- rbeta(1, se$shape1, se$shape2)
    freq <- table(cancer)
    cancer.calls <- rbinom(freq[2], 1, Se)
    is.cancer <- cancer==1
    calls[ is.cancer ] <- ifelse(cancer.calls==1, "+", "-")
    noncancer.calls <- rbinom(freq[1], 1, Sp)
    calls[ !is.cancer ] <- ifelse(noncancer.calls==1, "-", "+")
    calls
}

performance <- function(IU, testname){
    perf <- IU %>%
        rename(test_={{testname}}) %>%
        summarize(fn=sum(cancer==1 & test_=="-"),
                  tn=sum(cancer==0 & test_=="-"),
                  fp=sum(cancer==0 & test_=="+"),
                  tp=sum(cancer==1 & test_=="+"),
                  tpr=tp/(tp+fn), 
                  fpr=fp/(fp+tn),
                  fnr=fn/(fn+tp),
                  sens=tpr,
                  spec=1-fpr,
                  ## Since we are simulating a prospective study (not a case-control),
                  ## we can estimate these directly as well
                  ppv=tp/(tp+fp),
                  npv=tn/(tn+fn))
    perf
}
