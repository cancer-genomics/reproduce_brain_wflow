
library(tidyverse)
library(epiR)


source('functions.r')

args = commandArgs(trailingOnly = TRUE)


random.seed = as.numeric(args[1])
out.path = args[2]

#--------------------------------#
set.seed(random.seed)
#--------------------------------#
# learn hyperparameters for the distribution of sens and spec for delfi test
sens.params <- epi.betabuster(0.44, conf=0.85, imsure="greater than", x= 0.40)
spec.params <- epi.betabuster(0.99, conf=0.85, imsure="greater than", x=0.95)

# similarly, learn the hyper prarameters for imaging
imaging.sens <- epi.betabuster(0.985, conf=0.95, imsure="greater than", x=0.95)
imaging.spec <- epi.betabuster(0.985, conf=0.95, imsure="greater than", x=0.95)

# learn hyper parameters for clinician decision making
clinician.sens <- epi.betabuster(0.6, conf=0.95, imsure="greater than", x=0.4)
clinician.spec <- epi.betabuster(0.95, conf=0.95, imsure="greater than", x=0.8)
#--------------------------------#

n_soc = rbinom(n = 1, size = 10e6, prob = 0.125)
n_frag =  10e6 - n_soc


# simulated intended use population given the size and disease prevalence
IU <- intended_use(n_frag, 0.00045)

#  detection status in the intended use population
IU$delfi <- test(IU, sens.params, spec.params)

# performance of delfi test
perf <- performance(IU, "delfi")

# simulate the triage
go.to.imaging <- filter(IU, delfi=="+")
go.home <- filter(IU, delfi=="-")
n.image <- nrow(go.to.imaging)

go.to.imaging$radiology <- test(go.to.imaging, imaging.sens, imaging.spec)
perf2 <- performance(go.to.imaging, "radiology")

cancers.detected = perf2$tp

#------------------------------#
# within the branch that get the delfi test
out = data.frame(n_frag = n_frag, 
                 n_frag.cancer = sum(IU$cancer),
		 n_frag.scans = sum(IU$delfi == '+'),
		 n_frag.detected_cancer = cancers.detected)

#------------------------------#
# within the branch that get direct diagnostic CT

direct <- intended_use(n_soc, 0.0015)
direct$radiology <- test(direct, imaging.sens, imaging.spec)
perf3 = performance(direct, 'radiology')

out = cbind(out, data.frame(n_soc = n_soc,
                            n_soc.cancer = sum(direct$cancer),
			    n_soc.detected_cancer = perf3$tp))
			    
#------------------------------#
# summing stats for the combined approach

out$total.scans = out$n_soc + out$n_frag.scans
out$total.detected_cancer = out$n_soc.detected_cancer + out$n_frag.detected_cancer


#------------------------------#
write.table(out, file = out.path, quote = FALSE, row.names = FALSE, sep = '\t')




