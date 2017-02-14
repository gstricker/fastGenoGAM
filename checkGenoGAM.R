library(devtools)
library(BiocCheck)
library(futile.logger)
flog.threshold(ERROR)
document("~/workspace/fastGenoGAM")
test("~/workspace/fastGenoGAM")
check("~/workspace/fastGenoGAM")
BiocCheck("~/workspace/fastGenoGAM")
install("~/workspace/fastGenoGAM")

##TODO put capture = TRUE in logging of complex data structures
