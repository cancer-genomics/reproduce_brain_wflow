### ARTEMIS-DELFI Model

The scripts for training and evaluation of the ARTEMIS-DELFI model.

```
sbatch --array=1-484 A_LOO.sh

sbatch --array=1 B_Lock.sh

Rscript C_scores.r

Rscript D_Model_Importance.r
```
