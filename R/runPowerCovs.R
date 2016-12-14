# Generated on 2016-06-14 09:41:25
source("/share/PI/manishad/PCORI/lib/prioritySwitch.R")
.sbatch_paths <- c(
# "/share/PI/manishad/PCORI/power/src/sbatch/0001_power.sbatch",
# "/share/PI/manishad/PCORI/power/src/sbatch/0002_power.sbatch",
"/share/PI/manishad/PCORI/power/src/sbatch/0003_power.sbatch",
"/share/PI/manishad/PCORI/power/src/sbatch/0004_power.sbatch",
"/share/PI/manishad/PCORI/power/src/sbatch/0005_power.sbatch",
"/share/PI/manishad/PCORI/power/src/sbatch/0006_power.sbatch",
"/share/PI/manishad/PCORI/power/src/sbatch/0007_power.sbatch",
"") 
prioritySwitch(.sbatch_paths)