# SCCM


## Example data

1. First load some libraries.

```R
library(tidyverse)
library(openxlsx)
library(sccm)
```

2. Next download the example data by Mayhugh et al. (2018) from [DRYAD](https://datadryad.org/stash/dataset/doi:10.5061/dryad.p63d200) and extract the zip file.

3. Read the 29th participant from data.

```R
dat <- read.xlsx("Level1_RandomUpSleep_1365_PlosOne_Mayhugh_etal.xlsx") %>%
    as_tibble() %>%
    filter(Anom_ID == 29) %>%
    select(d_abstained, Stress, Crave)
```

## Example using code

```R
mod <- sccm(
    X = "d_abstained",
    M = "Stress",
    Y = "Crave",
    dat = dat,
    lag = 1,
    permutation = T,
    perm_reps = 1000
)

# get a summary
summary(mod)

# plot a path diagram (DAG)
plot(mod)

# print lavaan summary
print(mod)

# print lavaan model syntax
cat(mod$syntax)
```

## Example using the GUI

1. Write the data to a .csv file.

```R
write.csv(dat, "sccm.csv", row.names = F)
```

2. Start the GUI.

```R
sccm_gui()
```

## References

Mayhugh, R. E., Rejeski, W. J., Petrie, M. R., Laurienti, P. J., & Gauvin, L. (2018). Differing patterns of stress and craving across the day in moderate-heavy alcohol consumers during their typical drinking routine and an imposed period of alcohol abstinence (Z. C. Lin, Ed.). *PLOS ONE, 13*(4), e0195063. https://doi.org/10.1371/journal.pone.0195063
