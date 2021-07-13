


# psidR: make building panel data from PSID easy


[![Build Status](https://travis-ci.org/floswald/psidR.svg)](https://travis-ci.org/floswald/psidR)

[![Rdoc](http://www.rdocumentation.org/badges/version/psidR)](http://www.rdocumentation.org/packages/psidR)


This R package provides a function to easily build panel data from PSID raw data.

>**Warning**: the wealth-supplement setup has changed on the PSID system. wealth variables are now part of the family files for waves 1999 onwards. The `wealth=TRUE` option has therefore been removed from the package. See [this issue](https://github.com/floswald/psidR/issues/34) for more details.
>

## How to install this package

The package is on CRAN, so just type

```r
install.packages('psidR')
```

Alternatively to get the up-to-date version from this repository,

```r
install.packages('devtools')
install_github("psidR",username="floswald")
```

### PSID

The [Panel Study of Income Dynamics](http://psidonline.isr.umich.edu/) is a publicly available dataset. 

* you can use the [data center](http://simba.isr.umich.edu/default.aspx) to build simple datasets
* not workable for larger datasets
  * some variables don't show up (although you know they exist)
  * the ftp interface gets slower the more periods you are looking at
  * the click and scroll exercise of selecting the right variables in each period is extremely error prone. 
* merging the data manually is non-trivial.

### psidR

This package attempts to help the task of building a panel dataset. The user directly downloads ASCII data from the PSID server into `R`, **without the need** for any other software like stata or sas. To build the panel, the user must then specify the variable names in each wave of the questionnaire in a data.frame `fam.vars`, as well as the variables from the individual index in `ind.vars`. The helper function `getNamesPSID` is helpful in finding different variable names across waves - see examples below.


### Quick Start

```R
> library(psidR)

> build.psid(datadr = "~/data/PSID", small = TRUE)  # directory `datadr` must exist!
INFO [2021-07-13 10:34:26] Will download missing datasets now
INFO [2021-07-13 10:34:26] will download family files: 2013, 2015
INFO [2021-07-13 10:34:26] will download latest individual index: IND2019ER
This can take several hours/days to download.
 want to go ahead? give me 'yes' or 'no'.yes
please enter your PSID username: *****
please enter your PSID password: *****
INFO [2021-07-13 10:34:41] downloading file ~/data/PSID/FAM2013ER
INFO [2021-07-13 10:34:56] now reading and processing SAS file ~/data/PSID/FAM2013ER into R
INFO [2021-07-13 10:40:06] downloading file ~/data/PSID/FAM2015ER          
INFO [2021-07-13 10:40:22] now reading and processing SAS file ~/data/PSID/FAM2015ER into R
INFO [2021-07-13 10:45:34] downloading file ~/data/PSID/IND2019ER          
INFO [2021-07-13 10:46:39] now reading and processing SAS file ~/data/PSID/IND2019ER into R
INFO [2021-07-13 11:15:04] finished downloading files to ~/data/PSID/       
INFO [2021-07-13 11:15:04] continuing now to build the dataset
INFO [2021-07-13 11:15:04] psidR: Loading Family data from .rda files
INFO [2021-07-13 11:15:12] psidR: loaded individual file: ~/data/PSID/IND2019ER.rda
INFO [2021-07-13 11:15:12] psidR: total memory load in MB: 1538
INFO [2021-07-13 11:15:12] psidR: currently working on data for year 2013
INFO [2021-07-13 11:15:12] full 2013 sample has 82573 obs
INFO [2021-07-13 11:15:12] you selected 34856 obs belonging to SRC
INFO [2021-07-13 11:15:12] dropping non-heads leaves 5450 obs
INFO [2021-07-13 11:15:14] psidR: currently working on data for year 2015
INFO [2021-07-13 11:15:14] full 2015 sample has 82573 obs
INFO [2021-07-13 11:15:14] you selected 34856 obs belonging to SRC
INFO [2021-07-13 11:15:14] dropping non-heads leaves 5318 obs
INFO [2021-07-13 11:15:16] End of build.panel
```

### Usage

First present a real world example building a full 1968-2017 panel. Then we show some tests.


### Real World Example: With Missing Variables

* You want a `data.table` with the following columns: `PID,year,income,wage,age,educ` and some more variables.
* You went to the [PSID variable search](https://simba.isr.umich.edu/VS/s.aspx) to look up the relevant variable names in each year in either the `individual-level` or `family-level` datasets.
* You created a list of those variables as I did in [`inst/psid-lists`](inst/psid-lists) of this package
* You noted that there is **NO EDUCATION** variable in the individual index file in 1968 and 1969
    * Instead of the variable name for `EDUC` in 1968 and 1969 you want to put `NA`
* You noted that there is **NO HOURLY WAGE** variable in the family index file in 1993
    * Instead of the variable name for `HOURLY WAGE` in 1993 you want to put `NA`

```R
# Build panel with income, wage, age, education and several other variables
# [this is the body of the function build.psid()]
library(psidR)
library(data.table)
r = system.file(package="psidR")
f = fread(file.path(r,"psid-lists","famvars.txt"))
i = fread(file.path(r,"psid-lists","indvars.txt"))

> i
                           dataset year variable                  label   name
  1: PSID Individual Data by Years 1968  ER30019   INDIVIDUAL WEIGHT 68 weight
  2: PSID Individual Data by Years 1969  ER30042   INDIVIDUAL WEIGHT 69 weight
  3: PSID Individual Data by Years 1970  ER30066   INDIVIDUAL WEIGHT 70 weight
  4: PSID Individual Data by Years 1971  ER30090   INDIVIDUAL WEIGHT 71 weight
  5: PSID Individual Data by Years 1972  ER30116   INDIVIDUAL WEIGHT 72 weight
 ---                                                                          
143:    PSID Individual Data Index 2009  ER34020 HIGHEST GRADE FINISHED   educ
144:    PSID Individual Data Index 2011  ER34119 HIGHEST GRADE FINISHED   educ
145:    PSID Individual Data Index 2013  ER34230 HIGHEST GRADE FINISHED   educ
146:    PSID Individual Data Index 2015  ER34349 HIGHEST GRADE FINISHED   educ
147:    PSID Individual Data Index 2017  ER34548 HIGHEST GRADE FINISHED   educ

> f
                   dataset year variable                     label            name
  1: PSID Main Family Data 1968      V47 HD ANN HRS WORKED LAST YR           hours
  2: PSID Main Family Data 1969     V465 HD ANN HRS WORKED LAST YR           hours
  3: PSID Main Family Data 1970    V1138 HD ANN HRS WORKED LAST YR           hours
  4: PSID Main Family Data 1971    V1839 HD ANN HRS WORKED LAST YR           hours
  5: PSID Main Family Data 1972    V2439 HD ANN HRS WORKED LAST YR           hours
 ---                                                                              
609:     PSID Family-level 2009  ER42139  A52 LIKELIHOOD OF MOVING likelihood_move
610:     PSID Family-level 2011  ER47447  A52 LIKELIHOOD OF MOVING likelihood_move
611:     PSID Family-level 2013  ER53147  A52 LIKELIHOOD OF MOVING likelihood_move
612:     PSID Family-level 2015  ER60162  A52 LIKELIHOOD OF MOVING likelihood_move
613:     PSID Family-level 2017  ER66163  A52 LIKELIHOOD OF MOVING likelihood_move

# alternatively, use `getNamesPSID`:
# cwf <- read.xlsx("http://psidonline.isr.umich.edu/help/xyr/psid.xlsx")
# Suppose you know the name of the variable in a certain year, and it is
# "ER17013". then get the correpsonding name in another year with
# getNamesPSID("ER17013", cwf, years = 2001)  # 2001 only
# getNamesPSID("ER17013", cwf, years = 2003)  # 2003
# getNamesPSID("ER17013", cwf, years = NULL)  # all years
# getNamesPSID("ER17013", cwf, years = c(2005, 2007, 2009))   # some years

# next, bring into required shape:

i = dcast(i[,list(year,name,variable)],year~name, value.var = "variable")
f = dcast(f[,list(year,name,variable)],year~name, value.var = "variable")

> head(i)
   year     age    educ empstat  weight
1: 1968 ER30004 ER30010    <NA> ER30019
2: 1969 ER30023    <NA>    <NA> ER30042        # NOTICE THE NA for educ HERE!!
3: 1970 ER30046 ER30052    <NA> ER30066
4: 1971 ER30070 ER30076    <NA> ER30090
5: 1972 ER30094 ER30100    <NA> ER30116
6: 1973 ER30120 ER30126    <NA> ER30137

> head(f)
   year age_youngest_child debt empstat_ faminc hours hvalue ...
1: 1968               V120 <NA>     V196    V81   V47     V5 ...
2: 1969              V1013 <NA>     V639   V529  V465   V449 ...
3: 1970              V1243 <NA>    V1278  V1514 V1138  V1122 ...
4: 1971              V1946 <NA>    V1983  V2226 V1839  V1823 ...
5: 1972              V2546 <NA>    V2581  V2852 V2439  V2423 ...
6: 1973              V3099 <NA>    V3114  V3256 V3027  V3021 ...

# call the builder function

d = build.panel(datadir=datadr,fam.vars=f,ind.vars=i, heads.only = TRUE,sample="SRC",design="all")

# d contains your panel

save(d,file="~/psid.Rds")
```

Here are some tests:

```R
# one year test, no ind file
# call function `small.test.noind()`
# get var names from cross walk
cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003))
# create family vars data.frame
famvars = data.frame(year=c(2003),age=head_age_var_name)
# call function
build.panel(fam.vars=famvars,datadir=dd)

# one year test, ind file
# call function `small.test.ind()`

cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003))
educ = getNamesPSID("ER30323",cwf,years=2003)
famvars = data.frame(year=c(2003),age=head_age_var_name)
indvars = data.frame(year=c(2003),educ=educ)
build.panel(fam.vars=famvars,ind.vars=indvars,datadir=dd)


# three year test, ind file
# call function `medium.test.ind()`

cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003,2005,2007))
educ = getNamesPSID("ER30323",cwf,years=c(2003,2005,2007))
famvars = data.frame(year=c(2003,2005,2007),age=head_age_var_name)
indvars = data.frame(year=c(2003,2005,2007),educ=educ)
build.panel(fam.vars=famvars,ind.vars=indvars,datadir=dd)

# etc for
medium.test.noind()

# example output:

INFO [2018-10-10 10:58:23] Will download missing datasets now
INFO [2018-10-10 10:58:23] will download family files: 2003, 2005, 2007
This can take several hours/days to download.
 want to go ahead? give me 'yes' or 'no'.yes
please enter your PSID username: *******
please enter your PSID password: *******
INFO [2018-10-10 10:58:46] downloading file ~/psid/FAM2003ER
INFO [2018-10-10 10:58:50] now reading and processing SAS file ~/psid/FAM2003ER into R
INFO [2018-10-10 11:07:02] downloading file ~/psid/FAM2005ER               
INFO [2018-10-10 11:07:05] now reading and processing SAS file ~/psid/FAM2005ER into R
INFO [2018-10-10 11:14:44] downloading file ~/psid/FAM2007ER               
INFO [2018-10-10 11:14:48] now reading and processing SAS file ~/psid/FAM2007ER into R
INFO [2018-10-10 11:28:25] finished downloading files to ~/psid/           
INFO [2018-10-10 11:28:25] continuing now to build the dataset
INFO [2018-10-10 11:28:25] psidR: Loading Family data from .rda files
INFO [2018-10-10 11:28:34] psidR: loaded individual file: ~/psid/IND2015ER.rda
INFO [2018-10-10 11:28:34] psidR: total memory load in MB: 1252
INFO [2018-10-10 11:28:34] 
INFO [2018-10-10 11:28:34] psidR: currently working on data for year 2003
INFO [2018-10-10 11:28:36] 
INFO [2018-10-10 11:28:36] psidR: currently working on data for year 2005
INFO [2018-10-10 11:28:37] 
INFO [2018-10-10 11:28:37] psidR: currently working on data for year 2007
INFO [2018-10-10 11:28:39] balanced design reduces sample from 97377 to 89571
INFO [2018-10-10 11:28:39] End of build.panel
> x
       age interview ID1968 pernum sequence relation.head     pid year
    1:  92         1    848      2        1            10  848002 2003
    2:  64         2   1173      1        1            10 1173001 2003
    3:  48         3   1866     32        2            30 1866032 2003
    4:  48         3   1866    171        1            10 1866171 2003
    5:  48         3   1866    175        0             0 1866175 2003
   ---                                                                
89567:  49      8332   6069      4        2            20 6069004 2007
89568:  49      8332   6069     30        0             0 6069030 2007
89569:  49      8332   6069    171        3            33 6069171 2007
89570:  49      8332   6069    173        1            10 6069173 2007
89571:  49      8332   6069    174        0             0 6069174 2007

# etc for 
medium.test.ind.NA()
```








### Example Usage

the main function in the package is `build.panel` and it has a reproducible example which you can look at by typing

```r
require(psidR)
example(build.panel)
```

### Supplemental Datasets

The PSID has a wealth of add-on datasets. Once you have a panel those are easy to merge on. The panel will have a variable `interview`, which is the identifier in the supplemental dataset. 


## Citation

If you use `psidR` in your work, please consider citing it. You could just do 

```R
> citation(package="psidR")

To cite the 'psidR' package in publications use:

  Florian Oswald (2021). psidR: Build Panel Data Sets from PSID Raw Data. R package version
  2.1.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {psidR: Build Panel Data Sets from PSID Raw Data},
    author = {Florian Oswald},
    year = {2021},
    note = {R package version 2.1},
    url = {https://github.com/floswald/psidR},
  }
```

Thanks!
