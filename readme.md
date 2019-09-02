


# psidR: make building panel data from PSID easy


[![Build Status](https://travis-ci.org/floswald/psidR.svg)](https://travis-ci.org/floswald/psidR)

[![Rdoc](http://www.rdocumentation.org/badges/version/psidR)](http://www.rdocumentation.org/packages/psidR)


This R package provides a function to easily build panel data from PSID raw data.

### PSID

The [Panel Study of Income Dynamics](http://psidonline.isr.umich.edu/) is a publicly available dataset. You have to register and agree to terms and conditions, but there are no other strings attached. 

* you can use the [data center](http://simba.isr.umich.edu/default.aspx) to build simple datasets
* not workable for larger datasets
  * some variables don't show up (although you know they exist)
  * the ftp interface gets slower the more periods you are looking at
  * the click and scroll exercise of selecting the right variables in each period is extremely error prone. 
* merging the data manually is tricky.

### psidR

This package attempts to help the task of building a panel data. The user directly downloads ASCII data from the server to and into `R`, without the need for any other software like stata or sas. To build the panel, the user must then specify the variable names in each wave of the questionnaire in a data.frame `fam.vars`, as well as the variables from the individual index in `ind.vars`. The helper function `getNamesPSID` is helpful in finding different variable names across waves - see examples below.

### Usage

check out those example calls and output of function `medium.test.noind()`:

```R
# one year test, no ind file
# call function `small.test.noind()`
# get var names from cross walk
cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003))
# create family vars data.frame
famvars = data.frame(year=c(2003),age=head_age_var_name)
# call function
build.panel(fam.vars=famvars,datadir=dd)

# one year test, ind file
# call function `small.test.ind()`

cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003))
educ = getNamesPSID("ER30323",cwf,years=2003)
famvars = data.frame(year=c(2003),age=head_age_var_name)
indvars = data.frame(year=c(2003),educ=educ)
build.panel(fam.vars=famvars,ind.vars=indvars,datadir=dd)


# three year test, ind file
# call function `medium.test.ind()`

cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003,2005,2007))
educ = getNamesPSID("ER30323",cwf,years=c(2003,2005,2007))
famvars = data.frame(year=c(2003,2005,2007),age=head_age_var_name)
indvars = data.frame(year=c(2003,2005,2007),educ=educ)
build.panel(fam.vars=famvars,ind.vars=indvars,datadir=dd)

# etc for
medium.test.noind()

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
medium.test.ind.NA.wealth()

```



### Real World Example: Missing Variables

* You want a `data.table` with the following columns: `PID,year,income,wage,age,educ`.
* You went to the [PSID variable search](https://simba.isr.umich.edu/VS/s.aspx) to look up the relevant variable names in each year in either the `individual-level` or `family-level` datasets.
* You created a list of those variables as I did in `inst/psid-lists` of this package
* You noted that there is **NO EDUCATION** variable in the individual index file in 1968 and 1969
    * Instead of the variable name for `EDUC` in 1968 and 1969 you want to put `NA`
* You noted that there is **NO HOURLY WAGE** variable in the family index file in 1993
    * Instead of the variable name for `HOURLY WAGE` in 1993 you want to put `NA`

```R
# Build panel with income, wage, age and education
# this is the body of the function build.psid()
library(psidR)
r = system.file(package="psidR")
f = data.table::fread(file.path(r,"psid-lists","famvars.txt"))
i = data.table::fread(file.path(r,"psid-lists","indvars.txt"))

# alternatively, use getNamesPSID:
# cwf <- read.xlsx("http://psidonline.isr.umich.edu/help/xyr/psid.xlsx")
# Suppose you know the name of the variable in a certain year, and it is
# "ER17013". then get the correpsonding name in another year with
# getNamesPSID("ER17013", cwf, years = 2001)  # 2001 only
# getNamesPSID("ER17013", cwf, years = 2003)  # 2003
# getNamesPSID("ER17013", cwf, years = NULL)  # all years
# getNamesPSID("ER17013", cwf, years = c(2005, 2007, 2009))   # some years

# add a group identifier
f[1:38,vgroup := "wage"]
f[39:76,vgroup := "earnings"]
setkey(f,vgroup)

i[1:38,   vgroup := "age"]
i[39:76,  vgroup := "educ"]  # caution about 2 first years: no educ data
i[77:114, vgroup := "weight"]
setkey(i,vgroup)

> head(f)
                 dataset year variable                label   vgroup
1: PSID Main Family Data 1968      V81        FAM MONEY INC earnings
2: PSID Main Family Data 1969     V529       TOTAL FU $ INC earnings
3: PSID Main Family Data 1970    V1514 TOT FU MON INC OV414 earnings
4: PSID Main Family Data 1971    V2226       TOT FU MON INC earnings
5: PSID Main Family Data 1972    V2852       TOT FU MON INC earnings
6: PSID Main Family Data 1973    V3256       TOT FU MON INC earnings
> head(i)
                         dataset year variable                label vgroup
1: PSID Individual Data by Years 1968  ER30004 AGE OF INDIVIDUAL 68    age
2: PSID Individual Data by Years 1969  ER30023 AGE OF INDIVIDUAL 69    age
3: PSID Individual Data by Years 1970  ER30046 AGE OF INDIVIDUAL 70    age
4: PSID Individual Data by Years 1971  ER30070 AGE OF INDIVIDUAL 71    age
5: PSID Individual Data by Years 1972  ER30094 AGE OF INDIVIDUAL 72    age
6: PSID Individual Data by Years 1973  ER30120 AGE OF INDIVIDUAL 73    age


# create ind and fam data.tables
ind = cbind(i[J("age"),list(year,age=variable)],
            i[J("educ"),list(educ=variable)],
            i[J("weight"),list(weight=variable)])
fam = cbind(f[J("wage"),list(year,wage=variable)],
            f[J("earnings"),list(earnings=variable)])

> head(ind)
   year     age    educ  weight
1: 1968 ER30004      NA ER30019
2: 1969 ER30023      NA ER30042
3: 1970 ER30046 ER30052 ER30066
4: 1971 ER30070 ER30076 ER30090
5: 1972 ER30094 ER30100 ER30116
6: 1973 ER30120 ER30126 ER30137
> head(fam)
   year  wage earnings
1: 1968  V337      V81
2: 1969  V871     V529
3: 1970 V1567    V1514
4: 1971 V2279    V2226
5: 1972 V2906    V2852
6: 1973 V3275    V3256

# caution: this step will take many hours the first time.
d = build.panel(datadir="~/data",fam.vars=fam,
          ind.vars=ind,
          heads.only = TRUE,
          sample="SRC",
          design=2)
```



### How to install this package

The package is on CRAN, so just type

```r
install.packages('psidR')
```

Alternatively to get the up-t-date version from this repository,

```r
install.packages('devtools')
install_github("psidR",username="floswald")
```


### Example Usage

the main function in the package is `build.panel` and it has a reproducible example which you can look at by typing

```r
require(psidR)
example(build.panel)
```

### Supplemental Datasets

The PSID has a wealth of add-on datasets. Once you have a panel those are easy to merge on. The panel will have a variable `interview`, which is the identifier in the supplemental dataset:

```R
medium.test.ind.NA.wealth <- function(dd=NULL){
    cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
    head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2005,2007))
    educ = getNamesPSID("ER30323",cwf,years=c(2005,2007))
    educ[2] = NA
    r = system.file(package="psidR")
    w = fread(file.path(r,"psid-lists","wealthvars-small.txt"))
    famvars = data.frame(year=c(2005,2007),age=head_age_var_name)
    indvars = data.frame(year=c(2005,2007),educ=educ)
    build.panel(fam.vars=famvars,ind.vars=indvars,wealth.vars=w,datadir=dd,loglevel = DEBUG)
}
```

## Citation

If you use `psidR` in your work, please consider citing it. You could just do 

```R
> citation(package="psidR")

To cite package ‘psidR’ in publications use:

  Florian Oswald (2018). psidR: Build Panel Data Sets from PSID Raw
  Data. R package version 1.7. https://github.com/floswald/psidR

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {psidR: Build Panel Data Sets from PSID Raw Data},
    author = {Florian Oswald},
    year = {2018},
    note = {R package version 1.7},
    url = {https://github.com/floswald/psidR},
  }
```

Thanks!
