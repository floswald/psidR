


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

this package attempts to help the task of building a panel data. the user can either

1. download ASCII data from the server to disk and process with Stata or SAS to generate .dta or .csv files as input; or
2. `[RECOMMENDED]` use the option to directly download into an R data.frame via the `SAScii` package. You download only once.

To build the panel, the user must specify the variable names in each wave of the questionnaire in a data.frame `fam.vars`, as well as the variables from the individual index in `ind.vars`. The helper function `getNamesPSID` is helpful in finding different variable names across waves. 


### Real World Example: Missing Variables

* You want a `data.table` with the following columns: `PID,year,income,wage,age,educ`.
* You went to the [PSID variable search](https://simba.isr.umich.edu/VS/s.aspx) to look up the relevant variable names in each year in either the `individual-level` or `family-level` datasets.
* You created a list of those variables as I did in `inst/psid-lists` of this package
* You noted that there is **NO EDUCATION** variable in the individual index file in 1968 and 1969
    * Instead of the variable name for `EDUC` in 1968 and 1969 you put `NA`
* You noted that there is **NO HOURLY WAGE** variable in the family index file in 1993
    * Instead of the variable name for `HOURLY WAGE` in 1993 you put `NA`

```R
# Build panel with income, wage, age and education
# this is the body of the function build.psid()
library(psidR)
r = system.file(package="psidR")
f = fread(file.path(r,"psid-lists","famvars.txt"))
i = fread(file.path(r,"psid-lists","indvars.txt"))

# alternatively, use getNamesPSID:
# cwf <- read.xlsx("http://psidonline.isr.umich.edu/help/xyr/psid.xlsx")
# Suppose you know the name of the variable in a certain year, and it is
# "ER17013". then get the correpsonding name in another year with
# getNamesPSID("ER17013", cwf, years = 2001)  
# getNamesPSID("ER17013", cwf, years = 2003)
# getNamesPSID("ER17013", cwf, years = NULL)
# getNamesPSID("ER17013", cwf, years = c(2005, 2007, 2009))

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
          SAScii = TRUE, 
          heads.only = TRUE,
          sample="SRC",
          design=2)
```

### Usage

#### In case you go for psidR option 1 

* download the zipped family data from [http://simba.isr.umich.edu/Zips/ZipMain.aspx](http://simba.isr.umich.edu/Zips/ZipMain.aspx)
  * run any of the contained program statements in each of the downloaded folders
* download the cross-year individual file
* the user can set some sample design options
* subsetting criteria
* if some variables are not measured in a given wave for whatever reason, the package takes care of that (after you tell it which ones are missing. see examples in package).

#### If you go for psidR option 2

You don't have to prepare anything: just enough time (you should think about leaving your machine on over night/the weekend, depending on how many waves you want to use. The individual index file is very big).

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

#### Usage Outline

Suppose the user wants to have a panel with variables "house value", "total income" and "education" covering years 2001 and 2003. Steps 1 and 2 are relevant only for **option 1**, **option 2** requires only step 3 and 4:

1. Download the zipped family files and cross-period individual files from [http://simba.isr.umich.edu/Zips/ZipMain.aspx](http://simba.isr.umich.edu/Zips/ZipMain.aspx), best into the same folder. This folder will be the function argument `datadir`.
2. inside each downloaded folder, run the stata, sas or spss routine that comes with it. Fixes the text file up into a rectangular dataset. Save the data as either .dta or .csv. The default of the package requires that you use file names **FAMyyyy.dta** and **IND2009ER.dta** (case sensitive). 
3. Supply a data.frame **fam.vars** which contains the variable names for each wave from the family file.
4. Supply a data.frame **ind.vars** which contains the variable names for each wave from the individual index file.

```r
myvars <- data.frame(year=c(2001,2003),
                       house.value=c("ER17044","ER21043"),
                       total.income=c("ER20456","ER24099"),
                       education=c("ER20457","ER24148"))
indvars1 = data.frame(year=c(2001,2003),longitud.wgt=c("ER33637","ER33740"))
```

5. call the function, with `SAScii=TRUE` or `SAScii=FALSE` depending on your choice:

```r
option.1 <- build.panel(datadir=mydir,fam.vars=myvars,ind.vars=indvars,SAScii=FALSE)
option.2 <- build.panel(datadir=mydir,fam.vars=myvars,ind.vars=indvars,SAScii=TRUE)
```


Stata users may recognize this syntax from module [psiduse](http://ideas.repec.org/c/boc/bocode/s457040.html), which is similar. The names are up to you ("house.value" is your choice), but the rest is not, i.e. there must be a column "year". Notice if you knew house.value was missing in year 2001, you could account for that with 

```r
fam.vars <- data.frame(year=c(2001,2003),
                       house.value=c(NA,"ER21043"),
                       total.income=c("ER20456","ER24099"),
                       education=c("ER20457","ER24148"))
```

The function will then keep NA as the value of the variable in year 2001 and you can fix this later on. This functionality was needed because NAs have a generic meaning, i.e. a person who does not participate in a given year is kept in the register, but has no replies in the family file, so has NA in all variables of the family file after merging.


### Supplemental Datasets

The PSID has a wealth of add-on datasets. Once you have a panel those are easy to merge on. The panel will have a variable `interview`, which is the identifier in the supplemental dataset. 

### Additional Info

* Please check out [the R survey package](http://cran.r-project.org/web/packages/survey/index.html) for analyzing complex survey's with R. 
* Also go to [http://www.asdfree.com/](http://www.asdfree.com/) for a range of tutorials and tips for using survey data with R.


### Future Developments

* allow more complex panel designs, like accounting for wider family structure (i.e. using the family splitoff indicator to follow households that split up).


