
# psidR: make building panel data from PSID easy

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

this package attempts to help the task of building a panel data. the user has 
* to specify the variable names in each wave of the questionnaire
* download the zipped family data from [http://simba.isr.umich.edu/Zips/ZipMain.aspx](http://simba.isr.umich.edu/Zips/ZipMain.aspx)
  * run any of the contained program statements in each of the downloaded folders
* download the cross-year individual file
* the user can set some sample design options
* subsetting criteria
* if some variables are not measured in a given wave for whatever reason, the package takes care of that (after you tell it which ones are missing. see examples in package).

### Example Usage

Suppose the user wants to have a panel with variables "house value", "total income" and "education" covering years 2001 and 2003. Here are the steps to take:

1. Download the zipped family files and cross-period individual files from [http://simba.isr.umich.edu/Zips/ZipMain.aspx](http://simba.isr.umich.edu/Zips/ZipMain.aspx), best into the same folder.
2. Supply a data.frame **fam.vars** which contains the variable names for each wave from the family file.

```r
fam.vars <- data.frame(year=c(2001,2003),
                       house.value=c("ER17044","ER21043"),
                       total.income=c("ER20456","ER24099"),
                       education=c("ER20457","ER24148"))
```

Stata users may recognize this syntax from module [psiduse](http://ideas.repec.org/c/boc/bocode/s457040.html), which is similar.

3. Specify options for the panel, like *design* or *heads.only*
4. call the function **build.panel**
5. the result is a wide data.table where the id colums are *pid* (person identifier) and *year*. 


### Future Developments

* allow to merge [supplemental datasets](http://simba.isr.umich.edu/Zips/zipSupp.aspx), not only family files.


