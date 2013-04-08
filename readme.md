
# psidR: make building panel data from PSID easy

This R package provides a function to easily build panel data from PSID raw data.

## PSID

The [Panel Study of Income Dynamics](http://psidonline.isr.umich.edu/) is a publicly available dataset. You have to register and agree to terms and conditions, but there are no other strings attached. 

* you can use the [data center](http://simba.isr.umich.edu/default.aspx) for simple datasets
* not workable for larger datasets
  * some variables don't show up (although you know they exist)
  * the ftp interface gets slower the more periods you are looking at
  * the click and scroll exercise of selecting the right variables in each period is extremely error prone. 
* merging the data manually is tricky.

## psidR

this package attempts to help the task of building a panel data. the user has 
* to specify the variable names in each wave of the questionnaire
* download the zipped family data from [http://simba.isr.umich.edu/Zips/ZipMain.aspx](http://simba.isr.umich.edu/Zips/ZipMain.aspx)
  * run any of the contained program statements in each of the downloaded folders
* download the cross-year individual file
* **all to the same location**


