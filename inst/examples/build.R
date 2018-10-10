

# build PSID panel data
library(psidR)
f = fread("famvars.txt")
i = fread("indvars.txt")

f[1:37,vgroup := "wage"]
f[38:74,vgroup := "earnings"]
setkey(f,vgroup)

i[1:38, vgroup := "age"]
i[39:76, vgroup := "educ"]  # caution about 2 first years: no educ data
i[77:114, vgroup := "weight"]
setkey(i,vgroup)

ind = cbind(i[J("age"),list(year,age=variable)],i[J("educ"),list(educ=variable)],i[J("weight"),list(weight=variable)])
fam = cbind(f[J("wage"),list(year,wage=variable)],f[J("earnings"),list(earnings=variable)])

d = build.panel(fam.vars=fam,ind.vars=ind,heads.only = TRUE,sample="SRC",design=2)
