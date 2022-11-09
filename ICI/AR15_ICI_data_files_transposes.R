### After Review Khyati ICI data transpose

mydf92 = dcast(hd2,Sample+Disease.status~ASV,value.var = "ICI")
write_csv(mydf92, "Age_2_5_ici_transpose.csv")

mydf95 = dcast(hd5,Sample+Disease.status~ASV,value.var = "ICI")
write_csv(mydf95, "Age_5_ici_transpose.csv")
