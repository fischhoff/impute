DF <- read.csv("../../functions/pima-indians-diabetes.csv")
DF$id_field_vals = seq(1,dim(DF)[1])
nrounds_test = 2000
nrounds = nrounds_test
min_trees = 0
k_split = 0.8
label = "X1"
eta = 0.5
max_depth = 1
vars_pred = setdiff(names(DF), label)
vars_pred = setdiff(vars_pred, "id_field_vals")
vars = vars_pred
file_label = "test"
bootstrap = 1
n.minobsinnode = 2
id_field = "id_field_vals" 
i = 1
distribution = "bernoulli"
cv.folds = 5
OUT <-  bootstrapGBM(DF = DF, label = label, vars = vars_pred, k_split = k_split, distribution = "bernoulli", eta = eta, max_depth = 1, nrounds = nrounds_test, nruns =2, bootstrap = "observed", method = "cv", cv.folds = 5,
                     n.minobsinnode=2,file_label = "test", id_field= id_field)

temp = OUT[[3]]
str(temp)
