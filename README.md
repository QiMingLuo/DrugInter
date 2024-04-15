
## A demo using DrugInter
DrugInter is an R package dissecting pair-wise variable interaction strength from higher order combination data based on random forest, the predictive performance of ranking drug pair has been evaluated on the High Order Drug Combination data set(HODC). Here is a example using DrugInter to rank potent interactive variable pair.

For users ready to perform wet lab combination experiments for a certain initial drug set, the high order combination index could be guided by an orthogonal array returned by design_high_comb() function.
```{r eval=FALSE}
# Please install DrugInter package first.
# devtools::install_github("QiMingLuo/DrugInter")
# library(DrugInter)
InitialDrugNames = c("A","B","C","D","E","F","G")
CombArray = design_high_comb(InitialDrugNames)
print(CombArray)
```
design_high_comb() function use a vector of initial drug names as input. For a initial set with drug number of n, ideally this step reduce wet lab test number from all possible pair n*(n-1)/2 to (n+1).

After obtaining web lab data, users could build random forest model using RandomForest R package. Customized training is allowed, users could choose croos-validation folds and hyperparameters according to their own experiments. As far as the final model is a randomforest object returned by RandomForest::randomforest() function. Here, taking HODC dataset as an example:

```{r eval=FALSE}
# HODC contain 7 groups of high order drug combination
# to see the details please refer to our manuscript
data(HODC);names(HODC)
head(HODC$HODC_1)

# We use higher order combination wet lab data as input to build model
TrainSet = HODC$HODC_1[which(HODC$HODC_1$Order >= 3),]

feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
# You can customize the training procedure, use the most proper hyperparameters for your own data
rf_mtry6 <- randomForest::randomForest(f , data = TrainSet, importance = TRUE,mtry = 6);rf_mtry6

# NOTE: 
# HODC contain true single and pair wise biological phenotype data
# If you want to find out new predictive coefficients
# this part of data should be kept out of training procedure 
# and be used as test set to calculate true interaction strength 
# TestSet = HODC$HODC_1[which(HODC$HODC_1$Order  <= 2),]
```

Once random forest model is built, pair-wise interaction coefficients could be directly obtained by RF2TI() function 

```{r eval=FALSE}
# This may take a while
All_pair_TI = RF2TI(rf_mtry6)
str(All_pair_TI )
# List of 2

#  $ TI_df:'data.frame':	91102 obs. of  11 variables:
#   ..$ x_1_depth       : num [1:91102] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ x_2_depth       : num [1:91102] 2 2 0 2 3 3 3 3 0 0 ...
#   ..$ co_occ_log      : num [1:91102] 1 1 0 1 1 1 1 1 0 0 ...
#   ..$ x1_occ_log      : num [1:91102] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ x2_occ_log      : num [1:91102] 1 1 0 1 1 1 1 1 0 0 ...
#   ..$ var_diff_len    : num [1:91102] 1 1 0 1 2 2 2 2 0 0 ...
#   ..$ var_diff_len_min: num [1:91102] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ pair_name       : chr [1:91102] "navitoclax_trametinib" "navitoclax_trametinib" "navitoclax_trametinib" "navitoclax_trametinib" ...
#   ..$ tree_depth      : int [1:91102] 5 5 5 5 5 5 5 5 5 5 ...
#   ..$ routesNum       : int [1:91102] 8 8 8 8 8 8 8 8 8 8 ...
#   ..$ Tree_Num        : int [1:91102] 1 1 1 1 1 1 1 1 1 1 ...
#  $ KL_df:'data.frame':	30 obs. of  2 variables:
#   ..$ pair_name: chr [1:30] "navitoclax_trametinib" "navitoclax_PF04217903" "navitoclax_ribociclib" "navitoclax_alpelisib" ...
#   ..$ KL_value : num [1:30] 0 0 0 0 0 ...
```

And then ranked by rank_drug_pair() function.
```{r eval=FALSE}

pair_TI_rank = rank_drug_pair(All_pair_TI$KL_df)
head(pair_TI_rank)
# Higher KL_value suggests stronger interaction

# RowNumber           pair_name             KL_value
# kullback-leibler13  ribociclib_erlotinib  0.3798053
# kullback-leibler10  PF04217903_alpelisib  0.3453138
# kullback-leibler11  PF04217903_erlotinib  0.3181804
# kullback-leibler9   PF04217903_ribociclib 0.2079004
# kullback-leibler28  erlotinib_ribociclib  0.1769805
# kullback-leibler26  erlotinib_PF04217903  0.1246209

```

A demo ranking drug pair by interactive coefficients is done. For developers, the logic of calculating interactive coefficient is to dig model structure first. Forest structure could be extracted by get_route_depth() function.
```{r eval=FALSE}
# This may take a while
Basic_Coeff_List = get_route_depth(rf_mtry6)

# get_route_depth function return a big list
#
# The structure contains 3 levels
# level 1|--- A tree level
# level 2|----- A variable-pair level under a tree
# level 3|------- A detailed distribution level for variable pair in every decision route under a tree
#
# for every decision tree, return a list containing details for every possible variable pair 
# for every possible variable pair in a decision tree, return 5 lists
# List of 5
# $ conds       :List of 13
# $ routes      :List of 8
# $ tree_coff_df:'data.frame':	8 obs. of  7 variables:
# $ corr_var    : chr [1:2] "navitoclax" "trametinib"
# $ tree_depth  : int 5

# List 1 : conds
# Containing every single decision condition in the decision tree
# $ conds       :List of 13
# ..$ : chr "ribociclib < 0.5  &  trametinib < 0.5  &  navitoclax < 0.5  =>  0.0375558418153902"
# ..$ : chr "ribociclib > 0.5  &  trametinib < 0.5  &  navitoclax < 0.5  =>  0.0448454958588662"
# ..$ : chr "alpelisib < 0.5  &  trametinib > 0.5  &  navitoclax < 0.5  =>  0.0663340773164321"
# ..$ : chr "ribociclib < 0.5  &  PF04217903 < 0.5  &  navitoclax > 0.5  =>  0.620057518675292"
# ..$ : chr "ribociclib > 0.5  &  PF04217903 < 0.5  &  navitoclax > 0.5  =>  0.825206059642883"

# List 2 : routes
# Containing the variables split route in the decision tree
# $ routes      :List of 8
# ..$ : chr [1:3] "navitoclax" "trametinib" "ribociclib"
# ..$ : chr [1:3] "navitoclax" "trametinib" "alpelisib"
# ..$ : chr [1:3] "navitoclax" "PF04217903" "ribociclib"
# ..$ : chr [1:4] "navitoclax" "trametinib" "alpelisib" "erlotinib"

# List 3 : tree_coff_df
# A variable pair interaction frame containing the variables' depth, 
# co-occurrence, single occurrence, depth difference, min depth difference in every routes.
# NOTE: This is a main core data.frame
# Users can calculate new coefficients based on this data.frame
# x_1_depth        describe the variable's position (depth) in the routes for the first appearance,and so is x_2_depth
# x1_occ_log       describe whether a specific variable appears in a decision route, 1 means appear and 0 means not. So is       x2_occ_log
# co_occ_log       describe whether these two variable both appears in a decision route, 1 means two variable both appear.
# var_diff_len     describe the position (depth) difference of two variable in a decision route
# var_diff_len_min describe the min position (depth) difference of two variable for all routes in a tree
#
# $ tree_coff_df:'data.frame':	8 obs. of  7 variables:
# ..$ x_1_depth       : int [1:8] 1 1 1 1 1 1 1 1
# ..$ x_2_depth       : num [1:8] 2 2 0 2 3 3 3 3
# ..$ co_occ_log      : num [1:8] 1 1 0 1 1 1 1 1
# ..$ x1_occ_log      : num [1:8] 1 1 1 1 1 1 1 1
# ..$ x2_occ_log      : num [1:8] 1 1 0 1 1 1 1 1
# ..$ var_diff_len    : num [1:8] 1 1 0 1 2 2 2 2
# ..$ var_diff_len_min: num [1:8] 1 1 1 1 1 1 1 1

# List 4 : corr_var
# The variable names
# $ corr_var    : chr [1:2] "navitoclax" "trametinib"

# List 5 : tree_depth
# The depth of a tree, representing the max length of a route in this tree
# $ tree_depth  : int 5
```

And then use get_pair_inter() function to transfer structure details into interaction coefficients.
```{r eval=FALSE}
KL_df = get_pair_inter(Basic_Coeff_List )
```

Finally, rank drug pair.
```{r eval=FALSE}
rank_drug_pair(KL_df)
```

Notably, DrugInter could handle omics data and be used to find interactive genes, metabolites and microbes.
