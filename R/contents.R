#' HODC
#' @description
#' a dataset from literature, containing 7 groups higher-order drug combination test with 22 unique drugs and 5 different cell lines for 441 tests
#' @source description
#'
HODC = data(HODC)

# https://stats.stackexchange.com/questions/41443/how-to-actually-plot-a-sample-tree-from-randomforestgettree

#' find the previous decision conditions in a tree
#' @importFrom "randomForest" "getTree"
#' @param tree A single decision tree from a randomforest model
#' @param i The row number of a terminal node in a tree dataframe
#'
#' @return A list containing the decision condition and the row number of the previous node
#' @export
#'
#' @examples
#' \dontrun{
#' feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
#' f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
#' Train = HODC$HODC_1[which(HODC$HODC_1$Order>=2),]
#' rf_mtry8 <- randomForest(f , data = Train, importance = TRUE,mtry = 8)
#' tree_1 = as.data.frame(getTree(rf_mtry8,1))
#' prevCond(tree_1,which(tree_1$status== -1)[1])
#' }
prevCond<-function(tree,i){
  if(i %in% tree$right_daughter){
    id<-which(tree$right_daughter==i)
    cond<-paste(tree$split_var[id],">",tree$split_point[id])
  }
  if(i %in% tree$left_daughter){
    id<-which(tree$left_daughter==i)
    cond<-paste(tree$split_var[id],"<",tree$split_point[id])
  }
  return(list(cond=cond,id=id))
}


#' extract decision routes from a tree
#' @importFrom "randomForest" "getTree"
#' @param tree A single decision tree from a randomforest model
#'
#' @return A list of decision routes from the tree
#' @export
#'
#' @examples
#' \dontrun{
#' feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
#' f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
#' Train = HODC$HODC_1[which(HODC$HODC_1$Order>=2),]
#' rf_mtry8 <- randomForest(f , data = Train, importance = TRUE,mtry = 8)
#' tree_1 = as.data.frame(getTree(rf_mtry8,1))
#' tree_1_conds = getConds(tree_1)
#' }
getConds<-function(tree){
  #store all conditions into a list
  colnames(tree)<-sapply(colnames(tree),collapse)
  conds<-list()
  #start by the terminal nodes and find previous conditions
  id.leafs<-which(tree$status==-1)
  j<-0
  for(i in id.leafs){
    j<-j+1
    prevConds<-prevCond(tree,i)
    conds[[j]]<-prevConds$cond
    while(prevConds$id>1){
      prevConds<-prevCond(tree,prevConds$id)
      conds[[j]]<-paste(conds[[j]]," & ",prevConds$cond)
    }
    if(prevConds$id==1){
      conds[[j]]<-paste(conds[[j]]," => ",tree$prediction[i])
    }

  }

  return(conds)
}


#' Replace " " with "_" in a string
#'
#' @param x A string
#'
#' @return A string replaced  " " by "_"
#' @export
#'
#' @examples
#' x="A B"
#' collapse(x)
collapse<-function(x){
  x<-sub(" ","_",x)
  return(x)
}


#' find split variables in a single decision route
#' @importFrom "randomForest" "getTree"
#' @param x A single decision route
#' @param tree A single decision tree from a randomforest model
#'
#' @return characters of split variables in a single decision route
#' @export
#'
#' @examples
#' \dontrun{
#' feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
#' f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
#' Train = HODC$HODC_1[which(HODC$HODC_1$Order>=2),]
#' rf_mtry8 <- randomForest(f , data = Train, importance = TRUE,mtry = 8)
#' tree_1 = as.data.frame(getTree(rf_mtry8,1))
#' tree_1_conds = getConds(tree_1)
#' tree_1_sv = conds_to_route(tree_1_conds[[1]],tree_1)
#' }
conds_to_route = function(x,tree){
  unlist(strsplit( x," "))[unlist(strsplit( x," ")) %in% unique(na.omit(as.character(tree[,3])))]
}

#' Find tree decision routes and split variables inside routes, for a given variable pair, return depth of variable nodes, co-occurrence, and depth difference
#' @importFrom "randomForest" "getTree"
#' @param tree A single decision tree from a randomforest model
#' @param x_pair_wanted A character string of two node variable connected with "_"
#'
#' @return A list of tree decision routes,split variables inside routes, dataframe of tree coeffients for a variable pair, names of the variable pair and tree depth
#' @export
#'
#' @examples
#' \dontrun{
#' feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
#' f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
#' Train = HODC$HODC_1[which(HODC$HODC_1$Order>=2),]
#' rf_mtry8 <- randomForest(f , data = Train, importance = TRUE,mtry = 8)
#' tree_1 = as.data.frame(getTree(rf_mtry8,1))
#' tree_intercoff(tree_1 ,"3_6")
#' }
tree_intercoff = function(tree, x_pair_wanted){
  x_pair_wanted = unlist(strsplit(x_pair_wanted,"_"))
  x_1= x_pair_wanted[1] ;  x_2= x_pair_wanted[2]

  # 这个函数被嵌套在tree_to_allpair_coff里面，已经colnames被整过一次了，所以这里不需要
  # colnames(tree)<-sapply(colnames(tree),collapse)
  conds<-getConds(tree)
  routes = (lapply(conds, function(x) rev(conds_to_route(x,tree)) ))
  # 可以除去重复的
  routes = routes[!duplicated(routes)]
  # 这个树一共有6个路径
  # 找最大子树
  routes_tree_depth = lapply(routes,length)
  tree_depth = max(unlist(routes_tree_depth))
  # 给两个组合，然后查共发生没有，以及相差的深度
  # unique(na.omit(as.character(tree$split_var)))
  x_1_depth = lapply(routes, function(x) which( x %in% x_1))
  x_1_depth =lapply(x_1_depth, function(x) ifelse(identical(x, integer(0))==T,0,x))
  x_1_depth = unlist(x_1_depth)

  x_2_depth = lapply(routes, function(x) which(x%in% x_2  ))
  x_2_depth = lapply(x_2_depth, function(x) ifelse(identical(x, integer(0))==T,0,x))
  x_2_depth = unlist(x_2_depth)

  # 判断在每一条Route里是否同时出现
  co_occ_log = lapply(routes, function(x) ifelse((x_1 %in% x ) & (x_2 %in% x ),1,0))
  co_occ_log = unlist(co_occ_log)
  # 可以直接sum，也可以只取一个，只要全部非0就可以

  x1_occ_log = lapply(routes, function(x) ifelse((x_1 %in% x ) ,1,0))
  x1_occ_log = unlist(x1_occ_log )

  x2_occ_log = lapply(routes, function(x) ifelse((x_2 %in% x ) ,1,0))
  x2_occ_log = unlist(x2_occ_log )

  var_diff_len = ifelse(co_occ_log == 0 , 0,x_2_depth-x_1_depth)
  # 注意这里可能有正有负
  # 有可能全部都是0 要写一个判断，避免全是0的时候出现Inf
  var_diff_len_min = ifelse(sum(  var_diff_len) <= 0 ,0,min((var_diff_len[which(var_diff_len > 0)])))

  # 分有0 没0 两种情况
  #
  # ifelse(
  #   length(var_diff_len[which(var_diff_len == 0)]) == 0,
  #   min( var_diff_len[which(var_diff_len>0)])
  # )
  tree_coff_df =data.frame(
    x_1_depth =  x_1_depth ,
    x_2_depth =  x_2_depth  ,
    co_occ_log =co_occ_log   ,
    x1_occ_log = x1_occ_log ,
    x2_occ_log  = x2_occ_log,
    var_diff_len =var_diff_len,
    var_diff_len_min  = var_diff_len_min
  )
  corr_var=c(x_1,x_2)
  # x_1_depth
  # x_2_depth
  # co_occ_log
  # x1_occ_log
  # x2_occ_log
  # var_diff_len
  # var_diff_len_min
  return(tree_coff = list(conds = conds,
                          routes = routes,
                          tree_coff_df =  tree_coff_df,
                          corr_var =  corr_var,
                          tree_depth = tree_depth
  ))
}


#' Given a tree, for all possible variable pair, extract tree decision routes and split variables inside routes, depth of variable nodes, co-occurrences, and depth differences
#' @importFrom "randomForest" "getTree"
#'
#' @param tree A single decision tree from a randomforest model
#'
#' @return A big list containing coefficients for all possible variable pairs, namely the tree decision routes and split variables inside routes, depth of variable nodes, co-occurrences, and depth differences
#'
#' @export
#'
#' @examples
#' \dontrun{
#' feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
#' f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
#' Train = HODC$HODC_1[which(HODC$HODC_1$Order>=2),]
#' rf_mtry8 <- randomForest(f , data = Train, importance = TRUE,mtry = 8)
#' tree_1 = as.data.frame(getTree(rf_mtry8,1))
#' tree_to_allpair_coff(tree_1)
#' }
tree_to_allpair_coff =function(tree){
  colnames(tree)<-sapply(colnames(tree),collapse)
  # print(  colnames(tree_ana))
  var_pool <- unique(na.omit(as.character(tree[,3])))
  # x_pair = lapply(as.data.frame(combn(1:length(var_pool ), 2)), function(x) var_pool[x])
  # 再反过来
  # 加了这一段，得到很多空值
  # 原因是忘加了上面两句
  x_pair = lapply(as.data.frame( cbind(combn(1:length(var_pool ), 2),
                                       combn(1:length(var_pool ), 2)[c(2,1),]
  ) ), function(x) var_pool[x])
  names(x_pair) = unlist(lapply(x_pair ,function(x)(paste0(x,collapse = "_"))))
  tree_allpair_coff = lapply(x_pair, function(x_pair_wanted) tree_intercoff(tree, x_pair_wanted))
  return(tree_allpair_coff)
}

#' For a varaiable pair calculate TBRF coefficients based on tree coefficients list
#' @importFrom "randomForest" "getTree"
#' @param tree_coff_list A tree coefficients list returned by tree_to_allpair_coff()
#' @param Tar_pair A character string of two node variable connected with "_"
#'
#' @return A data frame containing basic coefficients to calculated TBRF score
#' @export
#'
#' @examples
#' \dontrun{
#' feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
#' f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
#' Train = HODC$HODC_1[which(HODC$HODC_1$Order>=2),]
#' rf_mtry8 <- randomForest(f , data = Train, importance = TRUE,mtry = 8)
#' tree_1 = as.data.frame(getTree(rf_mtry8,1))
#' tree_cofflist = tree_to_allpair_coff(tree_1 )
#' coff_TBRFA (tree_cofflist,"3_6")
#' }
coff_TBRFA = function(tree_coff_list, Tar_pair){
  O_ab =  sum(tree_coff_list[(names(tree_coff_list)  == Tar_pair)][[1]]$tree_coff_df$ co_occ_log)
  O_ab = ifelse(!O_ab == 0 , 1 ,0)

  O_a =  sum(tree_coff_list[(names(tree_coff_list)   == Tar_pair)][[1]]$tree_coff_df$ x1_occ_log)
  O_a = ifelse(!O_a == 0 , 1 ,0)

  O_b =  sum(tree_coff_list[(names(tree_coff_list)   == Tar_pair)][[1]]$tree_coff_df$ x2_occ_log)
  O_b = ifelse(!O_b == 0 , 1 ,0)


  D_T  = tree_coff_list[(names(tree_coff_list)  == Tar_pair)][[1]]$tree_depth

  # 这里有出现负数，就会导致D_ab出现无限值，需要重新定义一哈
  D_ab = ifelse(O_ab == 0 | (!length(which( tree_coff_list[(names(tree_coff_list)  == Tar_pair)][[1]]$tree_coff_df$ var_diff_len<0))==0),
                0,
                min(tree_coff_list[(names(tree_coff_list)  == Tar_pair)][[1]]$tree_coff_df$ var_diff_len[
                  which((tree_coff_list[(names(tree_coff_list)  == Tar_pair)][[1]]$tree_coff_df$ var_diff_len > 0))
                ]))

  TBRFA_df = data.frame(Tar_pair = Tar_pair , O_ab=O_ab,O_a=O_a,O_b=O_b, D_T= D_T, D_ab= D_ab)
  return(TBRFA_df = TBRFA_df)
}

#' Given a tree list of randomforest model,calculated the TBRF score of every variable pair
#' @importFrom "randomForest" "getTree"
#' @param rf_tree_list A list of trees from a random forest model
#'
#' @return A data.frame containing TBRF score of all variable pair for a randomforest model
#' @export
#'
#' @examples
#' \dontrun{
#' feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
#' f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
#' Train = HODC$HODC_1[which(HODC$HODC_1$Order>=2),]
#' rf_model <- randomForest(f , data = Train, importance = TRUE,mtry = 8)
#' rf_ntree = rf_model$ntree
#' find_k_tree = function(x) {tree_to_allpair_coff(getTree(rf_model, k=x, labelVar=TRUE))}
#' rf_tree_list = lapply(1: rf_ntree, find_k_tree)
#' tree_to_TBRFA(rf_tree_list)
#' }
tree_to_TBRFA =function(rf_tree_list){
  TBRFA_df =  do.call(rbind,lapply(1:length(rf_tree_list),function(y)
    do.call(rbind,lapply(1:length(rf_tree_list [[y]]),
                         function(x) coff_TBRFA(rf_tree_list [[y]],
                                                names(rf_tree_list [[y]])[x]))
    )
  ))
  TBRFA_df$Tree_Num = unlist(lapply(1:length(rf_tree_list), function(x)rep(x,
                                                                           length(names(rf_tree_list[[x]]) )
  )))
  return(  TBRFA_df =  TBRFA_df)
}


#' Given a tree list of randomforest model,calculated the TI score of every variable pair
#' @importFrom "randomForest" "getTree"
#' @param rf_tree_list A list of trees from a random forest model
#'
#' @return A data.frame containing TI score of all variable pair for a randomforest model
#' @export
#'
#' @examples
#' \dontrun{
#' feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
#' f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
#' Train = HODC$HODC_1[which(HODC$HODC_1$Order>=2),]
#' rf_model <- randomForest(f , data = Train, importance = TRUE,mtry = 8)
#' rf_ntree = rf_model$ntree
#' find_k_tree = function(x) {tree_to_allpair_coff(getTree(rf_model, k=x, labelVar=TRUE))}
#' rf_tree_list = lapply(1: rf_ntree, find_k_tree)
#' tree_to_TI(rf_tree_list)
#' }
tree_to_TI = function(rf_tree_list){
  TI_df =  do.call(rbind,lapply(1:length(rf_tree_list),function(y) #length(rf_tree_list)
    do.call(rbind,lapply(1:length(rf_tree_list [[y]]),
                         function(x) {
                           temp_coff = rf_tree_list [[y]][x][[1]]$tree_coff_df
                           temp_coff$pair_name  = names(rf_tree_list [[y]][x])
                           temp_coff$tree_depth = rf_tree_list [[y]][x][[1]]$tree_depth
                           temp_coff$routesNum  = nrow( rf_tree_list [[y]][x][[1]]$tree_coff_df)
                           temp_coff$Tree_Num   = y
                           rbind( temp_coff)
                         }
    )

    )
  ))
  # unique( TI_df $Tree_Num)

  return(  TI_df )
}

#' Calculate TI score for a randomforest model
#' @importFrom "randomForest" "getTree"
#' @importFrom "philentropy" "KL"
#' @importFrom "stats" "na.omit"
#' @importFrom "utils" "combn"
#' @param rf_model A randomforest model
#'
#' @return A data.frame containing TI score for all variable pair
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' feats <- (colnames(HODC$HODC_1)[4:9]);f <- paste(feats,collapse=' + ')
#' f <- paste('PhenoRespon ~',f);f <- as.formula(f);f
#' Train = HODC$HODC_1[which(HODC$HODC_1$Order>=2),]
#' rf_model <- randomForest(f , data = Train, importance = TRUE,mtry = 8)
#' RF2TI(rf_model)
#' }
RF2TI = function(rf_model){
  rf_ntree = rf_model$ntree
  rf_tree_list = lapply(1: rf_ntree, function(x) tree_to_allpair_coff(getTree(rf_model, k=x, labelVar=TRUE)))
  TI_df =  tree_to_TI(rf_tree_list)
  pair_KL_l = lapply(1:length(unique(TI_df$pair_name)), function(y) {
    pair_name = unique(TI_df$pair_name)[y]
    TI_df_temp = TI_df[which(TI_df$pair_name == unique(TI_df$pair_name)[y]),]
    # 有时候长度不相等
    Ori_dist = table(TI_df_temp [-which(TI_df_temp $x_1_depth==0),] $x_1_depth)/sum(table(TI_df_temp[-which(TI_df_temp $x_1_depth==0),]  $x_1_depth))
    Opt_dist = table(TI_df_temp[which(TI_df_temp$co_occ_log == 1),] $x_1_depth)/sum(table(TI_df_temp[which(TI_df_temp$co_occ_log == 1),] $x_1_depth))
    if ( !length(Ori_dist) == length(Opt_dist)){
      if (length(Ori_dist) > length(Opt_dist)){
        Opt_dist [(length(Opt_dist)+1):length(Ori_dist)] = 0
      }
      if (length(Ori_dist) < length(Opt_dist)){
        Ori_dist [(length(Ori_dist)+1):length(Opt_dist)] = 0
      }
    }
    x = rbind(Ori_dist,
              Opt_dist
    )
    # KL(x, unit='log')
    # KL(x, test.na = TRUE, est.prob = "empirical")
    KL_pair = KL(x)
    KL_df = data.frame(pair_name = pair_name ,KL_value=KL_pair)
    return(KL_df )
  })
  pair_KL_df = do.call(rbind,pair_KL_l)

  return( TL_list = list(
    TI_df     = TI_df ,
    KL_df     = pair_KL_df
  ))
}


