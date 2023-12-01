# Creating a stump for a tree
stump <- function(data,
                  j){

  # Creating the base node
  node <- list()


  node[["node0"]] <- list(
    # Creating the node number
    node_number = 0,
    j = j,
    inter = NA,
    isRoot = TRUE,
    # Creating a vector with the tranining index
    train_index = 1:nrow(data$x_train),
    test_index = 1:nrow(data$x_test),
    depth_node = 0,
    node_var = NA,
    node_cutpoint_index = NA,
    left = NA,
    right = NA,
    parent_node = NA,
    ancestors = NA,
    terminal = TRUE,
    betas_vec = NULL
  )

  # Returning the node
  return(node)

}

# Get all the terminal nodes
get_terminals <- function(tree){

  # Return the name of the termianl nodes
  return(names(tree)[unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)])
}

# Get nog terminal nodes
get_nogs <- function(tree){

  # Return the name of the termianl nodes
  non_terminal <- names(tree)[!unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)]

  # In case there are non nonterminal nondes
  if(length(non_terminal)==0){
    return(non_terminal)
  }

  bool_nog <- vector("logical",length = length(non_terminal))
  for(i in 1:length(bool_nog)){
    # Checking if both children are terminal
    if( tree[[tree[[non_terminal[i]]]$left]]$terminal & tree[[tree[[non_terminal[i]]]$right]]$terminal) {
      bool_nog[i] <- TRUE
    }
  }

  return(  non_terminal[bool_nog])
}

# Getting the maximum node index number
get_max_node <- function(tree){

  # Return the name of the termianl nodes
  return(max(unlist(lapply(tree, function(x){x$node_number}),use.names =  TRUE)))
}




# A function to calculate the loglikelihood
nodeLogLike <- function(curr_part_res,
                        j_,
                        index_node,
                        data){

  # Subsetting the residuals
  curr_part_res_leaf <- curr_part_res[index_node]

  # Getting the number of observationsin the terminal node
  n_leaf <- length(index_node)
  d_basis <- length(j_)
  ones <- matrix(1,nrow = n_leaf)
  # Getting the index of all covariates for each basis from each ancestors
  D_subset_index <- unlist(data$basis_subindex[j_])


  # Getting the p_{tell} i.e: number of betas of the current terminal node
  d_basis <- length(D_subset_index)
  D_leaf <- data$D_train[index_node,D_subset_index, drop = FALSE]

  if(NCOL(D_leaf)==0){
    stop(" Node Log-likelihood: No variables")
  }

  # Using the Andrew's approach I would have
  mean_aux <- rep(0,length(curr_part_res_leaf))
  diag_tau_beta_inv <- diag(x = 1/unique(data$tau_beta), nrow = NCOL(D_leaf))


  if(data$dif_order==0){
    # P_aux <- data$P[D_subset_index,D_subset_index]
    # for(jj in 1:length(j_)){
    #   P_aux[data$basis_subindex[[jj]],data$basis_subindex[[jj]]] <- data$tau_beta[jj]*data$P[data$basis_subindex[[jj]],data$basis_subindex[[jj]]]
    # }
      cov_aux <- diag(x = (data$tau^(-1)),nrow = n_leaf) + D_leaf%*%tcrossprod(diag_tau_beta_inv,D_leaf)
  } else {
      P_aux <- data$P
      for(jj in 1:length(j_)){
        P_aux[data$basis_subindex[[jj]],data$basis_subindex[[jj]]] <- data$tau_beta[jj]*data$P[data$basis_subindex[[jj]],data$basis_subindex[[jj]]]
      }
      P_aux <- data$P[D_subset_index,D_subset_index]

      cov_aux <- diag(x = (data$tau^(-1)),nrow = n_leaf) + D_leaf%*%solve(P_aux,t(D_leaf))

  }


  result <- mvnfast::dmvn(X = curr_part_res_leaf,mu = mean_aux,
                          sigma = cov_aux ,log = TRUE)


  return(c(result))

}



# Grow a tree
grow <- function(tree,
                 curr_part_res,
                 data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0

  # acceptance_grid <- numeric(100)
  # for(kk in 1:100){
  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2

    # Sample a split var
    # ===== Uncomment this line below after ========
    p_var <- sample(1:NCOL(data$x_train),size = 1)
    # ==============================================
    # p_var <- 7

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[g_node$train_index,p_var])

    # Case of invalid range
    if(length(valid_range_grow)==0){
      return(tree)
    }

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)
    # sample_cutpoint <- valid_cutpoint[kk]

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% g_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% g_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% g_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% g_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(g_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(g_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }


    # === Uncomment those lines after
  if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
    # Getting out of the while
    break
  } else {

    # Adding one to the counter
    valid_count = valid_count + 1

    # Stop trying to search for a valid cutpoint
    if(valid_count > 2) {
      valid_terminal_node = FALSE
      return(tree)
    }
  }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }

  # Calculating loglikelihood for the grown node, the left and the right node
  # Recover the g_node index

  if(!any(is.na(g_node$inter))){
    node_index_var <- c(g_node$j,which( names(data$basis_subindex) %in% paste0(g_node$j,sort(g_node$inter))))
  } else {
    node_index_var <- g_node$j
  }

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = node_index_var,
                           index_node = g_node$train_index,
                           data = data)


  left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               j_ = node_index_var,
                               index_node = left_index,
                               data = data)

  right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                j_ = node_index_var,
                                index_node = right_index,
                                data = data)

  # Calculating the prior
  prior_loglike <- log(data$alpha*(1+g_node$depth_node)^(-data$beta)) + # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+g_node$depth_node+1)^(-data$beta)) - # plus the prior of the two following nodes being terminal
    log(1-data$alpha*(1+g_node$depth_node)^(-data$beta)) # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_nog_nodes+1))-log(0.3/n_t_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+left_loglike+right_loglike+prior_loglike+log_trasition_prob)
  # acceptance_grid[kk] <- acceptance



  # par(mfrow=c(1,2))
  # plot(data$xcut_m[,2],(acceptance_grid), main = "Acceptance to split on X2", xlab = "X2", ylab = "Prob. Acceptance")

  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

        if(any(is.na(g_node$ancestors))){
          new_ancestors <- p_var
        } else {
          new_ancestors <- c(g_node$ancestors,p_var)
        }

        left_node <- list(node_number = max_index+1,
                          j = g_node$j,
                          inter = g_node$inter,
                          isRoot = FALSE,
                          train_index = left_index,
                          test_index = left_test_index,
                          depth_node = g_node$depth_node+1,
                          node_var = p_var,
                          node_cutpoint_index = sample_cutpoint,
                          left = NA,
                          right = NA,
                          parent_node = g_node_name,
                          ancestors = new_ancestors,
                          terminal = TRUE,
                          betas_vec = rep(0,ncol(data$D_train)))

        right_node <- list(node_number = max_index+2,
                           j = g_node$j,
                           inter = g_node$inter,
                           isRoot = FALSE,
                           train_index = right_index,
                           test_index = right_test_index,
                           depth_node = g_node$depth_node+1,
                           node_var = p_var,
                           node_cutpoint_index = sample_cutpoint,
                           left = NA,
                           right = NA,
                           parent_node = g_node_name,
                           ancestors = new_ancestors,
                           terminal = TRUE,
                           betas_vec = rep(0,ncol(data$D_train)))

    # Modifying the current node
    tree[[g_node_name]]$left = paste0("node",max_index+1)
    tree[[g_node_name]]$right = paste0("node",max_index+2)
    tree[[g_node_name]]$terminal = FALSE

    tree[[paste0("node",max_index+1)]] <- left_node
    tree[[paste0("node",max_index+2)]] <- right_node


  } else {

    # Do nothing

  }

  # Return the new tree
  return(tree)
}


# Adding interaction
add_interaction <- function(tree,
                 curr_part_res,
                 data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  # g_node_name <- "node2"
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0

  # Maybe need to change this when we start to work with continuous variables
  interaction_candidates <- (1:NCOL(data$x_train))[-g_node$j]

  while(length(interaction_candidates)!=0){

      # Sample a split var
      # ===== Uncomment this line below after ========
      p_var <- sample(interaction_candidates,size = 1)
      p_var <- 2
      # ==============================================
      # Getting the interaction name
      int_name_ <- paste0(sort(c(p_var,g_node$j)),collapse = '')

      # Making sure that not selecting a interaction that's already in the terminal node
      if(any(is.na(g_node$inter))){
        break
      } else {
        if(p_var %in% g_node$inter){
          index_candidate <- which(interaction_candidates %in% p_var)
          interaction_candidates <- interaction_candidates[-index_candidate]
        } else {
          break # Considering the case that the variable isnt there
        }
      }

  }

  # All interactions were already included.
  if(length(interaction_candidates)==0){
    return(tree)
  }

  # Calculating loglikelihood for the grown node, the left and the right node

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_= g_node$j, # Here j_ refers to the main effect from that tree
                           index_node = g_node$train_index,
                           data = data)

  # Adding the interaction term on the list
  new_j <- c(which(names(data$basis_subindex) %in% g_node$j),which(names(data$basis_subindex) %in% int_name_))
  # new_j <- c(10)
  new_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               j_ = new_j,
                               index_node = g_node$train_index,
                               data = data)



  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+new_loglike)
  # acceptance_grid[kk] <- acceptance



  # par(mfrow=c(1,2))
  # plot(data$xcut_m[,2],(acceptance_grid), main = "Acceptance to split on X2", xlab = "X2", ylab = "Prob. Acceptance")

  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  if(g_node$j %in% interaction_candidates){
    Stop("Cannot have the main effect variable into the candidates of the model")
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Modifying the current node (case is the first interaction)
    if(any(is.na(tree[[g_node_name]]$inter))){
      tree[[g_node_name]]$inter <- p_var
    } else {
      tree[[g_node_name]]$inter = c(tree[[g_node_name]]$inter,p_var) # Case of previous interactions
    }

  } else {
    # Do nothing

  }

  # Return the new tree
  return(tree)
}


# Pruning a tree
prune <- function(tree,
                  curr_part_res,
                  data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)

  # Just in case to avoid errors
  if(n_nog_nodes==0){
    return(tree)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(nog_nodes,size = 1)
  p_node <- tree[[p_node_name]]

  # Getting the indexes from the left and right children from the pruned node
  children_left_index <- tree[[p_node$left]]$train_index
  children_right_index <- tree[[p_node$right]]$train_index
  children_left_ancestors <- tree[[p_node$left]]$ancestors
  children_right_ancestors <- tree[[p_node$right]]$ancestors

  # Calculating loglikelihood for the grown node, the left and the right node

  if(!any(is.na(p_node$inter))){
    node_index_var <- c(p_node$j,which( names(data$basis_subindex) %in% paste0(p_node$j,sort(p_node$inter))))
  } else {
    node_index_var <- p_node$j
  }

  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           index_node = p_node$train_index,
                           j_ = node_index_var,
                           data = data)


  p_left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                 index_node =  children_left_index,
                                 j_ = node_index_var,
                                 data = data)

  p_right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                  index_node = children_right_index,
                                  j_ = node_index_var,
                                  data = data)

  # Calculating the prior
  prior_loglike <- log(1-data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the new terminal node
    log(data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+p_node$depth_node+1)^(-data$beta))  # plus the prior of the two following nodes being terminal
  # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_t_nodes))-log(0.3/n_nog_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(p_loglike-p_left_loglike-p_right_loglike+prior_loglike+log_trasition_prob)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[p_node$left]] <- NULL
    tree[[p_node$right]] <- NULL

    # Modifying back the pruned node
    tree[[p_node_name]]$left <- NA
    tree[[p_node_name]]$right <- NA
    tree[[p_node_name]]$terminal <- TRUE

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

# Pruning a tree
prune_interaction <- function(tree,
                  curr_part_res,
                  data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)

  t_with_inter <- names(which(sapply(terminal_nodes,function(node){!all(is.na(tree[[node]]$inter))})))

  if(length(t_with_inter)==0){
    # to avoid to waste an interaction
    std_prune <- prune(tree = tree,
                       curr_part_res = curr_part_res,data = data)
    return(std_prune)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(t_with_inter,size = 1)
  p_node <- tree[[p_node_name]]


  # Calculating loglikelihood for the grown node, the left and the right node
  if(!any(is.na(p_node$inter))){
    node_index_var <- c(p_node$j,which( names(data$basis_subindex) %in% paste0(p_node$j,sort(p_node$inter))))
    inter_index_ <- p_node$inter

      # Sampling the new interactions subset
      if(length(inter_index_)==1){
        p_inter_index <- p_node$inter
        new_p_inter <- NA
        new_node_index_var <- p_node$j
        # new_node_index_var <- c(p_node$j,which( names(data$basis_subindex) %in% paste0(p_node$j,sort(new_p_inter))))

      } else  {
        new_p_inter <- sort(sample(p_node$inter,size = length(p_node$inter)-1,replace = FALSE))
        new_node_index_var <- c(p_node$j,which( names(data$basis_subindex) %in% paste0(p_node$j,sort(new_p_inter)))) #
      }

  } else {
    stop('Prune interaction was called where there is no interaction')
  }

  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = node_index_var,
                           index_node = p_node$train_index,
                           data = data)


  new_p_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                 j_ = new_node_index_var,
                                 index_node = p_node$train_index,
                                 data = data)


  # Calculating the acceptance probability
  acceptance <- exp(-p_loglike+new_p_loglike)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
   tree[[p_node_name]]$inter <- new_p_inter

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

change_stump <- function(tree = tree,
             curr_part_res = curr_part_res,
             data = data){

  # Getting the stump
  c_node <- tree$node0

  # Proposing a change to the stump
  change_candidates <- which(!(1:NCOL(data$x_train) %in% c_node$j))

  # In case there's other proposal trees (only for 1-d case)
  if(length(change_candidates)==0){
    # Gettina grown tree
    grown_tree <- grow(tree = tree,
         curr_part_res = curr_part_res,
         data = data)
    return(grown_tree)
  }

  # Print
  new_ancestor <- sample(change_candidates,size = 1)


  stump_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                                     j_ = c_node$j,
                                     index_node = c_node$train_index,
                                     data = data)

  new_stump_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = new_ancestor,
                           index_node = c_node$train_index,
                           data = data)

  acceptance <- exp(-stump_loglikelihood+new_stump_loglikelihood)

  # Update the stump in case of changing the main var
  if(stats::runif(n = 1)<acceptance){
    # Modifying the node0
    tree$node0$j <- new_ancestor
  }


  # Returning the new tree
  return(tree)


}


# Change a tree
change <- function(tree,
                   curr_part_res,
                   data){

  # Changing the stump
  if(length(tree)==1){
    change_stump_obj <- change_stump(tree = tree,
                                 curr_part_res = curr_part_res,
                                 data = data)
    return(change_stump_obj)
  }

  # Sampling a terminal node
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  c_node_name <- sample(nog_nodes,size = 1)
  c_node <- tree[[c_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0


  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2
    # Sample a split var
    p_var <- sample(1:ncol(data$x_train),size = 1)

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[c_node$train_index,p_var])

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% c_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% c_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% c_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% c_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(c_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(c_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }

    # Avoiding having terminal nodes with just one observation
    if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }


  # Getting the node_index var
  if(!any(is.na(c_node$inter))){
    node_index_var <- c(c_node$j,which( names(data$basis_subindex) %in% paste0(c_node$j,sort(c_node$inter))))
  } else {
    node_index_var <- c_node$j
  }

  # Calculating loglikelihood for the new changed nodes and the old ones
  c_loglike_left <- nodeLogLike(curr_part_res = curr_part_res,
                                index_node = tree[[c_node$left]]$train_index,
                                j_ = node_index_var,
                                data = data)


  c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                  index_node = tree[[c_node$right]]$train_index,
                                  j_ =  node_index_var,
                                  data = data)

  # Calculating a new ancestors left and right
  old_p_var <- tree[[c_node$left]]$node_var

  # Storing new left and right ancestors
  new_left_ancestors <- tree[[c_node$left]]$ancestors
  new_left_ancestors[length(new_left_ancestors)] <- p_var

  new_right_ancestors <- tree[[c_node$right]]$ancestors
  new_right_ancestors[length(new_right_ancestors)] <- p_var


  new_c_loglike_left <-  nodeLogLike(curr_part_res = curr_part_res,
                                     index_node = left_index,
                                     j = node_index_var,
                                     data = data)

  new_c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                      index_node = right_index,
                                      j =  node_index_var,
                                      data = data)


  # Calculating the acceptance probability
  acceptance <- exp(new_c_loglike_left+new_c_loglike_right-c_loglike_left-c_loglike_right)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1,min = 0,max = 1)<acceptance){

    # Updating the left and the right node
    # === Left =====
    tree[[c_node$left]]$node_var <- p_var
    tree[[c_node$left]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$left]]$train_index <- left_index
    tree[[c_node$left]]$test_index <- left_test_index
    tree[[c_node$left]]$ancestors <- new_left_ancestors

    #==== Right ====
    tree[[c_node$right]]$node_var <- p_var
    tree[[c_node$right]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$right]]$train_index <- right_index
    tree[[c_node$right]]$test_index <- right_test_index
    tree[[c_node$right]]$ancestors <- new_right_ancestors

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

# Change interaction
change_interaction <-  function(tree,
         curr_part_res,
         data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)

  t_with_inter <- names(which(sapply(terminal_nodes,function(node){!all(is.na(tree[[node]]$inter))})))

  # For the cases where there's no interactions
  if(length(t_with_inter)==0){
    std_change <- change(tree = tree,
                         curr_part_res = curr_part_res,
                         data = data)
    return(std_change)
  }

  # Selecting a node to be pruned
  c_node_name <- sample(t_with_inter,size = 1)
  c_node <- tree[[c_node_name]]


  # Calculating loglikelihood for the grown node, the left and the right node
  if(!any(is.na(c_node$inter))){
    node_index_var <- c(c_node$j,which( names(data$basis_subindex) %in% paste0(c_node$j,sort(c_node$inter))))
    inter_index_ <- c_node$inter

    # Sampling the new interactions subset
    if(length(inter_index_)==1){
      c_inter_index <- c_node$inter
      new_interaction_candidates <- (1:NCOL(data$x_train))[-c(c_node$j,c_node$inter)]

      if(length(new_interaction_candidates)==0){
        stop('Change interaction error')
      }

      new_c_inter <-   sample(new_interaction_candidates,size = 1)
      new_node_index_var <- c(c_node$j,which( names(data$basis_subindex) %in% paste0(c_node$j,sort(new_c_inter))))

    } else  {
      inter_to_change <- sample(c_node$inter,size = 1)
      new_interaction_candidates <- which(!(1:NCOL(data$x_train)) %in% c(c_node$inter,c_node$j))
      if(length(new_interaction_candidates)==0){
        return(tree)
      }
      new_c_inter_single  <- sample(new_interaction_candidates,size = 1)
      inter_to_change_index <- which(c_node$inter %in% inter_to_change)
      new_c_inter <- c_node$inter
      new_c_inter[inter_to_change_index] <- new_c_inter_single
      new_c_inter <- sort(new_c_inter)
      new_node_index_var <- c(c_node$j,which( names(data$basis_subindex) %in% paste0(c_node$j,new_c_inter))) #
    }

  } else {
    stop('Prune interaction was called where there is no interaction')
  }

  c_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = node_index_var,
                           index_node = c_node$train_index,
                           data = data)


  new_c_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                j_ = new_node_index_var,
                                index_node = c_node$train_index,
                                data = data)


  # Calculating the acceptance probability
  acceptance <- exp(-c_loglike+new_c_loglike)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[c_node_name]]$inter <- new_c_inter

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}




# ============
# Update Betas
# ============
updateBetas <- function(tree,
                        curr_part_res,
                        data){


  # Getting the terminals
  t_nodes_names <- get_terminals(tree)


  for(i in 1:length(t_nodes_names)){


    # Select the current terminal node
    cu_t <- tree[[t_nodes_names[i]]]
    # THIS LINE IS COMPLETELY IMPORTANT BECAUSE DEFINE THE ANCESTEORS BY j only

    if(!any(is.na(cu_t$inter))){
      node_index_var <- c(cu_t$j,which( names(data$basis_subindex) %in% paste0(cu_t$j,sort(cu_t$inter))))
    } else {
      node_index_var <- cu_t$j
    }

    res_leaf <- matrix(curr_part_res[cu_t$train_index], ncol=1)

    # Creatinga  vector of zeros for betas_vec
    tree[[t_nodes_names[[i]]]]$betas_vec <- rep(0,ncol(data$D_train))

    # Selecting the actually parameters subsetting
    leaf_basis_subindex <- unlist(data$basis_subindex[unique(node_index_var)]) # Recall to the unique() here too
    basis_dim <- length(leaf_basis_subindex)
    D_leaf <- data$D_train[cu_t$train_index,leaf_basis_subindex, drop = FALSE]
    n_leaf <- length(cu_t$train_index)
    diag_leaf <- diag(nrow = n_leaf)
    diag_basis <- diag(nrow = basis_dim)


    #  Calculating the quantities need to the posterior of \beta
    b_ <- crossprod(D_leaf,res_leaf)
    data_tau_beta_diag <- rep(data$tau_beta[node_index_var], NCOL(D_leaf)) # Don't really use this
    U_ <- data$P
    for(k in 1:length(unique(node_index_var))){
      aux_P_indexes <- unlist(data$basis_subindex[node_index_var[k]])
      U_[aux_P_indexes,aux_P_indexes] <- U_[aux_P_indexes,aux_P_indexes]*(data$tau_beta[node_index_var[k]])
    }
    U_ <- U_[leaf_basis_subindex,leaf_basis_subindex, drop = FALSE]
    U_inv_ <- U_
    Q_ <- (crossprod(D_leaf) + data$tau^(-1)*U_inv_)
    Q_inv_ <- chol2inv(chol(Q_))
    # Q_inv_ <- solve(Q_)

    # tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex] <- c(keefe_mvn_sampler(b = b_,Q = Q_))
    tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex] <- mvnfast::rmvn(n = 1,mu = Q_inv_%*%b_,sigma = (data$tau^(-1))*Q_inv_)
  }

  # Returning the tree
  return(tree)

}


# =================
# Update \tau_betas
# =================
update_tau_betas_j <- function(forest,
                             data){


  # if(data$dif_order!=0){
  #   stop("Do not update tau_beta for peanalised version yet")
  # }

  # Setting some default hyperparameters
  # a_tau_beta <- d_tau_beta <- 0.1
  # Setting some default hyperparameters
  # a_tau_beta <- 0.1
  # d_tau_beta <- 0.1

  # Setting some default hyperparameters
  a_tau_beta <- data$a_tau_beta_j
  d_tau_beta <- data$d_tau_beta_j

  tau_b_shape <- 0.0
  tau_b_rate <- 0.0


  if(data$interaction_term){
    tau_b_shape <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
    tau_b_rate <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
    tau_beta_vec_aux <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
  } else{
    tau_b_shape <- numeric(NCOL(data$x_train))
    tau_b_rate <- numeric(NCOL(data$x_train))
    tau_beta_vec_aux <- numeric(NCOL(data$x_train))
  }

  # Iterating over all trees
  for(i in 1:length(forest)){

    # Getting terminal nodes
    t_nodes_names <- get_terminals(forest[[i]])
    n_t_nodes <- length(t_nodes_names)

    # Iterating over the terminal nodes
    for(j in 1:length(t_nodes_names)){

      cu_t <- forest[[i]][[t_nodes_names[j]]]
      cu_t$ancestors <- cu_t$j
      var_ <- cu_t$ancestors


            # Getting the interactions as well
            if(!any(is.na(cu_t$inter))){
              interaction_index <- cu_t$inter
              interaction_index <- sapply(interaction_index,function(x){sort(c(cu_t$j,x))})
              for(ii in 1:NCOL(interaction_index)){
                var_ <- c(var_,paste0(interaction_index[,ii],collapse = ""))
              }
              # for(var_ in 1:cu_t$ancestors){
              var_ <- which(names(data$basis_subindex) %in% var_)
            }


            # Getting ht leaf basis
            for(kk in 1: length(var_)){
              leaf_basis_subindex <- unlist(data$basis_subindex[var_[kk]]) # Recall to the unique() function here
              p_ <- length(leaf_basis_subindex)
              betas_mat_ <- matrix(cu_t$betas_vec[leaf_basis_subindex],nrow = p_)

              tau_b_shape[var_[kk]] <- tau_b_shape[var_[kk]] + p_
              tau_b_rate[var_[kk]] <- tau_b_rate[var_[kk]] + c(crossprod(betas_mat_,crossprod(data$P[leaf_basis_subindex,leaf_basis_subindex, drop = FALSE],betas_mat_)))

            }
      # }

    }


  }

  if(data$interaction_term){
      for(j in 1:(NCOL(data$x_train)+NCOL(data$interaction_list)) ){
        tau_beta_vec_aux[j] <- rgamma(n = 1,
                                   shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                   rate = 0.5*tau_b_rate[j] + d_tau_beta)

      }
  } else {
      for(j in 1:NCOL(data$x_train)){
        tau_beta_vec_aux[j] <- rgamma(n = 1,
                                      shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                      rate = 0.5*tau_b_rate[j] + d_tau_beta)

      }
  }

  return(tau_beta_vec_aux)

}


update_tau_betas <- function(forest,
                             data){

  if(data$dif_order!=0){
    stop("Do not update tau_beta for peanalised version yet")
  }

  # Setting some default hyperparameters
  a_tau_beta <- d_tau_beta <- 0.1
  tau_b_shape <- 0.0
  tau_b_rate <- 0.0


  # Iterating over all trees
  for(i in 1:length(forest)){

    # Getting terminal nodes
    t_nodes_names <- get_terminals(forest[[i]])
    n_t_nodes <- length(t_nodes_names)

    # Iterating over the terminal nodes
    for(j in 1:length(t_nodes_names)){

      cu_t <- forest[[i]][[t_nodes_names[j]]]
      leaf_basis_subindex <- unlist(data$basis_subindex[unique(cu_t$j)]) # Recall to the unique() function here

      if(!is.null(cu_t$betas_vec)){
        tau_b_shape <- tau_b_shape + length(leaf_basis_subindex)
        tau_b_rate <- tau_b_rate + c(crossprod(cu_t$betas_vec[leaf_basis_subindex]))
      }

    }


    tau_beta_vec_aux <- rgamma(n = 1,
                               shape = 0.5*tau_b_shape + a_tau_beta,
                               rate = 0.5*tau_b_rate + d_tau_beta)
  }

  return(tau_beta_vec_aux)

}


# ===================
# Updating the \delta
# ===================

# A function to get predictions
getPredictions <- function(tree,
                           data){

  # Creating the vector to hold the values of the prediction
  if(data$interaction_term){
    y_hat <- matrix(0, nrow = nrow(data$x_train), ncol = NCOL(data$x_train)+NCOL(data$interaction_list))
    y_hat_test <- matrix(0,nrow(data$x_test), ncol = NCOL(data$x_test)+NCOL(data$interaction_list))
  } else {
    y_hat <- matrix(0, nrow = nrow(data$x_train), ncol = ncol(data$x_train))
    y_hat_test <- matrix(0,nrow(data$x_test), ncol = ncol(data$x_test))
  }

  # Getting terminal nodes
  t_nodes <- get_terminals(tree = tree)
  n_t_nodes <- length(t_nodes)

  for(i in 1:n_t_nodes){


    # Getting the current terminal node
    cu_t <- tree[[t_nodes[[i]]]]
    leaf_train_index <- cu_t$train_index
    leaf_test_index <- cu_t$test_index
    # leaf_ancestors <- unique(tree[[t_nodes[[i]]]]$ancestors) # recall the unique() argument here

    if(!any(is.na(cu_t$inter))){
      node_index_var <- c(cu_t$j,which( names(data$basis_subindex) %in% paste0(cu_t$j,sort(cu_t$inter))))
    } else {
      node_index_var <- cu_t$j
    }

    leaf_ancestors <- node_index_var # here isnt really the ancestors, but the variables that are being used

    leaf_basis_subindex <- data$basis_subindex[leaf_ancestors]

    # Test unit
    if(length(leaf_ancestors)!=length(leaf_basis_subindex)){
      stop("Error on the getPredictions function")
    }

    # Only add the marginal effects if the variables are within that terminal node
    if(length(leaf_basis_subindex)!=0){
      for(k in 1:length(leaf_basis_subindex)){

        y_hat[leaf_train_index,leaf_ancestors[k]] <- y_hat[leaf_train_index,leaf_ancestors[k]] + data$D_train[leaf_train_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]
        y_hat_test[leaf_test_index,leaf_ancestors[k]] <- y_hat_test[leaf_test_index,leaf_ancestors[k]] + data$D_test[leaf_test_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]

      }
    }

  }

  # Returning both training and test set predictions
  return(list(y_train_hat = y_hat,
              y_hat_test = y_hat_test))

}

# Updating tau
update_tau <- function(y_train_hat,
                       data){

  # Sampling a tau value
  n_ <- nrow(data$x_train)
  tau_sample <- stats::rgamma(n = 1,shape = 0.5*n_+data$a_tau,rate = 0.5*crossprod((data$y_train-y_train_hat))+data$d_tau)

  return(tau_sample)

}


