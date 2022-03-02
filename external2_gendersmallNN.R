
######################   SUBPROGRAMS DISCRIMINATOR #########################

# loss function
Loss<-function(table_aux){
  miin = 0.01
  maax = 1-0.01
  pred_true<-sapply(table_aux[d ==1]$V1, function(x) min(max(x,miin),maax))
  pred_false<-sapply(table_aux[d ==0]$V1, function(x) min(max(x,miin),maax))
  N <- length(pred_true)
  H <- length(pred_false)
  loss_f<-1/N*sum(log(pred_true))+1/H*sum(log(1-pred_false))
  return (loss_f)
}



#define the network 
define_model <- function(){
  model_keras <- keras_model_sequential()
  
  model_keras %>%
    
	    # First hidden layer
    layer_dense(
      units              = 20,
      kernel_initializer = "uniform",
      activation         = "sigmoid",
      #kernel_regularizer = regularizer_l1_l2(l1 = 0, l2 = 0.1),
      input_shape        = c(20),
      use_bias=TRUE,
      name = "input") %>%
    
    # Dropout to prevent overfitting
    #layer_dropout(rate = 0.1,
    #             name = "dropbout1") %>%
    
    # Second hidden layer
    #layer_dense(
    #units              = 10,
    #kernel_initializer = "uniform",
    #activation         = "sigmoid",
    #kernel_regularizer = regularizer_l1_l2(l1 = 0, l2 = 0.1),
    #name = "hid1") %>%

    # Dropout to prevent overfitting
    #layer_dropout(rate = 0.1,
    #             name = "drop2") %>%
     
  # Third hidden layer
   #layer_dense(
  #   units              = 10,
  #   kernel_initializer = "uniform",
  #   activation         = "sigmoid",
  #   name = "hid2") %>%
  
      # Dropout to prevent overfitting
  # layer_dropout(rate = 0.5,
  #               name = "drop3") %>%
  
  # Second hidden layer
  #layer_dense(
  #  units              = 2,
  #  kernel_initializer = "uniform",
  #  activation         = "sigmoid",
  #  name = "hid3") %>%
  
  # Dropout to prevent overfitting
  #layer_dropout(rate = 0.1,
  #              name = "drop4") %>%
  
  # Output layer
  layer_dense(
    units              = 1,
    kernel_initializer = "uniform",
    activation         = "sigmoid",
    name = "hid4") %>%
    
    # Compile ANN
    compile(
      optimizer = 'adam',
      loss      = 'binary_crossentropy',
      metrics   = c('accuracy')
    )
  return(model_keras)
}
# 

train_model <- function(x_train, y_train,epoch_p){
  
  x_mat_train <- as.matrix(x_train)
  y_mat_train <- as.matrix(y_train)
  #start_time <- Sys.time()
  history <- fit(
    object           = model_keras, 
    x                = x_mat_train, 
    y                = y_mat_train,
    batch_size       = 240, 
    epochs           = epoch_p,
    verbose = 0,
    shuffle=TRUE,
    validation_split = .2
    )
  #end_time <-Sys.time()
  #print('Network Training time')
  #print(end_time-start_time)
  score <- model_keras %>% evaluate(x_mat_train, y_mat_train)
  #print(score$acc)
  return(score$acc)
  
}


train_model2 <- function(x_train, y_train, x_test, y_test, epoch_p,model_keras,file_name){
  
  x_mat_train <- as.matrix(x_train)
  y_mat_train <- as.matrix(y_train)
  x_mat_test <- as.matrix(x_test)
  y_mat_test <- as.matrix(y_test)
  
  
  #start_time <- Sys.time()
  #es = callback_early_stopping(monitor='loss', min_delta=0, patience=50, verbose=0, mode='auto', baseline=NULL,restore_best_weights=FALSE)
  #mcp = callback_model_checkpoint("/Experiments 9 parameters - Copy", monitor='acc', verbose=0, save_best_only=TRUE, save_weights_only=FALSE, mode='auto', period=1)
  
  history <- fit(
    object           = model_keras, 
    x                = x_mat_train, 
    y                = y_mat_train,
    batch_size       = 249, 
    epochs           = epoch_p,
    validation_data = list(x_mat_test, y_mat_test),
    callbacks        = list(callback_early_stopping(monitor = "loss",patience = 100,verbose = 0), callback_model_checkpoint(file_name,monitor='loss',save_best_only = TRUE,verbose = 0)),
#			    callback_lambda(on_train_end  = keras::k_clear_session())),
    verbose          = 0
  )
  #end_time <-Sys.time()
  #print('Network Training time')
  #print(end_time-start_time)
  #score <- model_keras %>% evaluate(x_mat_train, y_mat_train)
  #print(score$acc)
  #return()
#  keras::k_clear_session()
}


#define the network 
define_model_nohealth <- function(){
  model_keras <- keras_model_sequential()
  
  model_keras %>%
    
    # First hidden layer
    layer_dense(
      units              = 20,
      kernel_initializer = "uniform",
      activation         = "relu",
      input_shape        = c(13),
      use_bias=TRUE,
      name = "input") %>%
    
    # Dropout to prevent overfitting
    #layer_dropout(rate = 0.1,
    #              name = "dropbout1") %>%
    
    # Second hidden layer
    layer_dense(
      units              = 10,
      kernel_initializer = "uniform",
      activation         = "sigmoid",
      name = "hid1") %>%
    
    # Dropout to prevent overfitting
    #layer_dropout(rate = 0.1,
    #              name = "drop2") %>%
  
  # Third hidden layer
  #layer_dense(
  #  units              = 100,
  #  kernel_initializer = "uniform",
  #  activation         = "relu",
  #  name = "hid2") %>%
  
  # Dropout to prevent overfitting
  #layer_dropout(rate = 0.1,
  #              name = "drop3") %>%
  
  # Second hidden layer
  #layer_dense(
  #  units              = 2,
  #  kernel_initializer = "uniform",
  #  activation         = "sigmoid",
  #  name = "hid3") %>%
  
  # Dropout to prevent overfitting
  #layer_dropout(rate = 0.1,
  #              name = "drop4") %>%
  
  # Output layer
  layer_dense(
    units              = 1,
    kernel_initializer = "uniform",
    activation         = "sigmoid",
    name = "hid4") %>%
    
    # Compile ANN
    compile(
      optimizer = 'adam',
      loss      = 'binary_crossentropy',
      metrics   = c('accuracy')
    )
  return(model_keras)
}
# 




#define the network 
define_model_x <- function(){
  model_keras <- keras_model_sequential()
  
  model_keras %>%
    
    # First hidden layer
    layer_dense(
      units              = 50,
      kernel_initializer = "uniform",
      activation         = "relu",
      input_shape        = c(14),
      use_bias=TRUE,
      name = "input") %>%
    
    # Dropout to prevent overfitting
    layer_dropout(rate = 0.1,
                  name = "dropbout1") %>%
    
    # Second hidden layer
    #layer_dense(
    #  units              = 20,
    #  kernel_initializer = "uniform",
    #  activation         = "relu",
    #  name = "hid1") %>%
    
    # Dropout to prevent overfitting
    #layer_dropout(rate = 0.1,
    #              name = "drop2") %>%
  
  # Third hidden layer
  #layer_dense(
  #  units              = 100,
  #  kernel_initializer = "uniform",
  #  activation         = "relu",
  #  name = "hid2") %>%
  
  # Dropout to prevent overfitting
  #layer_dropout(rate = 0.1,
  #              name = "drop3") %>%
  
  # Second hidden layer
  #layer_dense(
  #  units              = 2,
  #  kernel_initializer = "uniform",
  #  activation         = "sigmoid",
  #  name = "hid3") %>%
  
  # Dropout to prevent overfitting
  #layer_dropout(rate = 0.1,
  #              name = "drop4") %>%
  
  # Output layer
  layer_dense(
    units              = 1,
    kernel_initializer = "uniform",
    activation         = "sigmoid",
    name = "hid4") %>%
    
    # Compile ANN
    compile(
      optimizer = 'adam',
      loss      = 'binary_crossentropy',
      metrics   = c('accuracy')
    )
  return(model_keras)
}
# 



train_model_new <- function(x_train, y_train,epoch_p,model_keras){
  
  x_mat_train <- as.matrix(x_train)
  y_mat_train <- as.matrix(y_train)
  #start_time <- Sys.time()
  history <- fit(
    object           = model_keras, 
    x                = x_mat_train, 
    y                = y_mat_train,
    batch_size       = 1000, 
    epochs           = epoch_p,
    verbose = 0
  )
  #end_time <-Sys.time()
  #print('Network Training time')
  #print(end_time-start_time)
  score <- model_keras %>% evaluate(x_mat_train, y_mat_train)
  #print(score$acc)
  return(score$acc)
  
}

