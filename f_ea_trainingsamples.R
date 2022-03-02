training_samples <- function(data_D,split_vec){
  
  data_D2 <- data_D[,split_vec:=split_vec]

  data_D_tr <- data_D2[split_vec==1]
  x_train_d1 <- data_D_tr[d==1,.(gender,PI,age_1996,asset_1996,asset_1998,asset_2000,asset_2002,asset_2004,asset_2006,alive_1998,alive_2000,alive_2002,alive_2004,alive_2006,health_1996,health_1998, health_2000, health_2002, health_2004,health_2006)]
  x_train_d0 <- data_D_tr[d==0,.(gender,PI,age_1996,asset_1996,asset_1998,asset_2000,asset_2002,asset_2004,asset_2006,alive_1998,alive_2000,alive_2002,alive_2004,alive_2006,health_1996,health_1998, health_2000, health_2002, health_2004,health_2006)]
  
  y_train_d1 <- as.data.table(data_D_tr[d==1]$d)
  y_train_d0 <- as.data.table(data_D_tr[d==0]$d)
  
  x_train <- rbind(x_train_d1,x_train_d0)
  y_train <- rbind(y_train_d1,y_train_d0)
  setnames(y_train, "d")
  
  data_D_val <- data_D2[split_vec==0]
  x_val_d1 <- data_D_val[d==1,.(gender,PI,age_1996,asset_1996,asset_1998,asset_2000,asset_2002,asset_2004,asset_2006,alive_1998,alive_2000,alive_2002,alive_2004,alive_2006,health_1996,health_1998, health_2000, health_2002, health_2004,health_2006)]
  x_val_d0 <- data_D_val[d==0,.(gender,PI,age_1996,asset_1996,asset_1998,asset_2000,asset_2002,asset_2004,asset_2006,alive_1998,alive_2000,alive_2002,alive_2004,alive_2006,health_1996,health_1998, health_2000, health_2002, health_2004,health_2006)]
  
  y_val_d1 <- as.data.table(data_D_val[d==1]$d)
  y_val_d0 <- as.data.table(data_D_val[d==0]$d)
  
  x_val <- rbind(x_val_d1,x_val_d0)
  y_val <- rbind(y_val_d1,y_val_d0)
  setnames(y_val, "d")
  
  
  outcome <-list(x_train= x_train, y_train=y_train, x_val=x_val, y_val=y_val)
  return(outcome)
}

training_samples_nohealthg <- function(data_D,split_vec){
  
  data_D2 <- data_D[,split_vec:=split_vec]
  
  data_D_tr <- data_D2[split_vec==1]
  x_train_d1 <- data_D_tr[d==1,.(PI,age_1996,asset_1996,asset_1998,asset_2000,asset_2002,asset_2004,asset_2006,alive_1998,alive_2000,alive_2002,alive_2004,alive_2006)]
  x_train_d0 <- data_D_tr[d==0,.(PI,age_1996,asset_1996,asset_1998,asset_2000,asset_2002,asset_2004,asset_2006,alive_1998,alive_2000,alive_2002,alive_2004,alive_2006)]
  
  y_train_d1 <- as.data.table(data_D_tr[d==1]$d)
  y_train_d0 <- as.data.table(data_D_tr[d==0]$d)
  
  x_train <- rbind(x_train_d1,x_train_d0)
  y_train <- rbind(y_train_d1,y_train_d0)
  setnames(y_train, "d")
  
  data_D_val <- data_D2[split_vec==0]
  x_val_d1 <- data_D_val[d==1,.(PI,age_1996,asset_1996,asset_1998,asset_2000,asset_2002,asset_2004,asset_2006,alive_1998,alive_2000,alive_2002,alive_2004,alive_2006)]
  x_val_d0 <- data_D_val[d==0,.(PI,age_1996,asset_1996,asset_1998,asset_2000,asset_2002,asset_2004,asset_2006,alive_1998,alive_2000,alive_2002,alive_2004,alive_2006)]
  
  y_val_d1 <- as.data.table(data_D_val[d==1]$d)
  y_val_d0 <- as.data.table(data_D_val[d==0]$d)
  
  x_val <- rbind(x_val_d1,x_val_d0)
  y_val <- rbind(y_val_d1,y_val_d0)
  setnames(y_val, "d")
  
  
  outcome <-list(x_train= x_train, y_train=y_train, x_val=x_val, y_val=y_val)
  return(outcome)
}