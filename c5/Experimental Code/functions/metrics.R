# RSQ (R-squared)
RSQfunc = function(predicted, actual){
  rsq = 1 - ( (sum((actual-predicted)^2))  / 
                (sum((actual-mean(actual))^2))
  ) 
  return(round(rsq*100,3))
}

# RMSE (Root mean squared error)
RMSEfunc = function(predicted, actual, days){
  z = (predicted - actual)^2
  return(sqrt(sum(z)/days))
}

# MAPE (Mean absolute percentage error)
MAPEfunc = function(predicted,actual){
  return(mean(abs((actual-predicted)/actual))*100)
}