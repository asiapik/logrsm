source('./logRSM.R')

csv = read.csv('./data.csv', header = FALSE)
class = csv[, 1]
data = csv[, 3:ncol(csv)]

test_simple = function() {
  reg = logRSM(class, data, m = 5, initial_weights = calculate_weighss_with_univ, stopControl = logRSMStop(B = 1))
  result = rev(colnames(data)[order(reg$scores)])
  print (result[1:10])
  print(reg$control$initial_weights[1:20])
  print(reg$stop_reason)
  print(reg$ns)
}

test_simple()
