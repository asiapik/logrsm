source('./logRSM.R')

csv = read.csv('./data.csv', header = FALSE)
class = csv[, 1]
data = csv[, 3:ncol(csv)]

test_simple = function() {
  reg = logRSM(class, data, m = 5, stopControl = logRSMStop(min_ns = 20, max_time = 2))
  result = rev(colnames(data)[order(reg$scores)])
  print (result[1:10])
}

test_simple()
