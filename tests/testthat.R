## TESTTHAT 
library(testthat)
test_check("openpopscr", filter = "scr")
test_check("openpopscr", filter = "cjs")
test_check("openpopscr", filter = "js")
