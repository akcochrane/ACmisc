

context("check_YeoJohn")

library(ACmisc)

### ###
test_that('YeoJohn de-skews a variable',{
  expect_equivalent({
    round(
      psych::skew(
        YeoJohn(ACmisc::dat_cochraneEtAl_2019_PLOSOne$enum)
    )
  ,2)},
    0)
})
