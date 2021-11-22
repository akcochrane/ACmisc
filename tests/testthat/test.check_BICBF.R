context("check_BICBF")

suppressWarnings({
library(ACmisc)
library(lme4)
})

### ###
test_that('BICBF finds delta-BIC BFs for an lmer',{
  expect_equivalent({
    
    m <- lmer(ravens ~ alert * corsi  + (1 | isAdult), dat_cochraneEtAl_2019_PLOSOne)
    
    m_bf <- BICBF(m)
    
    round(m_bf$BFlog3,2)},
    c(NA, 0.45, 6.63, 0.87) # entered manually on 2021-11-10
    )
})

### ###
test_that('BICBF finds delta-BIC BFs for a glmer',{
  expect_equivalent({
    
    suppressWarnings({
    m <- glmer(ravens ~ alert * corsi + (1 | isAdult)
               , dat_cochraneEtAl_2019_PLOSOne
               ,family = binomial)
    
    
    m_bf <- BICBF(m)
    })
    
    round(m_bf$BFlog3,2)},
    c(NA, -0.20,  3.19, -0.04) # entered manually on 2021-11-10
  )
})

### ###
test_that('BICBF finds delta-BIC BFs for a rlm',{
  expect_equivalent({
    
    suppressWarnings({
      m <- rlm(ravens ~ alert * corsi + isAdult
                 , dat_cochraneEtAl_2019_PLOSOne
                 )
      
      
      m_bf <- BICBF(m)
    })
    
    round(m_bf$BFlog3,2)},
    c(NA, -1.47,  5.12, 23.42, -1.00) # entered manually on 2021-11-10
  )
})

### ###
test_that('BICBF finds delta-BIC BFs for a glm',{
  expect_equivalent({
    
    suppressWarnings({
      m <- glm(ravens ~ alert * corsi + isAdult
               , dat_cochraneEtAl_2019_PLOSOne
               ,family = binomial)
      
      
      m_bf <- BICBF(m)
    })
    
    round(m_bf$BFlog3,2)},
    c(NA, 0,  0, 0, 0) # need to be entered after fixing the function
  )
})

### ###
test_that('BICBF finds delta-BIC BFs for a lm',{
  expect_equivalent({
    
    suppressWarnings({
      m <- lm(ravens ~ alert * corsi + isAdult
               , dat_cochraneEtAl_2019_PLOSOne
               )
      
      
      m_bf <- BICBF(m)
    })
    
    round(m_bf$BFlog3,2)},
    c(NA, -1.43,  5.20, 23.59, -1.02) # need to be entered after fixing the function
  )
})
