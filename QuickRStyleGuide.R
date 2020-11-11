###R Style Guide for Charlotte
###I'm sure this is alot of review for you, but here are some notes on differences between R and Python

################################# Basic R Functions ###############################################
# Indexing in R starts at 1, not 0. Idk why.
# '<-' is R's assignment operatior. However, '=' is also R's assignment operator, so use them when you want
# c(1, 2, 3) is the same as [1, 2, 3] in python
# ':' is a sequence operator -- 
1:3 == c(1, 2, 3)
# skipping a dimension in when indexing a matrix implies all elements in that dimension
# matrix[,c(1, 2)] returns the first two columns of the matrix

# R does not assign variables by reference, so while in python the following code would be true
i = 1
j = i
j = 3
i == 3
# In R, i would remain 1

################################# Data structures #################################################
#Matrices have to hold data of the same type
#Data Frames can have columns of different types of data
df <- data.frame(characters = c("One", "Two", "Three"), 
                 numeric = c(1, 2, 3))

#Named columns of a dataframe can be refereced using the $ operator
df$characters == c('One', 'Two', 'Three')


################################# Pacakge Managers in R ###########################################
#1) install.packages('your package') is in-built and generally my first go to
#2) Bioconductor is a bioinformatics focused repository of packages <-- many big bioinformatic tools are only listed in bioconductor

##To install bioconductor##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

##To install packages with bioconductor
BiocManager::install('splatter')

################################# Loading packages ###########################################
#If you want to access a function without loading a package, use :: operator
stats::median(c(1, 2, 3, 4, 5))

#R Version of Import is library
library(splatter)

#Once you import a package, all functions/methods are available in your name space
library(stats)
stats::median(c(1, 2, 3, 4, 5)) == median(c(1, 2, 3, 4, 5))

#If you have imported two packages each with a function with identical names, R will use the
#provided arguments to automatically select the correct funcion

##For example given functions 
#A::multiply(number1, number2, decimalPlace)
#B::multiply(number1, number2, returnMessage)

#The code multiply(1, 2, 'done!') will automatically evaluate B's multiply

#In addition, R will ~generally~ figure out your arguments without specifying the exact argument that your variable is passed to
#So, if you do multiply('done!', 1, 2), R will ~generally~ map the string 'done!' to variable returnMessage
#But it is good practice to be explicit in your variable mapping

#multiply(number1 = 1,
#         number2 = 2,
#         returnMessage = 'Done!')
