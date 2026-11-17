# Create Project Folder

dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")
# load the data set into your R environment

data<- read.csv("raw_data/patient_info.csv")

# Inspect the structure of the data set using appropriate R functions

View(data)
str(data)
head(data)

# Convert variables to appropriate data types where needed

data$gender= as.factor(data$gender)
data$diagnosis= as.factor(data$diagnosis)
data$smoker =ifelse(tolower( data$smoker)== "yes" , TRUE , FALSE)

# Create a new variable for smoking status as a binary factor:
# 1 for "Yes", 0 for "No"


data$smoker_binary = as.factor(ifelse(data$smoker == TRUE , 1 , 0)) 
str(data)

# Save the cleaned data set in your clean_data folder with the name patient_info_clean.csv

write.csv(data, "clean_data/patient_info_clean.csv")

getwd()