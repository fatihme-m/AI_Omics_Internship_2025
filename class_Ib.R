#### Tasks ####

# 1. Set Working Directory
# Create a new folder on your computer "AI_Omics_Internship_2025".

# 2. Create Project Folder

# In RStudio, create a new project named "Module_I" in your "AI_Omics_Internship_2025" folder.
getwd()  # to check current working directory

setwd("/Users/fatihmemaarawi/Desktop/AI_Omics_Internship_2025/Module_I")

# Inside the project directory, create the following subfolders using R code:
# raw_data, clean_data, scripts, results or Tasks, plots etc

dir.create("raw_data")    # For storing raw data files
dir.create("clean_data") # For storing cleaned data files
dir.create("scripts")   # For saving R scripts
dir.create("results")  # For saving analysis outputs
dir.create("Tasks")   # For saving tasks
dir.create("plots")  # For saving plots and graphs

# 3. Download "patient_info.csv" dataset from GitHub repository
setwd("/Users/fatihmemaarawi/Desktop/AI_Omics_Internship_2025")  

# load the dataset into your R environment
patient_info <- read.csv("patient_info.csv")
View(patient_info)

# Inspect the structure of the dataset using appropriate R functions
str(patient_info) # See structure like data types, column names, etc

head(patient_info) # See first few rows of the data

summary(patient_info) # summary of statistics for each column

names(patient_info) # column names

dim(patient_info) #dimensions: rows and columns

# View the first and last few rows
head(patient_info)
tail(patient_info)

colSums(is.na(patient_info)) # Count of missing values per column


# Identify variables with incorrect or inconsistent data types.
sapply(patient_info, class) #identify the class of each column

char_cols <- sapply(patient_info, is.character) #identify character columns


# Convert variables to appropriate data types where needed
patient_info$gender    <- as.factor(patient_info$gender)    #categorical
patient_info$diagnosis <- as.factor(patient_info$diagnosis) #categorical
patient_info$smoker    <- as.factor(patient_info$smoker)    #categorical 

patient_info$patient_id <- as.character(patient_info$patient_id) #keep patient_id as character (identifier)

patient_info$age <- as.integer(patient_info$age) #ensure age is numeric

patient_info$bmi <- as.numeric(patient_info$bmi) #ensure bmi is numeric

str(patient_info) #check the updated structure


# Create a new variable for smoking status as a binary factor:
# 1 for "Yes", 0 for "No"
patient_info$smoking_status <- ifelse(patient_info$smoker == "Yes", 1, 0)

table(patient_info$smoker, patient_info$smoking_status) #check the results


# Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv
if (!dir.exists("clean_data")) {
  dir.create("clean_data")
}

write.csv(patient_info, "clean_data/patient_info_clean.csv", row.names = FALSE)


# Save your R script in your script folder with name "class_Ib"
if (!dir.exists("script")) {
  dir.create("script")
}

save(patient_info, file = "clean_data/FatihmeMaarawi_Class_1b_Assignment.Rdata")


load("clean_data/patient_info_clean.Rdata")


# Upload "class_Ib" R script into your GitHub repository