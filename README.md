# very_basic_R

very basic R for newbies


soumya
Very basic R

Install R on Ubuntu (link)

sudo apt-get install r-base

other dependencies

From R:

install.packages('devtools')

If there are errors, sudo apt-get install <name debian package>

and try

Install RStudio on Ubuntu (link)

sudo apt-get install gdebi-core

wget https://download2.rstudio.org/server/xenial/amd64/rstudio-server-1.3.1093-amd64.deb sudo gdebi rstudio-server-1.3.1093-amd64.deb

1) R code for basic data munging

# code in R for basic data munging

# inspired by Jeff Leek lecture

# basics of data munging in R

# https://www.youtube.com/watch?v=8DGBOh6hJaE

prismData <- read.csv('combined_hlty_prism_shotgun.csv')

names(prismData)

# names have dot, now split based on dot

splitNames = strsplit(names(prismData),'\\.')

splitNames[[5]]

splitNames[[5]][1]

# apply (map) this function

# lambda x -> x[1] applied to splitNames

# function(x){x[1]}

firstElement <- function(x){x[1]}

sapply(splitNames,firstElement)

head(prismData,3)

# substitute names

tempData <- sapply(splitNames,firstElement)

gsub("_","",tempData)

sub("_","",tempData)

# example to show merge

riskData <- read.csv('combined_hlty_prism_shotgun_MOD.csv')

names(riskData)

head(riskData,3)

merge(prismData,riskData,by.x='cohort',by.y='cohort',all=TRUE)

# only those that match

merge(prismData,riskData,by.x='cohort',by.y='cohort')

# sort values

prismData$cohort

sort(prismData$cohort)

sort(prismData$cohort)[1:10]

# sort data frame

prismData$cohort[1:10]

sort(prismData$cohort[1:10])

2) Install and load packages

install.packages('smatr')

library('smatr')

# Check if package installed, if not then install it (from link)

if(!require("zoo"))

{install.packages("zoo")}

# Install locally

library("highr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")

3) Creating data frames

 log10crime <- read.csv('log10crime.csv')

 log10comm_population <- read.csv('log10comm_population.csv')

 # assign to data frame

 df = data.frame(log10crime,log10comm_population)

4) Extract first few (head) of data frame

head(df)

> log10crime log10comm

1     2.2355    4.0785

2     2.6828    4.3640

# Note: log10crime  and log10commare headers

5) Reference column of data frame by column name

df["log10crime"]

6) Number of rows in data frame

nrow(df)

7) Formula in R (note: sma is a package in SMATR package; load using library(smatr) )

sma(formula='log10crime~log10comm', data=df)

8) Good IDE (RStudio)

9) Help on something

?factor

10) List of numbers

period <- 2003:2012

11) Array index

dat = t(harborSealWA)   # transpose

years = dat[1,]    # 1st row of dat contains years

12) Transpose (from MARSS package)

dat = t(harborSealWA)

13) Access arrays (from MARSS package)

years = dat[1,]    # 1st row of dat contains years

n = nrow(dat)-1

dat = dat[2:nrow(dat),]   # data (log counts), without years

14) Concatenate and print (from MARSS package)

cat("Model 1 AIC:",kem1$AIC,"\n")  #show the AIC

15) Show estimated model coefficients

coef(kem1)

16) Number of columns and rows in array

library(MARSS)

dat = t(harborSealWA)

n = ncol(dat) 

m = nrow(dat)

17) Plotting in R (ggplot2). Good tutorial here and here

# use ggplot2

library(ggplot2)

18) Rename columns

Q_GDP_expenditure_approach <- read.csv("~/Box Sync/machine_learning_workup/forecasting/Code/Q_GDP_expenditure_approach.csv", na.strings="..")

names(Q_GDP_expenditure_approach)[1] <- "Country"

19) Plot numeric data

plot(as.numeric(Q_GDP_expenditure_approach[Q_GDP_expenditure_approach$Country == "Australia",]))

# Labels, titles, etc

plot.forecast(best_arima_forecast, xlab="Years", ylab="Port throughput (AUD, imports)", pch=18, col="blue",main = "Best fit ARIMA model")

20) Time series object

GDP_Australia <- ts(as.numeric(Q_GDP_expenditure_approach[Q_GDP_expenditure_approach$Country ==  "Australia",2:ncol(Q_GDP_expenditure_approach)]),

                                       start=c(1990, 1), frequency=4)

21) Plot multiple time series (use ts and ts.plot)

GDP_Australia <- ts(as.numeric(Q_GDP_expenditure_approach[Q_GDP_expenditure_approach$Country == 

                                                            "Australia",2:ncol(Q_GDP_expenditure_approach)]), start=c(1990, 1), frequency=4)

GDP_Germany <- ts(as.numeric(Q_GDP_expenditure_approach[Q_GDP_expenditure_approach$Country == "Germany",2:ncol(Q_GDP_expenditure_approach)]), start=c(1990, 1), frequency=4)

ts.plot(GDP_Australia, GDP_Germany,

        gpars=list(xlab="year", ylab="GDP", lty=c(1:2)))

title("GDP of Australia and Germany (deseasoned)")

22) Function in R

# function to standardize time series according to MARSS guide

standardize_timeseries <- function(timeseries_vector)

{

  y_bar = apply(t(timeseries_vector), 1, mean, na.rm=TRUE)

  Sigma = sqrt(apply(t(timeseries_vector), 1, var, na.rm=TRUE))

  stdized_GDP_Australia = (timeseries_vector - y_bar) * (1/Sigma)

  return(stdized_GDP_Australia)

}

# call function

TEMP <- as.numeric(Q_GDP_expenditure_approach[Q_GDP_expenditure_approach$Country == "Germany",])

stdized_GDP_Germany = standardize_timeseries(TEMP)

 

23) R code to perform SMA (semi-major axis regression) to test for allometric power-law (on bitbucket)

24) Load function or import function or import R script

source('~/Box Sync/machine_learning_workup/forecasting/Code/MARSS_example_1.R')

25) Concatenate file names and strings

data_path = "~/Documents/imports-exports-forecasting/socioeconomic_data_munging/"

Q_GDP_expenditure_approach <- read.csv(paste(data_path,"Q_GDP_expenditure_approach.csv",sep=""), na.strings="..")

26) Create an R package

package.skeleton()

27) Timeseries analysis using R (link)

28) Timeseries analysis using R (plotting, using scan function and plot.ts)

rain <- scan("http://robjhyndman.com/tsdldata/hurst/precip1.dat",skip=1)

rainseries <- ts(rain,start=c(1813))

plot.ts(rainseries)

29) if statement in R

  if (b_constant_trend == TRUE) {

    if (b_seasonality == FALSE) {

      holtwinters_fit_model <- HoltWinters(timeseries_data_handle, beta = param_beta, gamma = param_gamma)

      # beta = param_beta, gamma = param_gamma should all be FALSE

    }

    else {

      cat("Error: cannot have constant trend and also seasonality")

      b_error = 1

      return (b_error)

    }

# ALSO

if {

    # something

}

} else {

 load(file="imports_quarterly_timeseries_cleaned.RData")

}

# AND if else else if in R (from link)

if (x ==1){

    print('same')

} else if (x > 1){

    print('bigger')

} else {

    print('smaller')

}

30) View command in R (view data frame in a window)

# dat is a data frame

View(dat)

31) Concatenate two or more vectors to create a matrix (also works for time series objects)

combined_stdized_GDP <- cbind(stdized_GDP_Australia, stdized_GDP_China, stdized_GDP_US, 

                              stdized_GDP_Japan, stdized_GDP_Singapore, stdized_GDP_Germany)

32) Save and load data files RData (similar to mat files in MATLAB)

save(imports_quarterly_raw, file="imports_quarterly_raw.RData")

load(file="imports_quarterly_raw.RData")

# save all workspace objects

save.image(file = "all_objects.RData")

33) Load required package

# require(<package_name>)

require(forecast)

34) Return multiple arguments in R function (from stackoverflow)

foo <- 12

bar <- c("a", "b", "e")

newList <- list("integer" = foo, "names" = bar)

# pack all objects to be returned in a list

ret_list <- list("fit_model"=holtwinters_fit_model, "forecast_ts"=forecast_timeseries)

# access using $

35) Remove all workspace variables (similar to clear all in MATLAB)

rm(list = ls())

36) List all objects in current workspace

ls()

37) Time series object refer to one column

target_ts[,"Tot_Value_FOB_AUD"]

38) Shiny R for visualisation (link) (tutorial)

install.packages("shiny")

library(shiny)

runExample("01_hello")

runApp("my_app")

shiny::runApp()

# Another example from the documentation

# Only run these examples in interactive R sessions

if (interactive()) {

# A basic shiny app with a plotOutput

shinyApp(

  ui = fluidPage(

    sidebarLayout(

      sidebarPanel(

        actionButton("newplot", "New plot")

      ),

      mainPanel(

        plotOutput("plot")

      )

    )

  ),

  server = function(input, output) {

    output$plot <- renderPlot({

      input$newplot

      # Add a little noise to the cars data

      cars2 <- cars + rnorm(nrow(cars))

      plot(cars2)

    })

  }

)

39) R cheat sheets for data viz, data wrangling by R studio

40) Time series objects

library(timeSeries)

ts() # creates time series object

length(<time_series_object>)

start(<time_series_object>) # 1995

end(<time_series_object>) # 2014

window(<time_series_object>, start, end) # new time series object with new start and end

frequency(<time_series_object>)

ts.union # union of two different time series with different frequency and start and end dates

ts.intersect # intersection of two different time series with different frequency and start and end dates

interpNA(<time_series_object>, method=“linear”) # interpolate NA in time series with imputed values

<ts_object_1> - <ts_object_2> # can subtract two time series objects

removeNA(forecast_var_model$res) # remove NA in time series

train_start = start(target_ts)

train_end   = c(2007,4) # 2014, 3rd Quarter 

test_end = end(target_ts)

# split into train and test

target_ts_train = window(target_ts, start=train_start, end=train_end)

41) Create a matrix

matrix()

42) String comparison in R

any(grepl("error", (fit_arima_model))

# grepl for grep that returns logical TRUE or FALSE

str_temp = 'caribb'

grepl(pattern = '*bla*|*caribb*', x = str_temp)

regular expressions in grep (link)

# if try-error then TRUE else FALSE

# MORE ADVANCED

# use stringdist package

require(stringdist)

?amatch

43) try in R

fit_arima_model <- try( arima( timeseries_data_handle, order=c(p, d, q),

                                  silent = TRUE )

# also like try catch (from stackoverflow)

retval_try =

        try( ( <statement> ),

           silent = TRUE

         )

if ( inherits(retval_try, 'try-error') )

{

       <error handling statement>

}

44) Socio-economic data from Quandl (link)

install.packages('Quandl')

library(Quandl)

data <- Quandl("FRED/GDP")

head(data)

45) Find minimum element of vector

which(x == min(x)) 

which(x == min(x, na.rm = TRUE)) # if remove NaN

temp_min_index = which(x == min(x, na.rm = TRUE))

x[temp_min_index]

Find index of matches using which

idx_to_remove <- which(d$PC2 > 0.55)[1]

46) Remove elements from vector (stack exchange)

# First approach

a <- sample (1 : 10)

remove <- c(1,2,3,5,17)

a %in% remove

temp_ind = a %in% remove

a[!temp_ind]

# Second approach

setdiff(a, remove)

47) For loop iterating through a sequence of values ( seq command) (stackoverflow)

for (i_temp_counter in seq(1,num_variables)) {

}

48) Time series in R (link)

49) Function to find number of times series has to be differenced (link)

ndiffs

and inverse of differencing (discrete integration) (link)

diffinv

50) Refer to column 

mat[, 1:5]

# OR

mat[, "col_name"]

51) Plot multiple things on the same plot (stackexchange)

plot( x, y1, type="l", col="red" )

par(new=TRUE)

plot( x, y2, type="l", col="green" )

# OR

plot.ts(target_ts_test, xlab ="Years", ylab = "USA GDP",

        main = "Normal scale forecast using VAR model",col = "blue", t="p", pch=19) # plot with dots etc

# OR

plot( x, y1, type="l", col="red" )# plot the actual test time series:

lines(target_ts_test, col='red')

# OR

par(mfrow=c(1,2)) # multi panel plot multi-panel figure

52) Get residuals of model fit

res_model <- residuals(model)

53) Coefficients of model fit

coef(var_model)

54) Histogram

hist(x[,"interp_extr_ts"], 20)

rug(x[,"interp_extr_ts"]) # show the data also

abline(v = 12, col = "magenta", lwd = 4) # a line through the hist

55) Linear model (specify model relationship in R)

lm(y ~ x)

# read file

colitis_score <- read.csv('mouse_colitis_score_extracted.csv')

# fields of data frame

col_score = colitis_score$score

metagene_score <- c(  mean( c( sum(df_il23$d0_R1), sum(df_il23$d0_R2), sum(df_il23$d0_R3), sum(df_il23$d0_R4) ) ),

mean( c( sum(df_il23$d1_R1), sum(df_il23$d1_R2), sum(df_il23$d1_R3), sum(df_il23$d1_R4) ) ),

mean( c( sum(df_il23$d3_R1), sum(df_il23$d3_R2), sum(df_il23$d3_R3), sum(df_il23$d3_R4) ) ),

mean( c( sum(df_il23$d6_R1), sum(df_il23$d6_R2), sum(df_il23$d6_R3), sum(df_il23$d6_R4) ) ),

mean( c( sum(df_il23$d14_R1), sum(df_il23$d14_R2), sum(df_il23$d14_R3), sum(df_il23$d14_R4) ) ),

mean( c( sum(df_il23$d21_R1), sum(df_il23$d21_R2), sum(df_il23$d21_R3), sum(df_il23$d21_R4) ) )

)

# create data frame

df = data.frame(col_score,metagene_score)

library(ggplot2)

# linear model

s = lm(formula = 'col_score ~ metagene_score', df)

# summary statistics

summary(s)

# plot

plot(metagene_score, col_score, pch = 16, cex = 1.3, col = "blue", 

     main = "Metagene score vs. colitis score", 

     xlab = "Metagene score", ylab = "Colitis score")

# line

abline(lm(formula = 'col_score ~ metagene_score', df))

56) Check version of R

version

57) Academic citation for package

citation("shiny")

58) Online book on using R for forecasting (link)

59) Find out type of an R object

class(imports_quarterly_raw)

60) SQL commands in R (link)

library("sqldf")

# find top exports by commodity

output_df = sqldf("select sum(exports_quarterly_raw.Weight) as sum_weight, exports_quarterly_raw.code2, lookup_code.description from exports_quarterly_raw inner join lookup_code on exports_quarterly_raw.code2 = lookup_code.tariff_code group by exports_quarterly_raw.code2 order by sum_weight desc")

 

61) which command in R to find indices that match in array or data frame

lookup_code[which(lookup_code$tariff_code == unique_codes)]

62) Code to load data from tables, perform SQL queries and plot histograms 

library("sqldf")

library(ggplot2)

df = sqldf("select sum(imports_quarterly_raw.Weight) as sum_weight from imports_quarterly_raw where imports_quarterly_raw.Code2 = 39 and imports_quarterly_raw.Port_of_Discharge in ('Boston) and imports_quarterly_raw.mode_of_transport = 'SEA' group by imports_quarterly_raw.Year, imports_quarterly_raw.Qtr")

df_ts = ts(df, start=c(start_year, start_qtr), frequency=4)

plot.ts(df_ts, main = "Imports: Plastics And Articles Thereof",col = "blue")

# histogram/distribution of plastics

qplot(df, data = df, geom = 'histogram')

63) Save plots in R

pdf(file = "imports_plastics_yearly.pdf")

plot.ts(df_ts, main = "Imports: Plastics And Articles Thereof (aggregated yearly)",col = "blue")

dev.off()

# eps format 

postscript(file="heatmap.eps",horiz=TRUE,

           onefile=FALSE,width=8.5,

           height=11,paper=letter) 

dev.off()

# for ggplots

ggsave(filename = 'temp.pdf',  useDingbats=FALSE)

ggsave(filename = "boxplot_numpeptides_vs_age.eps", gp, device = "eps")#,  useDingbats=FALSE)

64) Names of columns

names()

colnames()

rownames()

65) Other useful commands

unique(), sort(), order()

66) Apply function to to data along axis

apply(), lapply(), tapply()

67) Run UNIX command from R

system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=finished.pdf *.pdf")

68) concatenate strings and return in R

paste()

69) Get first few characters of a string in R

substr(df_desc[1,], 1, 8)

70) Assign names of columns in data frame

colnames(output_df) <- c("Year", "Qtr", "Country_origin")

71) Split based on character (strsplit)

# Say start_date is of format 2015-03-30 and you want to extract year

as.character(start_date)

strsplit(as.character(start_date), "-")

unlist(strsplit(as.character(start_date), "-"))

72) Convert to character 

as.character()

73) unlist

74) Histogram in R

r_all = hist(output_df$sum_fob_amt_aud, main="Histogram for FOB value", 

     xlab="FOB Value (AUD)", 

     border="black", 

     col="blue", breaks = 100) # xbin = c(1,6)

# access all counts, frequencies and midpoints of bins

r_all$counts

r_all$mids

75) Convert weird data format into time series 

# Say date is in format 17/2/09

as.Date(exports_customs_short$intended_export_date, "%d/%m/%y")

# converts to R Date format

# Now can feed into SQL (sqldf) format

df = sqldf("select sum(fob_amt_aud) as sum_fob_amt_aud, intended_export_date from exports_customs_short where transport_mode_type_code = 'S' and loading_port_code = 'AUSYD' and cargo_id_type_code = 'C'  group by intended_export_date order by intended_export_date")

# this will sort by date

76) Using Jupyter with R (link)

Install miniconda

conda install -c r ipython-notebook r-irkernel

ipython notebook

77) More group by dates in R SQLDF   (METHOD: to coerce some field to be date or numeric etc in sqldf)

# Recast datetime columns as R Date object (makes SQL queries easier)

imports_customs_short$search_arrival_date <- (as.Date(imports_customs_short$search_arrival_date, "%d/%m/%y"))

df= sqldf("select sum(gross_weight) as sum_gross_weight, search_arrival_date from imports_customs_short group by search_arrival_date order by search_arrival_date ")

 

78) Dataframes can have elements with different types but ts (time series) object can have only one type: Date/time

79) Create data frame

data.frame()

OR

final_list_peptides_combined = data.frame(peptide = "peptide", stringsAsFactors=FALSE)

# NOTE: stringsAsFactors=FALSE will prevent columns being used as factors (stackoverflow)

# Characters by default become factors or categorical variables

# Examples

# create data frame

df_bic_values = as.data.frame( cbind(as.numeric(mcl.model_1$bic), as.numeric(mcl.model$bic), as.numeric(mcl.model_3$bic)), stringsAsFactors=FALSE )

# rename columns

colnames(df_bic_values) <- c("BIC 1 normal distribution", "BIC 2 normal distributions", "BIC 3 normal distributions")

80) expand.grid()

81) Hash table in R

list(name = var)

82) Generating random numbers 

sample(1:10, 100, replace=TRUE) # 10 samples between 1 and 10

 

83) Create empty data frame (from link)

date.frame() # see also link

84) Append to existing data frame (use rbind)

  temp_df = round(temp_random_norm * df$sum_container_count[iCount] * i_percentage_export) 

  

  # append to data frame

  demand_export_df = rbind(demand_export_df, t(temp_df))

85) Generate range of numbers with step 1

20:80

86) Write data frame to csv file

write.csv( demand_export_df, file = "demand_export_shortcustoms_60days.csv" )

write.csv( demand_export_df, file = "demand_export_shortcustoms_60days.csv" , row.names = FALSE, quote=FALSE) # if no row names and no quote around characters

87) For loop

# temp_code is an array

for (icol in 1:length(temp_code))

{

  icol

}

for (temp in temp_code)

{

  temp

}

88) Select rows from a data frame when a column matches some element in a list (like a SQL select * from where in ....) (link)

# Use %in%

# Select total sum of FOB value for all HS codes that match with selected UN code i.e. sum over all subcommodities

df_imports_sea_sydney_trim = df_imports_sea_sydney[ df_imports_sea_sydney$Code2 %in% list_hs_codes, ]

89) Stacked plots using ggplot2 (from stackoverflow)

# Create dummy data

set.seed(11)

df <- data.frame(a = rlnorm(30), b = 1:10, c = rep(LETTERS[1:3], each = 10))

library(ggplot2)

# geom_area function to fill in area

ggplot(df, aes(x = b, y = a, fill = c)) + geom_area(position = 'stack') # + ylab(str_ylabel) + xlab(str_xlabel) + ggtitle(str_title)

 

ggsave(filename = 'temp.pdf' , useDingbats=FALSE)

ggsave(filename = "boxplot_numpeptides_vs_age.eps", gp, device = "eps")#,  useDingbats=FALSE)

90) List all variables in the workspace

ls()

91) Barplots in R

# X is some matrix with columns being separate bars in the barplot

X <- as.matrix(cbind(df_short_cont_trim_summed$totalsum_fob_amt_aud, df_short_noncont_trim_summed$totalsum_fob_amt_aud))

barplot(X,

              names.arg = c("Containerized","Bulk"), ylab = "")

title(main = "Bar chart of Containerized vs. Bulk", font.main = 4)

92) Data munging (dplyr package)

      

dplyr

93) Concatenate vectors in R

vector_1 = c(1,2)

vector_2 = c(1,2)

c(vector_1, vector_2)

94) R code to test for allometric power law relationship (on bitbucket)

95) R code to perform forecasting and SQL like queries for a road accident forecasting and data exploration project (on bitbucket) (deployed on shinyapps)

     Deployed web application to perform data exploration using SQL-like queries and perform machine learning analysis (on shinyapps)

96) R function for generating stacked plots using ggplot (on bitbucket)

97) R code to test for allometric power law relationship (on bitbucket)

98) feather -   fast on disk format for data frames in R and python

99) Iterate over lists using apply, mapply (link1) (link2 with examples)

100) Call a function on each element of multiple lists

mapply

# can also be done on multiple cores using library('parallel')

# install.packages('parallel')

  list_un_commodity_code = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21)

  list_datasource = c(1)

  list_piechart = c(1)

  list_fob_weight = c(1)

  # generate_save_all_plots is a function that takes 4 parameters (scalars)

  mapply(generate_save_all_plots, list_un_commodity_code, list_datasource, list_piechart, list_fob_weight)

101) Datetime library (link)

library(lubridate)

year(as_datetime("2015-03-03"))

month(as_datetime("2015-03-03"))

lubridate::year(lubridate::today())

lubridate::make_date(  year = lubridate::year(df_epic_labtests_withage$date_of_birth_approx),

                                                                        month = lubridate::month(df_epic_labtests_withage$date_of_birth_approx),

                                                                        day = lubridate::day(df_epic_labtests_withage$date_of_birth_approx)

                                                                      )

102) Fit a gamma distribution to data

library(MASS)

x.gam <- rgamma(100, shape = 1, scale = 1/10)

fitdistr(x.gam, "gamma")

103) String replacement

# gsub, grep, grepl, regexpr 

# replace commas, which cause problem with SQL sum()

gsub(",", "", df_temp$total_cval) 

# replace every character after _ALL with none  (using .*) (courtesy Maria Jose Gomez)

gsub('_ALL.*', '',  temp_file)

104) Create list in R

foo <- 12

bar <- c("a", "b", "e")

newList <- list("integer" = foo, "names" = bar)

105) Function to check NaN, NA, Inf, - Inf  (link)

is.finite()

106) assert() function for testing 

install.packages('testit')

library(testit)

assert("one equals one", 1 == 1)

assert("one equals one", 1 == 2)

# check/assert if numeric using is.finite()

assert("check is numeric \n", is.finite(variable) )

# generic function for asset checking

assert_generic_function <- function(var_input, i_ignore)

{

    ###############################################################################

    # Testing function

    ###############################################################################

  

    # var_input: data frame etc

    # i_ignore: if TRUE then ignore and if FALSE then assert

  

    # Assert statement

    if (i_ignore == FALSE)

    {

        assert("check is numeric \n", is.finite(var_input) )      

    }

}

107) List of very handy packages in R (list)

has lubridate, sqldf, plyr, stringr, qcc (for quality control)

108) Package that is good replacement for ggplot

ggraptR (link)

109) R package for rapid data mining and data exploration (like Weka) (courtesy Yuriy Tyshetskiy)

rattle

110) Dynamic reports in R (like ipython notebook)

knitr (link)

making multi-chapter reports with references in knitr (R markdown) (link)

R markdown notebooks and shiny R (link)

chunk options

with 

```r { }

```

    include = FALSE prevents code and results from appearing in the finished file. R Markdown still runs the code in the chunk, and the results can be used by other chunks.

    echo = FALSE prevents code, but not the results from appearing in the finished file. This is a useful way to embed figures.

    message = FALSE prevents messages that are generated by code from appearing in the finished file.

    warning = FALSE prevents warnings that are generated by code from appearing in the finished.

    fig.cap = "..." adds a caption to graphical results.

# Insert value of variable dynamically in report

`r toString(nreplicates)`

# table of values using kable

df_bic_values = as.data.frame( cbind(as.numeric(mcl.model_1$bic), as.numeric(mcl.model$bic), as.numeric(mcl.model_3$bic)), stringsAsFactors=FALSE )

colnames(df_bic_values) <- c("BIC 1 normal distribution", "BIC 2 normal distributions", "BIC 3 normal distributions")

knitr::kable( df_bic_values, caption = "Model selection."  )

# Multi-panel plots in knitr (link)

```{r fig_sicseq_deseq2_metacluster2_3, fig.cap="Identification of genes",   echo=FALSE}

library(png)

library(grid)

library(gridExtra)

img1 <-  rasterGrob(as.raster(readPNG("figures/deseq2_inflam_cluster_scseq_metacluster2.png")), interpolate = FALSE)

img2 <-  rasterGrob(as.raster(readPNG("figures/deseq2_cluster_scseq_metacluster3.png")), interpolate = FALSE)

grid.arrange(img1, img2, ncol = 2)

```

# run knitr and rmarkdown from command line

R

library(knitr)

library(rmarkdown)

knitr::knit('~/periphery_project/combined_analysis_v4.rmd')

link

rmarkdown::render("test.Rmd", "html_document")

Rscript -e 'library(rmarkdown); rmarkdown::render("/path/to/test.Rmd", "html_document")'

# A typical header of a R markdown file will look like

---

title: "Analysis and Writeup"

header-includes:

- \usepackage{placeins}

- \usepackage{float}

- \floatplacement{figure}{H}

output:

  pdf_document:

    fig_caption: yes

    keep_tex: yes

    latex_engine: xelatex

    number_sections: yes

  word_document: default

  html_document:

    df_print: paged

bibliography: Periphery_project.bib

urlcolor: blue

---

```{r include=FALSE}

knitr::opts_chunk$set(echo = TRUE)

knitr::opts_chunk$set(cache = TRUE)

knitr::opts_chunk$set(warning = FALSE)

knitr::opts_chunk$set(out.extra = '')

#knitr::opts_chunk$set(fig.pos = 'H')

```

\begin{centering}

\vspace{3 cm}

\Large

\normalsize

Soumya Banerjee, `r format(Sys.time(), "%b %d %Y")`

\vspace{3 cm}

\end{centering}

\setcounter{tocdepth}{2}

\tableofcontents

\newpage

```{r,include=FALSE}

library(knitr)

library(gridExtra)

library(rmarkdown)

# EQUATIONS in rmarkdown

$$ eGFR = eGFR_{0} + b_{before}*t_{before} $$

```

Italics in rmarkdown using *metafor*

Code can be rendered or shown in rmarkdown using

```

dsBaseClient::ds.summary(x='surv_object')

```

111) Bayesian structural time series in R

CausalImpact

112) C++ coding best practices and tips for using with R on Mac OS X (link)

113) Reading a csv file

file_strong_binders = read.csv('strong_binders.txt', sep = ' ', header = FALSE, stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

OR

df_monocyte_genes_IGP <- read.csv('IGP_SuppTable_MOD.csv', 

                                sep = ',', header = TRUE, 

                                stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

# quote="" # for dealing with ' etc in files

NOTE: stringsAsFactors=FALSE will prevent columns being used as factors (stackoverflow)

NOTE: quote="" # for dealing with ' etc in files

114) Find number of rows in a dataframe

nrow()

115) Debugging in R

function like matlab keyboard command is browser()

see link1 and link2

116) Loop through all files in directory (link)

files <- list.files(path=str_results_directory, pattern="*.csv", recursive=FALSE)

117) Adding elements to a list in a while loop (stackoverflow)

list_val = list()

i_temp_counter = 1

for (temp_lambda in c(0.01, 0.1, 0.15,  0.2, 0.3,  0.5, 0.6))

{

  list_val[[i_temp_counter]] = temp_lambda

  

  i_temp_counter = i_temp_counter + 1

  

}

# get all elements in a nice vector

c(as.numeric (list_val) )

118) Sample from an empirical distribution (non-parametric distribution)

sample(empirical_distr, N_samples, replace=TRUE)

119) Concatenate columns of a dataframe and create one column with concatenated data

# Uses with() and paste0()

copy_numbers_AIREWT_file$V6 = with(copy_numbers_AIREWT_file,                             paste0(copy_numbers_AIREWT_file$V1,copy_numbers_AIREWT_file$V2,copy_numbers_AIREWT_file$V3,copy_numbers_AIREWT_file$V4,copy_numbers_AIREWT_file$V5))

120) Replace non-zero elements of a matrix (stackoverflow)

# find all non-zero elements

i_index_notzero = which(all_genes_withzeros_LOG10 != 0, arr.ind = TRUE)

# replace non-zero elements with log10

all_genes_withzeros_LOG10[i_index_notzero] = log10(all_genes_withzeros_LOG10[i_index_notzero])

OR if needs to be replaced with a single value

all_genes_withzeros_LOG10[all_genes_withzeros_LOG10==0] = NA

OR get index of column and remove that column

idx_to_remove <- which(d$PC2 > 0.55)[1]

IMPORTANT: how to remove column (courtesy Maria Jose Gomez)

# remove that column with [,c(-254)] notation

df_file_str_filename_RNASeq_withpath_NORM_PCAREMOVE = df_file_str_filename_RNASeq_withpath_NORM[,c(-idx_to_remove)]

# use grep to column position and then remove it

i_col_to_remove = grep("Comment_IBD_in_family", colnames(df_filename_scseq_infl_cluster))

df_filename_scseq_infl_cluster = df_filename_scseq_infl_cluster[,c(-i_col_to_remove)]

IMPORTANT: how to remove column (courtesy Stephen Sansom)

rownames(all_gene_matched_withprobeid_osm_agg) <- all_gene_matched_withprobeid_osm_agg$ENSEMBL

all_gene_matched_withprobeid_osm_agg$ENSEMBL <- NULL

121) Turn off scientific notation in R (no 1e-3 only 0.003)

options(scipen = 100)

122) append to a file and add a data frame

# write.table with append=TRUE

# will not work with write.csv

write.table(df_final_metadata, file=filename_metadata_FINAL,

                      row.names = FALSE, quote=FALSE, append = FALSE, sep = ",")  #, col.names = NA)

123) continue statement in R

next

124) create directory in R

dir.create("./data")

125) Download files from an URL and save it in a directory (courtesy Elsa Arcaute)

data_url <- “XX”

dir.create("./data")

dest_file="data/temp.zip"

downloader::download(data_url,destfile=dest_file, mode = "wb")

126) Unzip file

unzip(dest_file,exdir="data")

127) Report generation using knitr with caching of results

128) Updating R, R Studio and updating packages (link)

# run from R console

install.packages('devtools') #assuming it is not already installed

library(devtools)

install_github('andreacirilloac/updateR')

library(updateR)

updateR(admin_password = 'Admin user password')

# the update R Studio if required

# update.packages() from R console

129) Heatmaps in R (using the ComplexHeatmap package)

library(ComplexHeatmap)

library(circlize)

set.seed(123)

# mat is generic matrix

mat = cbind(rbind(matrix(rnorm(16, -1), 4), matrix(rnorm(32, 1), 8)),

            rbind(matrix(rnorm(24, 1), 4), matrix(rnorm(48, -1), 8)))

# permute the rows and columns

mat = mat[sample(nrow(mat), nrow(mat)), sample(ncol(mat), ncol(mat))]

rownames(mat) = paste0("R", 1:12)

colnames(mat) = paste0("C", 1:10)

Heatmap(mat)

# Also heatmaps using the heatmap.2 package (link) 

# Also machine learning chapter in Irizarry ( book)

# Excerpt here

library(RColorBrewer)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

Now, pick the genes with the top variance over all samples:

library(genefilter)

rv <- rowVars(e)

idx <- order(-rv)[1:40]

While a heatmap function is included in R, we recommend the heatmap.2 function from the gplots

package on CRAN because it is a bit more customized. For example, it stretches to fill the window.

Here we add colors to indicate the tissue on the top:

library(gplots) ##Available from CRAN

cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(tissue)]

head(cbind(colnames(e),cols))

heatmap.2(e[idx,], labCol=tissue,

trace="none",

ColSideColors=cols,

col=hmcol)

# save heatmap.2

png(filename = "heatmap.png")

heatmap.2(m)

dev.off()

130) Get column and row sums

rowSums()

colSums()

131) ggplot with labels

# TRICK CONCEPT: add a label to the dataframe with condition (stim)

library(ggplot2)

# ggplot magic to plot with labels

gp <-  ggplot(d$stim, aes(d$PC1, d$PC2, color=stim)) + geom_point(size=5)

print(gp)

132) How to append a string to the end of each entry of a column (using paste0)

file_all_target_mappings_orig$SAMPLE <- paste0(file_all_target_mappings_orig$SAMPLE,".CEL")

133) rownames() and colnames() to get names of rows and columns

134) METHOD: add a column in the dataframe for stimulation condition

d$stim = sqldf("select file_all_target_mappings_orig_LPS_LPSaIL10.STIMULATION as stim 

               from file_all_target_mappings_orig_LPS_LPSaIL10 

               inner join d 

               on file_all_target_mappings_orig_LPS_LPSaIL10.SAMPLE = d.sample")

135) sorting and ordering an array on a column (using   order() )

order(x[,1]) # sort on first column

idx = order(x[,1]) # get indices of sorted array

x[idx,] # feed these indices into array to get sorted

# all together (succinct)

x[order(x[,1]), ] 

136) Summary statistics of an object

summary(as.vector(mat))

137) Filtering out rows of an array

idx2 <- apply(rawdata, 1, max) > 10 # only take in rows that have max > 10

rawdata <- rawdata[idx2, ]   # use this to index on the array

138) table() and is.na() to summarize data

# how many are NA ?

is.na( select(x=hta20transcriptcluster.db, keys=c(rownames(all_genes_pca_loading_sorted)), columns = c("GENENAME")   ) )

table( is.na( select(x=hta20transcriptcluster.db, keys=c(rownames(all_genes_pca_loading_sorted)), columns = c("GENENAME")   ) )  )

139) METHOD: turn column names into a field

# add a column for sample name

# d is a dataframe

d$sample = rownames(d)

140) Data wrangling using dplyr (cheatsheet)

141) R markdown notebooks and shiny R (link)

chunk options

with

```r { }

```

    include = FALSE prevents code and results from appearing in the finished file. R Markdown still runs the code in the chunk, and the results can be used by other chunks.

    echo = FALSE prevents code, but not the results from appearing in the finished file. This is a useful way to embed figures.

    message = FALSE prevents messages that are generated by code from appearing in the finished file.

    warning = FALSE prevents warnings that are generated by code from appearing in the finished.

    fig.cap = "..." adds a caption to graphical results.

MORE LINKS (link)

142) NICE COMMANDS

subset(), with()

143) qplot() in ggplot2

library(ggplot2)

library(datasets)

data(mpg)

qplot(displ, hwy, data=mpg)

144) Chi-squared test with contigency table (from link) (also see link)

M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))

dimnames(M) <- list(gender = c("F", "M"),

                    party = c("Democrat","Independent", "Republican"))

View(M)

(Xsq <- chisq.test(M))  # Prints test summary

Xsq$observed   # observed counts (same as M)

Xsq$expected   # expected counts under the null

Xsq$residuals  # Pearson residuals

Xsq$stdres     # standardized residuals

Fisher's exact test  (link)

fisher.test(M)

145) Convert data frame column (severe) to a factor with meaningful labels and then plot in ggplot  TRICK

df_metagene_risk_score_final_joined_atNF_combined$severe = factor(df_metagene_risk_score_final_joined_atNF_combined$severe, c(0,1), c("Severe", "Not Severe"))

ggplot(df_metagene_risk_score_final_joined_atNF_combined, aes(x=metagene_score, colour = severe)) + geom_density()

146) ggplot multi-panel figure (courtesy Maria Jose Gomez Vazquez)

library(gridExtra)

gp1 <- qplot(metagene_score, data = df_metagene_risk_score_final_joined_atNF_combined, geom = 'histogram')

gp2 <- qplot(metagene_score, data = df_metagene_risk_score_final_joined_atNF_combined, geom = 'density')

grid.arrange(gp1, gp2, nrow=1)

147)  R max DLL reached problem

.Renviron file in home directory

R_MAX_NUM_DLLS=250

148) ggplot theme

  theme_set(theme_gray())

149) select annotations for microarrays

# http://bioconductor.org/packages/3.6/data/annotation/manuals/illuminaMousev2.db/man/illuminaMousev2.db.pdf

source("https://bioconductor.org/biocLite.R")

BiocInstaller::biocLite("illuminaMousev2.db")

library(illuminaMousev2.db)

select(x=illuminaMousev2.db, c('ILMN_1217075'), c("ENSEMBL","SYMBOL"))

# find all columns using 

AnnotationDbi::colnames(hta20transcriptcluster())

150) transpose a data frame

#transpose data matrix

matrix_all_genes_risk_scseq_array_weighted = (t(as.matrix(all_genes_risk_scseq_array_weighted[,1:321])))

# make new data frame

df_matrix_all_genes_risk_scseq_array_weighted = as.data.frame(matrix_all_genes_risk_scseq_array_weighted,

                                                              stringsAsFactors = FALSE)

# get new column names and rownames

colnames(df_matrix_all_genes_risk_scseq_array_weighted) <- all_genes_risk_scseq_array_weighted$gene_name

df_matrix_all_genes_risk_scseq_array_weighted$id = rownames(matrix_all_genes_risk_scseq_array_weighted)

151) Get indices of columns that have only a particular pattern (courtesy Maria Jose Gomez)

# get indices of columns that have only BAMC id

temp_idx_columns_with_risk = grep("BAMC", colnames(all_genes_risk_scseq_array_weighted))

152) Wilcoxon test (link)

boxplot(Ozone ~ Month, data = airquality)

wilcox.test(Ozone ~ Month, data = airquality,

            subset = Month %in% c(5, 8))

153) Use of %in% command to subset a data frame

# all IDs that are in another data frame

df_matrix_all_genes_risk_scseq_array_weighted_CONTROL =

df_matrix_all_genes_risk_scseq_array_weighted[df_matrix_all_genes_risk_scseq_array_weighted$id %in% df_risk_healthy$SampleID,]

# logical operation

# all IDs if in another data frame, set a new column to the logical result (TRUE/FALSE)

df_matrix_all_genes_risk_scseq_array_weighted$control = df_matrix_all_genes_risk_scseq_array_weighted$id %in% df_matrix_all_genes_risk_scseq_array_weighted_CONTROL$id

154) ggplot boxplot (link)

gp <- ggplot(df_num_metagene_score_NEW, aes(x=control,y=num_metagene_score_NEW, fill=control))

gp <- gp + geom_boxplot()

gp

155) reshape, cast and melt (link)

156) ggplot facet_wrap facet wrap (link)

      gp <- ggplot(df_combined_score_data_frame, aes(x=control,y=norm_metagene,fill=control))

      gp <- gp + geom_boxplot()

      gp <- gp + facet_wrap( ~ signature, ncol=2, scales = "free_y")

      gp <- gp + xlab("IBD vs. Control")

      gp <- gp + ylab("Normalised metagene")

      gp

# where control, norm_metagene and signature are columns of a data frame named df_combined_score_data_frame

157) Demonstration of how to load a data frame, melt it and plot it in ggplot using facet_wrap (courtesy Kevin Rue-Albrecht and Maria Jose Gomez)

t = 

  cbind(

    sum(df_il23$d0_R1),

    sum(df_il23$d0_R2),

    sum(df_il23$d0_R3),

    sum(df_il23$d0_R4)

  )

    

t = rbind(t,  cbind(

  sum(df_il23$d1_R1),

  sum(df_il23$d1_R2),

  sum(df_il23$d1_R3),

  sum(df_il23$d1_R4)

)

)

t = rbind(t,  cbind(

  sum(df_il23$d3_R1),

  sum(df_il23$d3_R2),

  sum(df_il23$d3_R3),

  sum(df_il23$d3_R4)

)

)

t = rbind(t,  cbind(

  sum(df_il23$d6_R1),

  sum(df_il23$d6_R2),

  sum(df_il23$d6_R3),

  sum(df_il23$d6_R4)

)

)

t = rbind(t,  cbind(

  sum(df_il23$d14_R1),

  sum(df_il23$d14_R2),

  sum(df_il23$d14_R3),

  sum(df_il23$d14_R4)

)

)

t = rbind(t,  cbind(

  sum(df_il23$d21_R1),

  sum(df_il23$d21_R2),

  sum(df_il23$d21_R3),

  sum(df_il23$d21_R4)

)

)

df_tt = data.frame(t, stringsAsFactors = FALSE)

df_tt$day = c("0","1","3","6","14","21")

df_tt$source = c("metacluster5","metacluster5","metacluster5","metacluster5","metacluster5","metacluster5")

df_tt_melted = melt(df_tt, id=c("day","source"))

# convert to numeric

df_tt_melted$day = as.numeric(df_tt_melted$day)

# make factors      

df_tt_melted$source = factor(df_tt_melted$source)

df_tt_melted$day  = factor(df_tt_melted$day)

gp <- ggplot(df_tt_melted, aes(x=day,y=value,fill=source))

gp <- gp + geom_point()

gp <- gp + xlab("Days") + ylab("Metagene score")

gp

theme_set(theme_gray())

gp <- ggplot(df_tt_melted, aes(x=day,y=value,color=source)) # fill=source

gp <- gp + geom_point()

gp <- gp + facet_wrap(~ source, ncol=2, scales = "free_y")

gp <- gp + xlab("Days") + ylab("Metagene score")

gp

158) Make a list unique OR rownames unique (courtesy Maria Jose Gomez)

mat_mod$gene_name <- make.unique(mat_mod$gene_name)

159) Assigning a Boolean value of an expression to a column (using %in%)

df_matrix_all_genes_risk_scseq_array_weighted$control = df_matrix_all_genes_risk_scseq_array_weighted$ccfaid %in% df_matrix_all_genes_risk_scseq_array_weighted_CONTROL$ccfaid

160) Other useful functions

merge, unique, subset, intersect, Reduce

intersect with multiple sets (link)

Reduce( intersect, list(a,b,c) )

161) ggplot point with line group by source (link) (courtesy Stephen Sansom)

      theme_set(theme_gray())

      gp <- ggplot(df_combined_mouse_metagene_withsource, aes(x=day,y=value,group=source)  )

      gp <- gp + geom_point(aes(color=source))

      gp <- gp + geom_line(aes(color=source))

      gp

162) Heatmap with ggplot (link)

geom_tile()

163) zscore 

scale()

164) remove a row due to a row having all zeros or all NAs (link) (stackoverflow)

# use complete.cases()

zscored_values <- zscored_values[(complete.cases(zscored_values)),]

# remove a row due to a column having NA or Inf

# relies on the fact that Inf * 0 is NA and NA * 0 is NA  (link)

complete.cases(c(1, Inf))

[1] TRUE TRUE

complete.cases(c(1, Inf*0))

[1]  TRUE FALSE

complete.cases(c(1, Inf*0, NA))

[1]  TRUE FALSE FALSE

complete.cases(c(1, Inf*0, NA*0))

[1]  TRUE FALSE FALSE

# final code to remove a row due to a column having NA or Inf

dt_tests_with_age <- dt_tests_with_age[complete.cases(dt_tests_with_age$egfr_recalc * 0), ]

# ALTERNATIVE remove NA or remove a row with NA

na.omit(d2)

165) cut() command to subset a dataframe

166) Rename one column of a data frame

df_combined_score_data_frame[(df_combined_score_data_frame$signature == "metacluster5"),]$signature = "ribosomal"

167) Select those rows that have a particular characteristic and subset an array based on that

#keep those IDs that have medicine usage in active disease

# select IDs

temp_idx_to_keep = which( colnames(m_to_be_plotted_withbar_subset_active) %in% df_subset_subset$id )

# subset array

m_to_be_plotted_withbar_subset_active_SUBSET = m_to_be_plotted_withbar_subset_active[,temp_idx_to_keep]

168) ggplot histogram (courtesy Maria Jose Gomez Vazquez) (other examples of ggplot histogram from R for Data Science)

library(gridExtra)

gp1 <- qplot(metagene_score, data = df_metagene_risk_score_final_joined_atNF_combined, geom = 'histogram')

gp2 <- qplot(metagene_score, data = df_metagene_risk_score_final_joined_atNF_combined, geom = 'density')

grid.arrange(gp1, gp2, nrow=1)

OR (courtesy Shanquan Chen)

NOTE: use of position="identity". CAUTION: default is position="stacked"  (link)

gp <- ggplot(df_temp_combined_ckd_grouped, aes(x=min_age_calculated_asof_test/365, fill=flag_responder_binary) )

gp <- gp + geom_histogram(bins=15, position = "identity", alpha=0.4) # stat = "count"

gp <- gp + xlab('Age of subject at time of Lithium cessation')

gp <- gp + ylab('Frequency')

gp

OR

# facet wrap and facet grid histogram

# where data frame df_late_inflamed_fc1_new_distribution_metagene_melted has score in field value and column named histology

theme_set(theme_gray())

gp2 <- ggplot(df_late_inflamed_fc1_new_distribution_metagene_melted, aes(x=value,fill=histology))

gp2 <- gp2 + geom_histogram()#(alpha=0.5)

gp2 <- gp2 + xlab("Metagene score") + ylab("Number of patients")

#gp <- gp + facet_wrap(histology~.)

gp2 <- gp2 + facet_grid(histology~.)

gp2

# facet wrap and facet grid

# geom_density with alpha makes it transparent

theme_set(theme_gray())

gp <- ggplot(df_late_inflamed_fc1_new_distribution_metagene_melted, aes(x=value,fill=histology))

gp <- gp + geom_density(alpha=0.5)

#gp <- gp + facet_wrap(histology~.)

gp <- gp + facet_grid(histology~.)

gp

print(gp)

# if you want to change the order in which variables are facetted TRICK (courtesy Kevin-Rue Albrecht)

df_combined_metagene_with_hist_melted_refactored = melt(df_combined_metagene_with_hist, 

                                                        id=c("id","histology"))

df_combined_metagene_with_hist_melted_refactored$histology = factor(df_combined_metagene_with_hist_melted_refactored$histology,

                                                                    c("macro", "micro", "control")

)

theme_set(theme_gray())

gp2 <- ggplot(df_combined_metagene_with_hist_melted_refactored, aes(x=value,fill=histology))

gp2 <- gp2 + geom_histogram()#(alpha=0.5)

gp2 <- gp2 + xlab("Metagene score") + ylab("Number of patients")

#gp2 <- gp2 + facet_wrap(histology~.)

gp2 <- gp2 + facet_grid(histology~.)

gp2

# if you want to relabel factors (e.g. if 1 = "General", 2 = "Academic")

df_counts_data$prog <- factor(df_counts_data$prog,

                              levels = c(1,2,3),

                              labels = c("General", "Academic", "Vocational")

                              ) 

169) Named list in R (like key value dict in Python)

mouse_late_early = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")

170) Run R from command line

nohup R --no-save < analysis_degeneracy_ageing.R

171) Colour schemes and palettes for R (link)

172) Get rid of factors in a data frame (courtesy Maria Jose Gomez)

# CONCEPT: factors are stored as numbers. They need to be converted to character first and then numbers

# the field active has become a factor; convert it into numeric

df_data$active = as.numeric(as.character(df_data$active))

173) remove rows in a matrix or dataframe with duplicated gene names (courtesy Kevin Rue-Albrecht)

# using duplicated function

all_genes = all_genes[!duplicated(all_genes$gene_name), ]

 

174) Excellent ggplot visualization website (link)

         Another excellent ggplot visualization link (link)

175) dplyr ordering (courtesy Maria Jose Gomez)

y_comb <- y_comb %>% group_by(source) %>% arrange(source, desc(V1)) %>% as.data.frame()

176) geom_path ggplot with facet wrap

# use dplyr to sort data frame helps with geom_path

y_comb <- y_comb %>% group_by(source) %>% arrange(source, desc(V1)) %>% as.data.frame()

#############################

# Plot facet

#############################

gp <- ggplot(y_comb, aes(x=y_comb$V1, y=y_comb$V2))

gp <- gp + geom_path()

gp <- gp + facet_wrap(~source, ncol=2)

gp <- gp + geom_hline(data=y_comb, aes(yintercept = threshold), linetype="dashed")

gp <- gp + xlim(0,1)

gp <- gp + ylim(0,1)

gp <- gp + xlab("Recall")

gp <- gp + ylab("Precision")

gp

177) Violin plot with jitter and median lines shown

theme_set(theme_classic())

gp2 <- ggplot(df_combined_metagene_with_hist_melted_refactored, aes(y=value))

gp2 <- gp2 + geom_violin(aes(x = disease_category, y=value, fill=disease_category), draw_quantiles = c(0.5))

gp2 <- gp2 + geom_jitter(height = 0, width = 0.1, aes(x=disease_category, y=value))

gp2 <- gp2 + xlab("Disease category") + ylab("Metagene score")

# gp2 <- gp2 + xlim(20,100)

# gp2 <- gp2 + ylim(10,90)

#gp <- gp + facet_wrap(histology~.)

# gp2 <- gp2 + facet_grid(histology~.)

gp2

178) Volcano plot using ggplot

  #######################################

  # prettier MA plot using ggplot

  #######################################

  

  res_df$gene_name <- rownames(res_df)

  data = res_df

  data$color = "default"

  

  # Invert sign of log fold change

  data$log2FoldChange = -1 * data$log2FoldChange

  

  data$log2_adjpvalue = -log2(data$padj)

  

  # Horizontal lines for fold change

  i_y_line_threshold_intercept_top = -log2(0.05)

  i_y_line_threshold_intercept_bottom = -1

  i_x_line_threshold_intercept_1 = log2(2)

  i_x_line_threshold_intercept_2 = log2(1/2)

  

  top = head(data, n=20) # by -pvalue 

  data$color[abs(data$log2FoldChange)> 1 & data$padj<0.05] <- "> 2 fold, padj < 0.05"

  

  cols <- c("default"="darkgrey", "> 2 fold, padj < 0.05"="red")

  

  # theme_set(theme_gray())

  theme_set(theme_classic())

   # gp <- ggplot(data, aes(baseMean, log2FoldChange, color=color)) + geom_point(alpha=0.8)

  gp <- ggplot(data, aes(log2FoldChange, log2_adjpvalue, color=color)) + geom_point(alpha=0.8)

  # gp <- gp + scale_x_continuous(trans="log2")

  #gp <- gp + geom_text(data=top, aes(label=gene_name), vjust=1,hjust=-0.2, color="black")

  gp <- gp + geom_text_repel(data=top, aes(label=gene_name), vjust=1, hjust=-0.2, color="black")

  gp <- gp + geom_hline(yintercept=i_y_line_threshold_intercept_top, linetype="dashed")

  gp <- gp + geom_vline(xintercept=i_x_line_threshold_intercept_1, linetype="dashed")

  gp <- gp + geom_vline(xintercept=i_x_line_threshold_intercept_2, linetype="dashed")

  #gp <- gp + geom_hline(yintercept=i_y_line_threshold_intercept_bottom, linetype="dashed")

  gp <- gp + scale_color_manual(values=cols)

  gp <- gp + ggtitle("DESeq2 analysis of monocytes across all donors within a metacluster")

  # gp <- gp + xlab("Fold change (log2)expression (baseMean)") + ylab("Fold change (log2)")

  gp <- gp + xlab("Log fold change") + ylab("-log (p-value)")

  print(gp)

  

  str_maplot_filename = paste0("volcanoplot_meta", i_meta_cluster_id, ".pdf")

  ggsave(filename = str_maplot_filename, gp, device = "pdf",  useDingbats=FALSE)

  ggsave(filename = "boxplot_numpeptides_vs_age.eps", gp, device = "eps")#,  useDingbats=FALSE)

  

179) Replace ALL missing values in a data frame with NA

idx_missing_values_in_cells = which(df_joined_cleaned_removed == '', arr.ind = TRUE)

df_joined_cleaned_removed[idx_missing_values_in_cells] = 'NA'

Find missing values (link)

which(!complete.cases(df_politeness))

180) grep to select those rows that do not have a pattern (courtesy Maria Jose Gomez)

yourdata[-grep(c('/'),yourdata$DOB),]

OR

yourdata[grep(c('/'),yourdata$DOB, invert=TRUE),]

181) geom_error_bar() tutorial (link1) (link2)

# example code using geom_errorbar() and geom_errorbarh()

theme_set(theme_classic())

theme_set(theme_gray())

i_y_line_threshold_intercept_atnf   = 0.3661417

i_y_line_threshold_intercept_active = 0.7118644

gp <- ggplot(df_summary_aupr_signatures_RE, aes(x=aTNF,y=active_disease, color=signature))

#gp <- ggplot(df_summary_aupr_signatures_RE, aes(x=aTNF,y=active_disease, shape=signature, color=signature))

gp <- gp + geom_point(aes(size=number_genes)) # df_summary_aupr_signatures_RE, size=number_genes

gp <- gp + xlab("AUPR aTNF")

gp <- gp + ylab("AUPR active disease")

gp <- gp + geom_vline(xintercept=i_y_line_threshold_intercept_atnf, linetype="dashed")

gp <- gp + geom_hline(yintercept=i_y_line_threshold_intercept_active, linetype="dashed")

gp <- gp + geom_text_repel(data=df_summary_aupr_signatures_RE, 

                           aes(label=signature), vjust=1,hjust=-0.2, min.segment.length = unit(1, "mm"),

                           color="black") #(data=show, aes(label=gene_name))

gp <- gp + geom_errorbar(data = df_summary_aupr_signatures_RE_2, aes(ymin=active_disease-1.96*sd_active/i_num_samples, ymax=active_disease+1.96*sd_active/i_num_samples), width=0.1  )

gp <- gp + geom_errorbarh(data = df_summary_aupr_signatures_RE_2, aes(xmin=aTNF-1.96*sd_atnf/i_num_samples, xmax=aTNF+1.96*sd_atnf/i_num_samples))#, width=0.1  )

# gp <- gp + xlim(0,1)

# gp <- gp + ylim(0.6,1)

gp <- gp + scale_size(range = c(2,10))

#gp <- gp + scale_x_continuous()

#gp <- gp + scale_y_continuous()

gp

ggsave(filename = "signature_genes.pdf", useDingbats=FALSE)

182) Converting weird date formats (adapted from stackoverflow and r-bloggers)

# date is Aug-08

x <- df_filename_scseq_infl_cluster$Last_colonoscopy[1]

x <- as.Date(paste("01-", x, sep = ""), format = "%d-%b-%y")

x <- format(x, "%d/%m/%Y")

df_filename_scseq_infl_cluster$Last_colonoscopy[1] <- as.character(x)

183) Summarize and aggregate data using the aggregate() command

#collapse expression values for genes with multiple probes

data$probe_id <- NULL

data$gene_symbol <- NULL

data.agg <- aggregate(.~gene_id, data, median)

rownames(data.agg) <- data.agg$gene_id

data.agg$gene_id <- NULL

184) Use of apply() to filter data for QC

# Expression filter

hist(data.matrix)

data.matrix <- data.matrix[apply(data.matrix,1,max)>4,]

print(dim(data.matrix))

# Variance filter

hist(apply(data.matrix,1,var), breaks=30)

data.matrix <- data.matrix[apply(data.matrix,1,var)>0.5,]

print(dim(data.matrix))

185) Barplot in ggplot

theme_set(theme_classic())

#theme_set(theme_gray())

gp <- ggplot(data=df_intercept, aes(x=gene_name, y=coeff))

gp <- gp + geom_bar(stat = "identity")

# gp <- gp + geom_bar()

# gp <- gp + facet_wrap(~gene_name)

gp <- gp + xlab("Gene name")

gp <- gp + ylab("Coefficient")

gp <- gp + coord_flip()

gp

186) parsing command line arguments (optparse)

library(optparse)

# Options ----

option_list <- list(

    make_option(c("--seuratobject"), default="begin.Robj",

                help="A seurat object after PCA"),

    make_option(c("--clusterids"), default="none",

                help="A list object containing the cluster identities"),

    make_option(c("--outdir"), default="seurat.out.dir",

                help="outdir")

    )

opt <- parse_args(OptionParser(option_list=option_list))

outPrefix <- file.path(opt$outdir,"markers.summary")

s <- readRDS(opt$seuratobject)

saveRDS(df_combined_all_features, file = "source_data/df_combined_all_features.rds")

df_combined_all_features <- readRDS(file = "source_data/df_combined_all_features.rds")

187) RColorBrewer diverging palettes (link)

188) Negation of %in% operator (from link) 

'%!in%' <- function(x,y)!('%in%'(x,y))

which(df_filename_scseq_infl_cluster$misc %!in% c('Y','N','NA'))

189) TRICK: Use as.numeric() will convert non-numeric to NA (link)

190) Check if a file exists

file.exists(str_file_name)

191) ggplot linear regression (link) (link)

        theme_set(theme_gray())

        #theme_set(theme_classic())

        gp <- ggplot(data = df_colscore_metagenescore_noDUPL, aes(x=metagene_score, y=col_score))

        gp <- gp + geom_point()

        gp <- gp + geom_smooth(method = "lm", level=0.95) #se=TRUE)

        # gp <- gp + geom_smooth(span=0.3)

        gp

        OR

        

ggplot(data1, aes(x=QUET, y=SBP))+

    geom_point()+

    geom_smooth(method=lm, se=TRUE)

192) PCA in R

pca <- prcomp(2^emat, scale=T)

#pca_probelevel <- prcomp(2^emat_probelevel, scale=T)

# summary of data

summary(as.vector(emat))

#summary(as.vector(emat_probelevel))

dim(emat)

#dim(emat_probelevel)

# PCA analysis on raw data

#pca = prcomp(log10(geneCore@assayData$exprs + 1))  

#pca = prcomp(log10(geneCore@assayData$exprs), scale = TRUE) 

#pca = prcomp(geneCore@assayData$exprs)#, scale. = TRUE)  

# pca = prcomp(2^geneCore@assayData$exprs, scale = TRUE)

dim(2^geneCore@assayData$exprs)

# get rotation

d <- as.data.frame(pca$rotation)

#d_probelevel <- as.data.frame(pca_probelevel$rotation)

# PCA plot (not pretty)

plot(d$PC1,d$PC2)

#plot(d_probelevel$PC1,d_probelevel$PC2)

# generate scree plot

(summary(pca))$importance["Proportion of Variance",]

head(pca$sdev)

plot( (summary(pca))$importance["Proportion of Variance",], 

      xlab="component",

      ylab="proportion of variance",

      main="scree plot")

################################################

# what are top loading genes

################################################

head(pca$x[order(pca$x[,1]),])

head(pca$x[order(pca$x[,1], decreasing = TRUE ),] , n =20 )

#head(pca_probelevel$x[order(pca_probelevel$x[,1]),])

# order(pca$x[,1])

# idx = order(pca$x[,1]) # get indices of sorted array

# pca$x[idx,] # feed these indices into array to get sorted

# head(pca$x[idx,]) # head of that

#head(pca_probelevel$x[order(pca_probelevel$x[,1], decreasing = TRUE ),] , n =20 )

#AnnotationDbi::select(x = hta20transcriptcluster.db, keys = c("TC04001102.hg.1"), columns = c("GENENAME"))

#AnnotationDbi::select(x = hta20transcriptcluster.db, keys = c("TC12001162.hg.1"), columns = columns(x=hta20transcriptcluster.db))

#top_pca_loaded_genes <- head(pca$x[order(pca$x[,1], decreasing = TRUE ),] , n = 100 )

#AnnotationDbi::select(x = hugene10sttranscriptcluster.db, keys = c(rownames(top_pca_loaded_genes)), columns = c("GENENAME"))

# now look at all genes PCA loading sorted

#all_genes_pca_loading_sorted <- pca$x[ order(pca$x[,1], decreasing = TRUE )  ,]

# what are these genes

#AnnotationDbi::select(x=hugene10sttranscriptcluster.db, keys=c(rownames(all_genes_pca_loading_sorted)), columns = c("GENENAME")   )

###############################################################################################

# Prettier plots using ggplot

###############################################################################################

# add a column for sample name OR stimulation condition etc (NOTE: d is a dataframe)

#d$sample = rownames(d)

# add a column in the dataframe for stimulation condition

# join by column name with original metadata

#d$stim = sqldf("select file_all_target_mappings_orig_LPS_LPSaIL10.STIMULATION as stim 

#               from file_all_target_mappings_orig_LPS_LPSaIL10 

#               inner join d 

#               on file_all_target_mappings_orig_LPS_LPSaIL10.SAMPLE = d.sample")

# ggplot magic to plot with labels

#gp <-  ggplot(d$stim, aes(d$PC1, d$PC2, color=stim)) + geom_point(size=5)

#print(gp)

193) Return multiple arguments from function

function <- test()

{

# code here

return(list(

          demographics = demographics,

          diagnosis = diagnosis

        ))

}

# call function

list[demographics, diagnosis] <- test()

194) Charts, different kinds of plots and network diagrams using the networkd3 package in R (link)

195) Data munging like filter using dplyr and magrittr (courtesy Shanquan Chen)

# table command displays summary statistics (like frequency)

# filter filters the data

library(sqldf)

library(magrittr)

library(dplyr)

library(lubridate)

library(stringr)

mdat <- all_patients %$% table( Diagnosis_Code )

mdat <- data.frame(mdat)

mdat %>% filter(Diagnosis_Code == "F20")

sqldf::sqldf("select * from mdat where Diagnosis_Code like 'F20' ")

195) (CONTINUED FROM ABOVE) Pattern matching using stringr  (courtesy Shanquan Chen)

mdat$Diagnosis_Code

mdat[stringr::str_detect(mdat$Diagnosis_Code, 'F20.*'),]

196) (CONTINUED FROM ABOVE) Data summarization using dplyr and lubridate (courtesy Shanquan Chen)

table command displays summary statistics (like frequency)

year in lubridate gets the year from a date

table(lubridate::year(mdat$Diagnosis_Last_date))

197) (CONTINUED FROM ABOVE)  Convert data to long format using dcast (courtesy Shanquan Chen) (link)

df_meds_summarized <- dcast(df_records_patients_using_topMeds, rid + Death_Flag + Gender_Code ~ drug, length)

Also 

dcast(drug_table, rid ~ drug, any)

 

Also

dcast(drug_table, rid ~ drug, fun.aggregate = function(d) {ifelse(any(d), 1, 0)})

Also see 

drug_table <- data.table(

    rid = c(rep(c('a', 'b', 'c'), each=3), 'a'),

   drug = c('citalopram', 'insulin', 'simvastatin',

            'clozapine', 'metformin', 'simvastatin',

            'fluoxetine', 'mirtazapine', 'omeprazole', 'simvastatin')

)

# add a column value

drug_table$value <- TRUE

library(reshape2)

patient_table <- dcast(drug_table, rid ~ drug, any)

 patient_table <- dcast(drug_table, rid ~ drug, fun.aggregate = function(d) {ifelse(any(d), 1, 0)})

 patient_table <- dcast(drug_table, rid ~ drug, fun.aggregate = function(d) {ifelse(length(d) > 0, 1, 0)}, value.var = 'drug')

198) (CONTINUED FROM ABOVE)  Alternative using data.table

# The problem:

# - to take a list of drugs, and map them to a list of patients,

answering the

#   question: for each patient, does the drug occur?

library(data.table)

drug_table <- data.table(

     rid = c(rep(c('a', 'b', 'c'), each=3), 'a'),

     drug = c('citalopram', 'insulin', 'simvastatin',

              'clozapine', 'metformin', 'simvastatin',

              'fluoxetine', 'mirtazapine', 'omeprazole', 'simvastatin')

)

patient_table <- data.table(

     rid = c('a', 'b', 'c'),

     dummy = c(1, 2, 3)

)

setkey(patient_table, rid)

drug_frequencies <- drug_table[, .(n = .N), by = drug]

top_drugs <- drug_frequencies[n >= 1, drug]  # silly

for (drugname in top_drugs) {

     cat("Thinking about drug:", drugname, "\n")

     count_this_drug_by_patient <- drug_table[drug == drugname, .(n =

.N), by = rid]

     patient_table[

         count_this_drug_by_patient, on = "rid",

         # c(drugname) := n  # number of occurences of the drug

         c(drugname) := ifelse(n >= 1, 1, 0)  # does the drug occur?

     ]

     drugval <- patient_table[, get(drugname)]

     patient_table[is.na(drugval), (drugname) := 0]  # force NA to zero

}

# Note re data.table:

# ... trailing [] to prevent "doesn't print first time" bug:

#

 # https://stackoverflow.com/questions/32988099/data-table-objects-not-printed-after-returned-from-function

# https://github.com/Rdatatable/data.table/blob/master/NEWS.md#bug-fixes-5

patient_table <- patient_table[]  # fix nonprinting bug; see above

print(patient_table)

199) Connect to databases using the RODBC package (link)

200) All kinds of classification trees, regression trees, decision trees and random forests in R (link)

201) Making indicator variables

# I() function makes an indicator variable here the response 0 or 1

# For example

gam_m_fit_simple_logistic <- gam::gam(formula = I(wage>300) ~ s(year,df=4) + s(age,df=4),

                                      family = binomial,

                                      data = Wage)

202) Boxplot in R with formula

boxplot(frequency ~ attitude*gender, data=df_politeness, col=c("white","lightgray") )

203) ggplot axis label tilt to help with visibility (link)

# Use theme()

gp2 <- ggplot(df_patients_joined_meds_topdrugs, aes(x=drug))#,fill=Gender_Code))

gp2 <- gp2 + geom_histogram(stat = "count")#(alpha=0.5)

gp2 <- gp2 + xlab("Drug type") + ylab("Count") + title("Breakdown of data by type and name")

gp2 <- gp2 + facet_wrap(Gender_Code~.)

gp2 <- gp2 + theme(axis.text.x=element_text(angle = -90, hjust = 0))

# gp2 <- gp2 + facet_grid(drug~.)

gp2

204) R script for installing multiple packages (link)

205) Evaluate using eval()  (link)

eval(parse(text = '5+5'))

# ANOTHER TWIST if you are calling it in a function and want this to persist in the environment

#  that called the function

eval( parse(text = '5+5'), envir = parent.frame() )

# CONCEPT: parent.frame()  gets the environment in the call stack above this call

see code here

206) Other dplyr functions

merge, mutate, select, arrange

207) tidy function in broom package to clean up or tidy data

library(broom)

tidy_lm_aq <- tidy(lm_aq)

tidy_lm_aq

208) The R package visdat can also be used for rapid data visualization (see link)

library(visdat)

vis_dat(airquality)

209) Use of apply() and paste() to concatenate or paste all columns of a data frame;

how to concatenate or paste all columns of a data frame (on stackoverflow)

apply(df_tcr_repertoire_age_1wk, 1, paste, collapse="")

210) Namespace like system to manage R functions (link)

For an example see (link)

        util = new.env()

        util$bgrep = function [...]

        util$timeit = function [...]

        while("util" %in% search())

              detach("util")

        attach("util")

211) Power calculation in R (link)

sample size

power.t.test

power.t.test(n = NULL, power = .95, sd = 5, alternative = "two.sided", sig.level = 0.001, delta = 0.1)

212) Create venn diagram and euler diagram (link) (link)

VennDiagram

eulerr

library(VennDiagram)

grid.newpage()

draw.triple.venn(area1 = 22, area2 = 20, area3 = 13, n12 = 11, n23 = 4, n13 = 5,

    n123 = 1, category = c("Dog People", "Cat People", "Lizard People"), lty = "blank",

    fill = c("skyblue", "pink1", "mediumorchid"))

213) Great list of functions for doing statistics and plotting (rlib)

# Load packages and settings

library(sqldf)

library(ggplot2)

library(knitr)

library(rmarkdown)

library(gplots)

library(RColorBrewer)

library(reshape2)

library(png)

library(grid)

library(gridExtra)

library(lme4)

library(lmerTest)

library(lubridate)

library(reshape2)

library(data.table)

library(pheatmap)

library(cba)

library(rstanarm)

library(shinystan)

library(rstantools) # make a rstanarm package

library(bayesplot)

library("factoextra") #install_github("kassambara/factoextra")

library(rpart)

REQUIRED_PACKAGES <- c(

  "car",  # may need: sudo apt install libgfortran3

  "data.table",

  "ggplot2",

  "lubridate",

  "RODBC",

  "shinystan"

  # "survival"

  # "timereg"

)

LIBRARY_PREFIX <- "https://egret.psychol.cam.ac.uk/rlib/"

#source(paste0(LIBRARY_PREFIX, "datetimefunc.R"))

source("datetimefunc.R")

source(paste0(LIBRARY_PREFIX, "dbfunc.R"))

source(paste0(LIBRARY_PREFIX, "listfunc.R"))

source(paste0(LIBRARY_PREFIX, "listassign.R"))

source(paste0(LIBRARY_PREFIX, "miscfile.R"))

source(paste0(LIBRARY_PREFIX, "misclang.R"))

source(paste0(LIBRARY_PREFIX, "miscstat.R"))

source(paste0(LIBRARY_PREFIX, "stanfunc.R"))

source(paste0(LIBRARY_PREFIX, "cris_common.R"))

misclang$library_or_install(REQUIRED_PACKAGES)

# As advised by Stan:

rstan::rstan_options(auto_write = TRUE)

options(mc.cores = parallel::detectCores())

# Other misc commands

BIGSEP <-   "===============================================================================\n"

SMALLSEP <- "-------------------------------------------------------------------------------\n"

214) ggplot with text within figure

          require(ggplot2)

          require(PRROC)

          require(grid)

          #################################################################################

          # convert to ggplot

          theme_set(theme_classic())

          y <- as.data.frame(prroc_object$curve)

          gp <- ggplot(y, aes(y$V1, y$V2))

          gp <- gp + geom_path()

          gp <- gp + xlim(0,1)

          gp <- gp + ylim(0,1)

          gp <- gp + geom_hline(yintercept=i_y_line_threshold_signif, linetype="dashed")

          gp <- gp + xlab("Recall")

          gp <- gp + ylab("Precision")

         

         

          ############

          # add text

          my_text <- paste("AUC = ", as.character(round(prroc_object$auc.integral, digits = 2)

                                                  ))

          my_grob = grid.text(my_text,  gp=gpar(col="black", fontsize=14, fontface="bold"))

          gp <- gp + annotation_custom(my_grob)

         

         

          gp

215) Create combinations using combn or expand.grid

combn(c('a','b','c'), 2) # all combinations of 3 letters taken 2 at a time

combn(letters[1:4], 2)

expand.grid( height = seq(60, 80, 5), weight = seq(100, 300, 50), sex = c("Male","Female") )

216) Create a new data.table

sample_data <- data.table(subject_id = 1:100)

OR

dt_demographics <- data.table(

  read.csv('demo.csv',

                            sep = ',', header = TRUE,

                            stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

)

setkey(dt_demographics, STUDY_SUBJECT_DIGEST)

OR

Create a new data.table and add columns

library(data.table)

library(lubridate)

library(tidyr)

START_YEAR <- 2016

END_YEAR <- 2017

# Define years of interest

years <- data.table(

    year = START_YEAR:END_YEAR

)

# add a new column

years[, year_start_date := as.Date(paste0(year, "-01-01"))]

years[, year_end_date := as.Date(paste0(year, "-12-31"))]

MALE <- "M"

FEMALE <- "F"

# create another data.table

sexes <- data.table(sex = c(MALE, FEMALE))

age_ranges <- data.table(

    age_range_start = c(0, 1, seq(5, 90, 5))

)

age_ranges[, age_range_end := age_range_start + 4]

 

OR

library(dplyr)

library(data.table)

df

   %>%

  mutate(date = as.Date(date)

   %>%

   as.data.table()

217) filter in data.table

# Filter to relevant years only:

subject_mortality_f20 <- subject_mortality_f20[

    year >= lubridate::year(subject_start_date) &

    year <= lubridate::year(subject_end_date)

]

218) Merge in data.table

# Merge in standardized mortality (create a new table, for clarity):

subject_and_population_mortality_f20 <- merge(

    subject_mortality_f20,

    population_mortality[, c("year", "sex", "age_range", "deaths_per_person_per_year")],

    by.x = c("year", "sex", "age_range_this_year"),

    by.y = c("year", "sex", "age_range")

)

219) Crossing or cross-product or all combinations in data.table

# Make a table with a row for each subject/year combination:

subject_mortality_f20 <- data.table(tidyr::crossing(sample_data_f20, years))

220) For each patient get minimum date and group by each patient (select min and group by in SQL) se using data.table (link)

# use .SD[1] construct

setkeyv(dt_tests_with_age, c("STUDY_SUBJECT_DIGEST", "ResultDate"))

# sort by age

dt_first_test_per_patient <- dt_tests_with_age[, .SD[1], by=c("STUDY_SUBJECT_DIGEST")] # *** check sorting

# .SD[1] gets first row

# set that to a new column

dt_first_test_per_patient$first_creatinine_datetime <- dt_first_test_per_patient$ResultDate

# then merge with original table

dt_egfr_with_age_withfirst <- data.table::merge.data.table(x = dt_tests_with_age,

                                                                                             y = dt_first_test_per_patient,

                                                                                             by = 'patid')

221) Create a summary of number of counts of a field (like select sum in SQL) (link)

# Create sum of tests for each patient using .N

dt_patients_with_two_creatinines <- dt_all_tests[

    ,

    .(n_creatinine_tests = .N),

    by = 'STUDY_SUBJECT_DIGEST'  # by subject

]

# using .N

222) Saving memory by using pass by reference (:=) instead of pass by value (<-) in data.table

# 1st alternative in plain R using <- 

dt_tests_with_age$age_calculated_asof_test <- as.numeric( dt_tests_with_age$ResultDate - dt_tests_with_age$date_of_birth_approx )  

# 2nd alternative in data.table using :=

dt_tests_with_age[

  ,

  age_calculated_asof_test := as.numeric(ResultDate, date_of_birth_approx)

]

223) Do inner join in data.table using merge  (other joins like left join link)

# set primary key in data.table using setkey

setkey(dt_tests_with_age,          STUDY_SUBJECT_DIGEST)

setkey(dt_first_test_per_patient, STUDY_SUBJECT_DIGEST)

dt_tests_with_age <- merge(dt_tests_with_age, dt_first_test_per_patient, by=c("STUDY_SUBJECT_DIGEST")) 

224) Rename columns in data.table

setnames(dt_tests_with_age

         ,c("ResultNumeric.x", "age_calculated_asof_test.x", "GENDER_DESC.x" , "ETHNIC_GROUP_DESC.x")

         ,c("ResultNumeric"  , "age_calculated_asof_test"  , "GENDER_DESC"   , "ETHNIC_GROUP_DESC"  )

         )

225) Get counts (like count(distinct(*)) in SQL) using uniqueN in data.table

dt_demographics <- data.table(

  read.csv('demo.csv',

                            sep = ',', header = TRUE,

                            stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)

)

setkey(dt_demographics, STUDY_SUBJECT_DIGEST)

n_patients <- uniqueN(dt_demographics$STUDY_SUBJECT_DIGEST)

# get counts and group by year of birth and gender (courtesy Rudolf Cardinal)

# like count(distinct()) ….. group by yob, gender

dt_merged_treatment[, .(n_lithium_patients = .N), by = .(yob, gender)] 

# can create a new column based on this using := 

dt_merged_treatment[,  n_lithium_patients := .(n_lithium_patients = .N), by = .(yob, gender)] 

# if any errors use the conflicted package

library(conflicted)

conflict_prefer("list", "gsubfn")

# now select only those with 2 lithium_patients i.e select count(n_lithium_patients) as n  …..

dt_selected <- dt_merged_treatment[ n_lithium_patients >= 2 ]

226) Union of two tables (like union in SQL) using rbindlist (link)

list_tables  <- list(dt_epic_labtests, dt_meditech_labtests)

dt_all_tests <- data.table::rbindlist(list_tables, use.names = TRUE, fill = TRUE)

225) Read excel in R (readxl)

# Read in real ONS data

library(readxl)

smr_table_ons <- readxl::read_xls(path = 'deathsummarytables2017final.xls', sheet=4)

226) Make or turn a matrix into a vector in row-major format (link)

vector_ons_mortality_2016_2017 <- as.vector(t(mat_ons_mortality_2016_2017))

227) ifelse command

df_time_to_event$time_to_event <- ifelse(df_time_to_event$Death_Flag,

                                         as.numeric(df_time_to_event$Date_Of_Death_lubridate - df_time_to_event$Date_Of_referral_lubridate)/DAYS_PER_YEAR,

                                         as.numeric(END_DATE_CONSTANT_lubridate - df_time_to_event$Date_Of_referral_lubridate)/DAYS_PER_YEAR

                                        )

228) Memory free up size in R (link) (other tips to free up memory)

memory.limit()

memory.limit(size = 4095) # in MB

gc()  # garbage collector

229) Randomly sample lines from a big csv file (link)

# randomly sample a few patients

set.seed(12)

i_NUM_PATIENTS_SAMPLE = 10000

temp_list_row_numbers = 1:nrow(dt_demographics)

# now sample some patients from them

i_temp_sample_patient_index = sample(temp_list_row_numbers, i_NUM_PATIENTS_SAMPLE)

# CAUTION: overwrite demographics table

dt_demographics <- dt_demographics[i_temp_sample_patient_index, ]

230) Dropfactors or droplevels (link)

df_temp_summary <- droplevels(df_temp_summary)

231) Remove duplicate columns (link)

df_time_to_event <- subset(df_time_to_event, select=which(!duplicated(names(df_time_to_event))))

232) Get the last 3 digits of a string (link) courtesy Shanquan Chen

# use stringr

x <- "some text in a string"

substrRight <- function(x, n){

  substr(x, nchar(x)-n+1, nchar(x))

} Convert string to lower case str_to_lower() courtesy Rudolf Cardinal 233) Two strategies to merge, join, inner join (courtesy Rudolf Cardinal) dt_case_control_matched <- data.table::merge.data.table(dt_merged_treatment[i_temp_index,],

                                                            dt_merged_control,

                                                            #allow.cartesian = TRUE,

                                                            by=c("yob","gender","region")

                                                            )

# more efficient

dt_case_control_matched <- dt_merged_control[yob == dt_merged_treatment[i_temp_index,]$yob

                                                  &

                                                gender == dt_merged_treatment[i_temp_index]$gender

                                                  &

                                                region == dt_merged_treatment[i_temp_index]$region 

                                                ]

234) Very advanced data.table (code) (stackoverflow)

Set column value using the first next row in the same group that meets a condition

# Adapted from

# https://stackoverflow.com/questions/55070842/set-column-value-using-the-first-next-row-in-the-same-group-that-meets-a-conditi?noredirect=1&lq=1

dt <- structure(list(

  id        = c(1L, 1L, 2L, 2L, 3L, 4L, 5L, 5L, 5L, 6L, 6L, 6L),

  code      = c("p", "f", "f", "p", "p", "<NA>", "f", "p", "p", "f", "p", "p"),

  date_down = structure(c(17897, 17898, 17898, 17899, 17900, 17901, 17903, 17903, 17905, 17906, 17906, 17906), class = "Date"),

  date_up   = structure(c(17898, 17899, 17898, NA, NA, 17901, 17904, 17904, 17905, 17906, 17906, 17907), class = "Date")),

  class     = c("data.table", "data.frame"),

  row.names = c(NA, -12L))

data.table::setDT(dt)  # to reinit the internal self ref pointer (known issue)

dt[, row := .I] # add row numbers (as in all the solutions)

data.table::setkey(dt, id, row, date_down)  #set key for dt

head(dt)

# TODO: do for each patient patid == temp_patient_id

# For all rows of dt, create found_date by reference :=

dt[, found_date :=

     # dt[code=='p'][dt,

     dt[code=='p'][.(.SD),   # our subset (or another data.table), joined to .SD (referring to original dt)

                   .( x.date_up ),

                   on = .(id==id, row > row, date_up > date_down), 

                   mult = "first"] ]

head(dt)

View(dt)

235) Remove rows in data.table

# remove those rows marked for removal

dt_egfr_with_age_merged_lithium <- dt_egfr_with_age_merged_lithium[ remove_row == 0,  ]

235) Leave one out cross validation in R using the loocv package (link)

236) Create an empty data frame initialise an empty data frame of a certain size (from stackoverflow)

mat_empty_matrix <- matrix(data = 0, nrow = 100, ncol = 21  )

df_empty <-  data.frame( mat_empty_matrix, stringsAsFactors = FALSE)

235) Convert a data frame to a matrix using data.matrix

236) Show fit of linear model lm() using ggplot

abline

link

237) Testing using testthat (VERY GOOD link) (link) (link)

# very simple example below

library(testthat)

library(devtools)

string <- "Testing is fun! hellooooo..."

expect_match(string, "Testing") 

expect_match(string, 'hello', ignore.case = TRUE)

# expect_output(as.character(string), 'hello')

238) Making a package in R using devtools (link)

devtools::test()

devtools::test(filter = 'test.R') # where test.R is below

library(testthat)

string <- "Testing is fun! hellooooo..."

expect_match(string, "Testing") 

expect_match(string, 'hello', ignore.case = TRUE)

239) Check where libraries are stored

.libPaths()

240) Read in and write Stata (.dta) files

library(haven)

library(tidyverse)

a <- haven::read_dta(file = 'test.dta')

haven::write_dta(data = b, path = '~/mod.dta')

241) R environments and calling scope and lexical scope (link)

242) Initialize a matrix with zeros (link)

matrix(data = 0, nrow = 10, ncol = 10)

243) Bookdown in R

Simple example (link)

244) Write the output of a model summary or any model to a log file

str_log_filename = 'log_model1_offals.csv'

sink(str_log_filename)

print(summary(meta_model))

 sink()

245) Append to list incrementally

list_input_covariate = c('BMI')

list_input_covariate = cbind(list_input_covariate, 'AGE')

list_input_covariate = cbind(list_input_covariate, 'GENDER')

OR

list_input_covariate = NULL

list_input_covariate = cbind(list_input_covariate, 'AGE')

list_input_covariate = cbind(list_input_covariate, 'GENDER')

246) How to make a package in R (link, github, youtube)

247) Generate synthetic data in R

simstudy package

https://cran.r-project.org/web/packages/simstudy/vignettes/simstudy.html library(simstudy) # generate definition def <- defData(varname = "age", dist = "normal", formula = 10, variance = 2) def <- defData(def, varname = "female", dist = "binary", formula = "-2 + age * 0.1", link = "logit") def <- defData(def, varname = "visits", dist = "poisson", formula = "1.5 - 0.2 * age + 0.5 * female", link = "log") # generate data dd <- genData(1000, def) dd

248)Plotting CDF cumulative distribution function using stat_ecdf() in ggplot

link

link

# age should be a factor

gp <- ggplot(data = df_all, x = energy, aes(x = energy, color = age_factor)

gp <- gp + geom_step(aes(y=..y..), stat = "ecdf")

gp


249) R profiler (link, link) (code)

Rprof

profviz


250) qgraph for plotting graphs in R (link)


251) bookdown

Recipe by Tom Bishop

Bookdown recipe:


  In RStudio, install the bookdown package

install.packages('bookdown')

  Start a new project, but choose the project to be a book.

  In Github create an empty project  (no README.md) with the same name as your book

  In RStudio go to Tools -> Project options -> Version control -> Git (it will ask you to reload the project)

On Git command line, add the remote and push: git remote add origin https://github.com/USERNAME/PROJECT.git git push -u origin master


  In the file _bookdown.yml, and the line output_dir: "docs"

  In the file _output.yml, delete the parts that refer to generating a pdf (I couldn’t get this to work):

bookdown::pdf_book: includes: in_header: preamble.tex latex_engine: xelatex citation_package: natbib keep_tex: yes bookdown::epub_book: default


  Check in local changes

In Github go to Settings -> Pages and set the source to the main branch and folder to docs



252) sessionInfo()

to get all details related to session, what version packages are installed, version, etc.


253) Display pretty formatted R code in rmarkdown


```{r eval=FALSE}

x = y

```


254) Object oriented programming class in R S3 and S4 class (link)


255) Create inst/CITATION file for citing R package (link)

library(usethis)

usethis::use_citation()


Page updated
Google Sites
Report abuse
