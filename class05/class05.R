#Class 05 R graphics intro

# My first boxplot
x <- rnorm(1000,0)
boxplot(x)

summary(x)
hist(x)

boxplot(x, horizontal = TRUE)


weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE, sep = "") 
weight$Age
weight$Weight

plot(weight$Age, weight$Weight, type = "o", pch = 15, cex = 1.5, lwd = 2, ylim = c(2,10), xlab = "Age (months)", ylab= "Weight (kg)", main= "Weight by Age")


#using barplot
barplot(VADeaths, beside = FALSE)


#input our feature count data
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE, sep = "\t")
barplot(mouse$Count, horiz = TRUE, names.arg=mouse$Feature, main = "Some title", las = 2)

#change margin so that we can see labels
par(mar = c(3.1, 11.1, 4.1, 2))
barplot(mouse$Count, names.arg=mouse$Feature, horiz = TRUE, ylab = "", main = "Number of features in the mouse GRCm38 genome", las = 1, xlim =c(0,80000))

# $ is used as [,column]

#add some color
barplot(mouse$Count, horiz = TRUE, names.arg = mouse$Feature, las =2, col = rainbow(12))

#histogram
hist(rnorm(10000), rnorm(10000)+4, breaks = 50)
mf <- read.table("bimm143_05_rstats/male_female_counts.txt", header = TRUE, sep = "\t")
#barplot always needs height - use the Count column for this
barplot(mf$Count, names.arg = mf$Sample, horiz = FALSE,  las = 2, col=c("red", "blue"))



#Up Down Expression data

e<- read.table("bimm143_05_rstats/up_down_expression.txt", header = TRUE, sep = "")
#how many genes
nrow(e)
#how many up, down and changing?
table(e$State)
plot(e$Condition1, e$Condition2, col = e$State)

#plot
palette(c("red", "lightgray", "blue"))
plot(e$Condition1, e$Condition2, col = e$State)






