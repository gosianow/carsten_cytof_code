
library("cytofCore")

cytofCore.updatePanel()





library("citrus")

citrus.launchUI()

############################################

# Example with a single experimental condition (unstimulated data)
# Where the data lives
dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")

# Create list of files to be analyzed
fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))

# List disease group of each sample
labels = factor(rep(c("Healthy","Diseased"),each=10))

# List of columns to be used for clustering
clusteringColumns = c("Red","Blue")

# Run citrusAnalysis
results = citrus.full(fileList,labels,clusteringColumns,dataDirectory)

# Should be used to plot results
plot(results, outputDirectory="")

############################################

# Example using multiple experimental conditions, including a
# comparison of stimulated to unstimulated damples

# Where the data lives
dataDirectory = file.path(system.file(package = "citrus"),"extdata","example2")

# Create list of files to be analyzed
fileList = data.frame(unstim=list.files(dataDirectory,pattern="unstim"),stim1=list.files(dataDirectory,pattern="stim1"))

# Disease group of each sample
labels = factor(rep(c("Healthy","Diseased"),each=10))

# Vector of parameters to be used for clustering
clusteringColumns = c("LineageMarker1","LineageMarker2")

# Vector of parameters to calculate medians for
functionalColumns = c("FunctionalMarker1","FunctionalMarker2")

# Form condition comparison matrix to tell Citrus which conditions to analyze compare
conditionComparaMatrix = matrix(T,nrow=2,ncol=2,dimnames=list(colnames(fileList),colnames(fileList)))
conditionComparaMatrix[2]=F

# Run citrusAnalysis
results = citrus.full(fileList,labels,clusteringColumns,dataDirectory,
  conditionComparaMatrix=conditionComparaMatrix,
  featureType="medians",
  medianColumns=functionalColumns)

# Result should be plotted
# plot(results,outputDirectory="/path/to/outputDirectory")








