# Step-by-step Guide to Perform Differential Expression Analysis using DESeq2
## Teodora Tockovska
## April 9, 2021

This tutorial is the first document out of 3 that will go into detail about performing transcriptome analyses using real data. Within this tutorial, you will be introduced to the R library `DESeq2`, which is a method for differential analysis of genes by using count data. Love et al. (2014) describe that the analysis of count data from RNA-seq is a fundamental task to understand the systematic changes across experimental conditions.  

In our case, the data that we will be using are real sample data from ducks. The birds were infected with ABBV-1, aquatic bornavirus. This virus is known to infect waterfowl and the purpose of the experiment was to understand how this virus affects Muscovy ducks (*Cairina moschata*) when they were inoculated through the intracranial routes. Ducks were divided in control and intracranial groups. The intracranial groups were infected with virus whereas the control groups were inoculated with control inoculums in the brain. At weeks 4 and 12, the birds from infected and control groups were euthanized and brain samples were harvested for RNA extractions. There were a total of 28 samples that were submitted for RNA sequencing with 7 intracranial inoculated week 4 samples, 7 control week 4 samples, 7 intracranial inoculated week 12 samples, and 7 control week 12 samples.  

Sample data was created by retaining 30% of the duck metadata for the purposes of this tutorial. Hence, there are 16 total samples. The breakdown of groups is as follows: 4 intracranial inoculated week 4 samples, 4 control week 4 samples, 4 intracranial inoculated week 12 samples, and 4 control week 12 samples per species. Because we have control and infected birds, the goal is to find differentially expressed genes (DEGs) and annotating those genes to understand how the virus affected the hosts at weeks 4 and 12. Please refer to the table of contents below to view the general breakdown of the tutorial.

## Table of Contents
* [Data Processing](#loading-libraries-importing-data-and-proceeding-with-data-processing)
    + Preparing data for `DESeq2`.
* Differential Expression Analysis and Exploratory Data Analysis on Entire Datasets
    + Learn the basics of `DESeq2` by performing differential expression analysis (DEA) and exploratory data analysis using the entire [duck](#analyses-on-the-entire-dataset) metadata.
    + Learn how to conduct DEA based on sampling times
    + Create visualizations using built-in functions to understand the data
* Discover Differentially Expressed Genes for duck data split up into their respective sampling time points.
    + [Week 4](#ducks-week-4-dea) Samples
    + [Week 12](#ducks-week-12-dea) Samples
* [Next Steps](#next-steps)
    + Explanation and breakdown of next steps after finding DEGs in samples
* [References](#references)

\newpage

# Loading Libraries, Importing Data, and Proceeding with Data Processing

In this section, you will learn how to load the required libraries, load in your sample datasets, and conduct data processing. Data processing refers to preparing your data to get it ready for analysis. In your case, you would be preparing the data for differential expression analysis so you must prepare it in a specific way so ensure that your analysis will be correct.

First, load in the 5 R libraries. `Tidyverse` is a large package that contains many R libraries for data analysis, data structure maintenance and handling, and plotting. `RColorBrewer` is used to select different colours for visualizations. `vsn` is required for `DESeq2` to normalize results.
```{r, echo=T, results='hide', message=F, warning=F}
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("vsn")
# library("hexbin")
```

The next step is to load in the datasets. The duck sample information is loaded and saved into the variable `duck_samples`. Likewise, the duck gene count data is loaded and saved into a variable called `duck_counts`. The function `read.csv()` is used to read in the CSV files and saved as dataframes. The file is separated by commas (,) which is selected as an argument in the function to separate the data into columns.
```{r}
duck_samples <- read.csv("../EDA_StringTie_RNAseq/sample_metadata_Duck_5Apr2021.csv", sep = ",")
duck_counts <- read.csv("../EDA_StringTie_RNAseq/sample_count_data_Duck_5Apr2021.csv", sep = ",")
```

To understand the data, print the duck_samples dataframe. There are 3 columns with titles "Sample", "Treatment", and "Time". There are 8 control duck samples which are labeled by "D8XX", where *XX* represent numbers. Similarly, there are 8 infected duck samples that are labeled by "D5XX". The sampling time periods are outlined in the "Time" column. For example, Duck D801 was a control duck which was sampled from week 12.
```{r}
duck_samples
```

The gene count data is required for differential expression analysis. There are 7785 rows of genes and there are 18 columns. The first 2 columns are "Gene", which are the Ensembl gene identifier (ID), and "Symbol", which represents the gene symbol or Ensembl ID if there is no gene symbol associated with the Ensembl ID. In total, there are 16 duck individuals. 
```{r}
dim(duck_counts)
colnames(duck_counts)
```

It is important to view the data to have a better understanding of what the dataframe holds. The first 5 rows and columns are displayed. The first 2 columns contain the gene identifiers as described above. Columns 3-5 represent ducks D501, D503, and D506. Each row contains the RNA-seq gene information. For example, duck D501 has 585 counts for gene *CLIPS* (Ensembl ID: ENSAPLG00000010136).
```{r}
duck_counts[1:5, 1:8]
```

The data has been correctly imported. The next step is to extract the columns to be the same as the row names as the count data. This is a required step prior to begin the differential expression analysis. Convert the "Sample" column values to row names using the function `column_to_rownames()`. The piping method (`%>%`) to keep the code organized and clean. Let's view the updated dataframe. As you can see, there are now 2 columns and the row names have been set.
```{r}
duck_samples <- duck_samples %>%
  column_to_rownames(var="Sample")
duck_samples[1:8,]
```

We need to do the same thing on `duck_counts` dataframe which represents the gene count data. Once again, the function `column_to_rownames()` is used to convert the column "Gene" values into row names. The function `select()` selects all columns except for "Symbol" to keep. Let's view the updated dataframe. It is in the format that is needed.
```{r}
duck_counts <- duck_counts %>% 
  column_to_rownames(var="Gene") %>% 
  select(-Symbol)
duck_counts[1:5,1:8]
```

Once the dataframes have been modified, check that all sample IDs in the duck sample data (`duck_samples`) are also in the count data (`duck_counts`) and match their orders. To do this, we must check of the sames of the duck samples are in the count matrix.
```{r}
all(rownames(duck_samples) %in% colnames(duck_counts))
```

The following statement allows us to rearrange the order of the duck samples in our count matrix to match the row sample IDs in our dataframe. This is important because it is required that the dataframe must have the columns of the count matrix in the same order. It will be further explained in the [differential expression analysis](#differential-expression-analysis-dea) section. Check that the orders are correct. 
```{r}
duck_counts <- duck_counts[, rownames(duck_samples)]
all(rownames(duck_samples) == colnames(duck_counts))
```

Lastly, check that the number of columns equals to the number of rows between the two datasets.
```{r}
ncol(duck_counts) == nrow(duck_samples)
```

Because both statements return *TRUE*, this means that the count matrix and sample dataframe are both formatted correctly.

The last thing to do is convert the two columns of interest ("Treatment" and "Time") to factor data objects. Factor objects are used to categorize data and store them as levels. Factors are easier to work with than character (string) or numeric (integer) variables. 

Firstly, convert the values in "Treatment" column by using the function `as.factor()` on the entire column. Because I wanted to make the "Time" values to work with downstream in the analysis, I used the function `gsub()` to replace the "Week" term with an empty character/string value such that only the numbers are left. After doing that, I converted the values to factor by using `as.factor()`. Lastly, check your work and ascertain that the columns are in the correct format. I used the function `map()` which transforms input by applying functions to each column and it returns a vector of information. After the checks, view the updated dataframe.
```{r}
duck_samples$Treatment <- as.factor(duck_samples$Treatment)
duck_samples$Time <- gsub("Week ", "", duck_samples$Time) %>% as.factor()

# Checking factor levels
duck_samples %>% map(~.x %>% levels)

# Checking column value types.
duck_samples %>% map(~.x %>% class) %>% unique() %>% unlist()

# Displaying the first 8 rows
duck_samples[1:8,]
```

The data is now prepared! 

# Analyses on the Entire Duck Dataset

The data has been prepared according to the `DESeq2` manual **Add reference to the link below!** and we can now proceed with the differential expression analysis. 

### Differential Expression Analysis (DEA)

I will be showing you the general steps to perform DEA. First, we will be designing the models which means setting the parameters prior to conducting the `DESeq2` functions to predict which genes could be differentially expressed. To begin, recognize the data types that you are working with because that will dictate how the data can be used. Recall the type of data that we have: gene count data and sample information data. 

Since we recognized the types of data that we have, we will first create a `DESeqDataSet` (DDS) object from the gene count matrix ("duck_counts") and sample information data ("duck_samples") using the function `DESeqDataSetFromMatrix()`. To properly use this function, you will need to provide the gene count matrix and the informations on the samples as a `data.frame` or `DataFrame` object. As it was explained in the previous section, the sample information dataframe must have the columns of the count matrix in the same order, which we proved they are. If the data wasn't properly formatted, then the downstream analysis would result in errors and inaccurate/false results. The data must be formatted correctly prior to this step. 

The third parameter `design` for the function indicates how to model the samples. In our case, we want to measure the effects of the treatment types (control or infected) of the ducks, not based on the time of sampling.
```{r}
dds <- DESeqDataSetFromMatrix(countData = duck_counts,
                              colData = duck_samples, 
                              design = ~ Treatment)
```

Once the DDS object is created, you have the option to pre-filter low reads. I chose to pre-filter low reads to increase the speed of the analyses and to reduce the memory size of the DDS object. Depending on the size of your data, the DDS object could be quite large which results in inconveniences. I followed the guidelines from the `DESeq2` manual to remove reads (gene counts) < 10.  
```{r}
keep_reads <- rowSums(counts(dds)) >= 10
dds <- dds[keep_reads,]
```

Finally, re-level the factors in order to keep the control samples as the reference condition. This step should not be ignored because the code tells the `DESeq2` functions which treatment level to compare against. Recall that there are 2 levels in the treatment group: "Control" and "Infected" ducks. In order to compare the gene counts of the infected ducks against the control ducks, you need to specify the reference level. Ultimately, there are 2 ways to set the reference to the control group by using the function `relevel()` or setting the comparison in the *results* object (shown later). In this case, I demonstrate how to specifiy the level by using the function `relevel()`. It is recommended to set the reference level prior to building the models.  
```{r}
levels(dds$Treatment)
dds$Treatment <- relevel(dds$Treatment, ref = "Control")
```

Now that the DDS object formatted correctly and the reference treatment as been established, you can proceed with building the models. The function `DESeq()` contains all the necessary differential expression analysis steps where the gene expression estimations are calculated statistically. For more information on this function, look at the manual page `?Deseq`. The model predictions can take a long time depending on the size of your data.
```{r, message = F, warning=F}
dds <- DESeq(dds)
```

Once the predictions are complete, you can view the results of the models by using the function `results()` on your DDS object. Recall above I had mentioned that there were 2 ways to re-level factors. I showed you the first method above, and now I am showing you the second method. Using the function, you pass the second parameter for `contrast` set to a vector that contains the column "Treatment" followed by the order of treatment types to compare against. In this case, we want to compare the "Infected" ducks against the "Control" ducks.
```{r}
res <- results(dds, contrast=c("Treatment", "Infected", "Control"))
```

Let's take a look at the results (details about the treatment comparison). The first line beginning with "log2 fold change" specifies the gene expression estimates are of the logarithmic fold change log2(infected/control). The second line indicates the statisic test used (Wald test) and the third line represents the number of genes (5738) with 6 columns. The Wald test p-values are shown in the "pvalue" column. The p-values have been adjusted (corrected) using the Benjamini-Hochberg (BH) method and those p-values can be seen in the "padj" column.
```{r}
res[1:4,]
colnames(res)
```

In the following section, I will show you some methods to explore the results and understand what the data is presenting.

### Exploratory Data Analysis (EDA)

```{r, eval = F}
vsd <- vst(dds_all_samples, blind=FALSE)
```

empty

# Ducks Week 4 DEA

empty

# Ducks Week 12 DEA

empty

# Next Steps

empty

# References

empty
