# Integrating Expert Knowledge and Recursive Bayesian Inference

The scripts available in this repository allow for the reproduction of the examples presented in the article "Integrating Expert Knowledge and Recursive Bayesian Inference: A Framework for Spatial and Spatio-Temporal Data Challenges". They are organized as follows:

1. [**Data integration with spatial change of support**](./scripts/runme_Spatial_ChangeOfSupport.R): Example of the integration of heterogeneous data sources, where scarce geostatistical data are combined with areal data. Using two different support structures for the areal data, we show how the approach to dealing with spatial change of support improves the spatial analysis of the data.
2. [**Data integration with categorical change of support**](./scripts/runme_Categorical_aggregation.R): Example of integrating two different data sources where a categorical explanatory variable is aggregated in one of the simulated datasets. In this example, we show that, under certain conditions, it is possible to handle categorical change of support (aggregation/disaggregation of categorical levels) and improve the overall estimation.
3. [**Spatio-Temporal Temperature Example with Recursive Inference**](./scripts/runme_TemperatureSpatiotemporal_Recursive.R): Example of a Big Data problem with temperature data. The dataset consists of 308 spatial locations per temporal node (month) and 480 temporal nodes in total. For this real dataset, we implement the recursive approach in two cases: one using only the first 120 temporal nodes, and a second using all 480 nodes. The first case is lower dimensional to allow comparison with the standard full-data analysis, which cannot be performed for the full dataset. The code reproduces the results and illustrates the potential of the recursive approach grounded in the INLA framework.

## Supplementary material 

The supplementary material of the article includes additional results from the **Spatio-Temporal Temperature Example with Recursive Inference** and an additional simulated example with spatio-temporal structure. The code to reproduce these results is provided at the following link:

1. [**Spatio-Temporal Simulation Example with Recursive Inference**](./scripts/runme_Simulated_Spatiotemporal_Recursive.R): Simple simulated example demonstrating an implementation of the recursive approach. This example allows tracking how the information changes as successive data partitions (new data arrivals) are incorporated.
