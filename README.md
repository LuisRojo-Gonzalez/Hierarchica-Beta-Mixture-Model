# Hierarchical Beta Mixture Model

This is the implementation of a statistical technique called Mixture models using a Beta distribution to find a meta-distribution composed by others via weighted sub-distributions.

In particular, we use the Beta distribution to handle proportions as its domain relies in [0, 1]. Basically, we want to look for a pseudo-joint distribution. The motivation of this project is to find the distribution of transition probabilities part of the dynamic system in obesity prevalence, i.e., answer the kind of questions of "what's the probability of a obese person of certain age and sex to be normal-weight, underweight or overweight during next period?"

Thus, we handle a situation as follows

![BetaIllustration](/Document/Figures/BetaIllustration.png)

where solid colored lines represent the sub-distributions and the dotted one is the mixture or weighted distribution, while the dashed line represents the real distribution to be fitted. Here, as we are handling proportions, we know that the meta-distribution is a Dirichlet with three parameters since there are three possible transitions: i) To decrease, ii) To remain, iii) To increase. Then, as we want to find the sub-distributions, as well as their weighting, we use the close relationship between the Dirichlet distribution and the Beta distribution, being the former a generalization of the latter. The main property we take advantage is that each realization of a Dirichlet distribution sums one, which is what we want over time.

Then, as we have a high nonlinear optimization problem, we use proper optimization algorithms to maximize the goodness-of-fit proposed as the common area between the mixture or weighted distribution and the meta-distribution (which turns out to be a Beta distribution). In other words, we have a Beta distribution that represents how the probability a person of certain sex will suffer a certain transition between two periods, and we aim to know how this decomponse for each nutritional status. In this regard, we perform several computational experiments and see that the selected scheme outperforms a single optimization setting experimentation, i.e., a single setting would not find the best solution for all the combination, as follows

![Trace](/Document/Figures/trace.png)

Finally, selecting the best fitting, we obtain the following mixture models for each combination of sex and transition

![Density](/Document/Figures/mix_density.png)

Preliminary results were presented at ALIO-EURO 2022 conference. The detailed description of the project is public under [this proceeding document](http://dx.doi.org/10.48786/alioeuro.2022.13).

For an illustration of the project, the presentation I made is freely available [here](https://drive.google.com/file/d/1094s6Y7ExIq4wqb-4iMAtmmN_ZKbSCWX/view).
