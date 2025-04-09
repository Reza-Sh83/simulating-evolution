# %% [markdown]
# # <span style="color: #2E86C1;">Simulating Evolution by Mutation, Selection, and Genetic Drift</span>
# 
# ---
# 
# ## <span style="color: #2E86C1;">What is Evolution?</span>
# 
# When biologists discuss evolution, they use a more precise and inclusive definition than we might in casual conversation.
# 
# > <span style="color: #E74C3C;">**Evolution in biology**</span> is defined as **change over time in the heritable characteristics of a species or population**.
# 
# It might not be very obvious just from the definition, but this differs a lot from how many people think of evolution.
# 
# Here are some key differences:
# 
# - <span style="color: #E74C3C;">**No Progress or Perfection**</span>: This definition says nothing about "progress" or "perfection." Any change in traits over time counts equally as evolution.
# - <span style="color: #E74C3C;">**No Complexity Requirement**</span>: Gain or loss of a trait both equally count as evolution.
# - <span style="color: #E74C3C;">**Timeframe**</span>: Evolution can happen over short or long periods of time. Both count equally as evolution.
# - <span style="color: #E74C3C;">**Population Size**</span>: Evolution can happen within small populations or in whole species.
# - <span style="color: #E74C3C;">**Mechanisms**</span>: Change in characteristics need not necessarily be due to natural selection to count as evolution (e.g., genetic drift).
# 
# ---
# 
# ## <span style="color: #2E86C1;">Main Mechanisms for Evolution</span>
# 
# Several mechanisms can produce evolution. The most common are:
# 
# ### <span style="color: #E74C3C;">1. Natural Selection</span>
# Change in a trait because it increased or decreased fitness (fitness is usually measured as the number of surviving offspring). A trait will evolve by natural selection if and only if **three criteria** are met:
# 1. <span style="color: #E74C3C;">**Variation**</span>: The trait varies in a population.
# 2. <span style="color: #E74C3C;">**Heritability**</span>: The trait is inherited from parents to offspring.
# 3. <span style="color: #E74C3C;">**Selection**</span>: The trait non-randomly alters the fitness of organisms that carry it.
# 
# > <span style="color: #E74C3C;">**Fitness**</span> is a measure of how well an individual can pass on its genes to the next generation.
# 
# - **Neutral Traits**: Traits that do not affect fitness are called **neutral traits** or **selectively neutral traits**. These cannot evolve by natural selection unless they become non-neutral due to environmental or organismal changes.
# 
# ---
# 
# ### <span style="color: #E74C3C;">2. Mutation</span>
# Random changes in the DNA sequence that can include:
# - Substitution of **single nucleotides**.
# - Insertions, deletions, or rearrangements of longer stretches of DNA.
# 
# Even if these changes do not affect *fitness*, mutation itself may still change some traits. Mutation rates in nature are generally very low, so mutation is considered a slow mechanism of evolution. However, mutation can combine with natural selection to generate rapid trait changes.
# 
# ---
# 
# ### <span style="color: #E74C3C;">3. Genetic Drift</span>
# If populations were infinitely large, chance would play no role in evolution. However, all real populations have a limited number of individuals, so **sampling effects** can influence evolution. The smaller the population, the more powerful genetic drift will tend to be.
# 
# > <span style="color: #E74C3C;">**Key Point**</span>: Genetic drift occurs in all non-infinite populations and is separate from migration.
# 
# ---
# 
# ### <span style="color: #E74C3C;">4. Gene Flow</span>
# If two populations interbreed and have different traits, gene flow can result in altered traits in both populations. This counts as evolution under our definition.
# 
# - **Horizontal Gene Transfer (HGT)**: In microorganisms (e.g., bacteria), genes can be transferred between vastly different organisms, altering characteristics and driving evolution.
# 
# ---
# 
# ## <span style="color: #2E86C1;">Haploid vs Diploid</span>
# 
# The terms **haploid** and **diploid** describe the number of sets of chromosomes in a cell and are fundamental concepts in genetics and biology.
# 
# ---
# 
# ### <span style="color: #E74C3C;">Haploid</span>
# - **Definition**: A haploid cell contains **one complete set of chromosomes** (n).
# - **Examples**:
#   - **Gametes (sperm and egg cells)** in humans and most animals are haploid.
#   - In humans, the haploid number is \( n = 23 \), meaning each gamete has 23 chromosomes.
# - **Function**: Haploid cells are crucial for sexual reproduction. When two haploid gametes fuse during fertilization, they form a diploid zygote.
# 
# ---
# 
# ### <span style="color: #E74C3C;">Diploid</span>
# - **Definition**: A diploid cell contains **two complete sets of chromosomes** (2n), one from each parent.
# - **Examples**:
#   - Most somatic (body) cells in humans and other animals are diploid.
#   - In humans, the diploid number is \( 2n = 46 \), meaning each somatic cell has 46 chromosomes (23 pairs).
# - **Function**: Diploid cells allow for genetic diversity through recombination during meiosis and maintain the organism's full genetic information.
# 
# ---
# 
# ### <span style="color: #E74C3C;">Contrast</span>
# 
# | Feature                | Haploid (n)                  | Diploid (2n)                  |
# |------------------------|------------------------------|-------------------------------|
# | **Chromosome Sets**    | One set (e.g., 23 in humans) | Two sets (e.g., 46 in humans) |
# | **Examples**           | Gametes (sperm, egg)         | Somatic cells (skin, muscle)  |
# | **Formation**          | Produced by meiosis          | Produced by mitosis           |
# | **Purpose**            | Sexual reproduction          | Growth, repair, and function  |
# 
# ---
# 
# ### <span style="color: #E74C3C;">Example in Life Cycle</span>
# 1. **Human Reproduction**:
#    - Haploid sperm (\( n = 23 \)) fertilizes a haploid egg (\( n = 23 \)).
#    - The resulting zygote is diploid (\( 2n = 46 \)).
# 2. **Plant Life Cycle** (Alternation of Generations):
#    - Haploid spores (\( n \)) grow into a gametophyte, which produces gametes.
#    - Gametes fuse to form a diploid zygote, growing into a sporophyte (\( 2n \)).
# 
# ---
# 
# ## <span style="color: #2E86C1;">Alleles</span>
# 
# ### <span style="color: #E74C3C;">Definition</span>
# An **allele** is a variant form of a gene. Genes are segments of DNA that code for specific traits, and alleles are the different versions of those genes. Each individual inherits two alleles for each gene—one from each parent.
# 
# ---
# 
# ### <span style="color: #E74C3C;">Key Points</span>
# - **Location**: Alleles are found at the same position (locus) on homologous chromosomes.
# - **Variation**: Differences in alleles arise due to mutations or genetic variation, leading to different traits.
# - **Dominant and Recessive**:
#   - A **dominant allele** expresses its trait even if only one copy is present.
#   - A **recessive allele** expresses its trait only when two copies are present.
# 
# ---
# 
# ### <span style="color: #E74C3C;">Examples</span>
# 1. **Human Eye Color**:
#    - The gene for eye color has multiple alleles, such as those for brown, blue, or green eyes.
#    - Brown is often dominant over blue or green.
# 2. **Blood Type**:
#    - The ABO blood group system is determined by three alleles: \( I^A \), \( I^B \), and \( i \).
#    - \( I^A \) and \( I^B \) are dominant, while \( i \) is recessive.
# 
# ---
# 
# ### <span style="color: #E74C3C;">Homozygous vs Heterozygous</span>
# - **Homozygous**: An individual has two identical alleles for a gene (e.g., AA or aa).
# - **Heterozygous**: An individual has two different alleles for a gene (e.g., Aa).
# 
# ---
# 
# ### <span style="color: #E74C3C;">Significance</span>
# Alleles are responsible for genetic diversity and the inheritance of traits. They play a key role in evolution through natural selection, as certain alleles may confer advantages in specific environments.
# 
# ---
# 
# ### <span style="color: #E74C3C;">Visual Example</span>
# 
# | Gene           | Alleles                     | Trait Expressed               |
# |-----------------|-----------------------------|--------------------------------|
# | Eye Color      | Brown (B), Blue (b)         | \( Bb \): Brown eyes (dominant) |
# | Blood Type     | \( I^A \), \( I^B \), \( i \) | \( I^A i \): Type A blood       |
# 
# ---
# 
# ## <span style="color: #2E86C1;">Selection on Haploids</span>
# 
# Let's start by simulating natural selection in a **haploid** organism.
# 
# Haploid organisms make a nice starting place for our simulations because we only have to think about selection alone, without having to also integrate an understanding of genetics (i.e., the probabilities involved with different combinations of gene copies being passed on to offspring).
# 
# Let's start with an example involving bacteria that have two different alleles (variants or versions) of a particular gene. Let's further assume that one of these alleles confers resistance to an antibiotic that would usually prevent our bacterium from reproducing.
# 
# We'll call those two alleles: allele A1 and allele A2. Because our organisms are haploid, each will *either* have allele A1 *or* allele A2, but not both. We'll further assume that these are the *only* two alleles at this spot or *locus* in the genome. We will represent the **frequency** of allele A1 with the variable A1 and the frequency of allele A2 with the variable A2. The frequency of an allele is simply the proportion of all alleles at that locus (e.g., in that gene) that it represents.
# 
# So \( A1 = 0.6 \) would indicate that allele A1 represents 60% of the alleles at that locus in the population. Since there are only these two alleles at this locus in the population, \( A1 + A2 = 1.0 \). Therefore, if \( A1 = 0.6 \), \( A2 \) would be \( 0.4 \).

# %% [markdown]
# # <span style="color: #2E86C1;">Explanation of Code</span>
# ## <span style="color: #E74C3C;">Import Necessary Functions</span>

# %%
from collections import defaultdict
import matplotlib.pyplot as plt
from random import random
import numpy as np
from matplotlib.cm import viridis,hsv,Spectral
from matplotlib.pyplot import style

#Make plots have a dark background
style.use('dark_background')

# %% [markdown]
# ## Define a function for how the allele frequencies $A_1$ and $A_2$ will change in a generation
# 
# Let's imagine what happens after 1 generation. Each $A_1$ individual will have some chance of surviving, and some average number of surviving offspring if it does survive to reproduce. Together, this defines the *fitness* of the $A_1$ genotype. $A_2$ individuals will have a different chance of surviving to reproduce and/or average number of surviving offspring.
# 
# In the next generation, the proportion of $A_1$ genotypes will depend on the ratio of $A_1$ genotypes in the previous generation, times the average fitness of each $A_1$ individual, and divided by the same for the total pool of all individuals (both $A_1$ and $A_2$).
# 
# The proportion of $A_1$ genotypes in the next generation $p'_{A_1}$ is given by:
# 
# $$
# p'_{A_1} = \frac{p_{A_1} \cdot w_{A_1}}{\overline{w}}
# $$
# 
# Where:
# - $p_{A_1}$: Proportion of $A_1$ genotypes in the current generation.
# - $w_{A_1}$: Average fitness of $A_1$ individuals.
# - $\overline{w}$: Average fitness of the entire population, calculated as:
#   $$
#   \overline{w} = p_{A_1} \cdot w_{A_1} + p_{A_2} \cdot w_{A_2}
#   $$
#   - $p_{A_2}$: Proportion of $A_2$ genotypes in the current generation.
#   - $w_{A_2}$: Average fitness of $A_2$ individuals.
# 
# 
# **Example:**
# $A_1$ = 0.90
# fitness $A_1$ = 0.0
# $A_2$ = 0.10
# fitness $A_2$ = 1.0
# 
# Here no $A_1$ alleles survive, so as long as the fitness of $A_2$ and it's starting abundance are >0, the next generation should be all $A_2$.
# 
# Expected output: 0.0
# 
# **Example 2:**
# $A_1$ = 0.60
# $A_2$ = 0.40
# fitness $A_1$ = 1.0
# fitness $A_2$ = 1.0
# 
# Here the alleles are in unequal proportion but have equal fitness. Since neither has a fitness advantage, and we are assuming an infinite population size (so no randomness), $A_1$ and $A_2$ will not change in the next generation.
# 
# 
# Here is how we might implement that in code:
# 

# %%
def selection_on_haploid(A1,A2,fitness_A1,fitness_A2):
    """Return the allele frequency of allele A1 in haploids \
    in the presence of selection

    A1 -- the frequency of the A1 allele as a float
    A2 -- the frequency of the A2 allele as a float
    fitness_A1 -- the relative fitness of A1 (float).
    fitness_A2 -- the relative fitness of A2 (float).

    Since we're dealing with proportions, the scale of the fitnesses
    is arbitrary, but the ratio matters.
    """

    numerator = A1 * fitness_A1
    denominator = A1 * fitness_A1 + A2 * fitness_A2

    return numerator/denominator

def simulate_generations(initial_conditions, fitness_A1, fitness_A2, generations):
    """
    Simulate allele frequency changes over multiple generations for multiple initial conditions
    and plot the results.

    Parameters:
    - initial_conditions: List of tuples containing (A1_initial, A2_initial) for each condition.
    - fitness_A1: Relative fitness of A1 (float).
    - fitness_A2: Relative fitness of A2 (float).
    - generations: Number of generations to simulate (int).
    """
    plt.figure(figsize=(10, 6))
    colors = ["blue", "green", "red", "yellow", "white"]  # Colors for different initial conditions

    for idx, (A1_initial, A2_initial) in enumerate(initial_conditions):
        A1_frequencies = [A1_initial]  # Store A1 frequencies over generations
        A2_frequencies = [A2_initial]  # Store A2 frequencies over generations

        for _ in range(generations):
            A1_new = selection_on_haploid(A1_frequencies[-1], A2_frequencies[-1], fitness_A1, fitness_A2)
            A2_new = 1 - A1_new  # Since A1 + A2 = 1
            A1_frequencies.append(A1_new)
            A2_frequencies.append(A2_new)

        # Plot the results for this initial condition
        plt.plot(range(generations + 1), A1_frequencies, label=f"A1 (Initial = {A1_initial})", color=colors[idx], linestyle="-")
        plt.plot(range(generations + 1), A2_frequencies, label=f"A2 (Initial = {A2_initial})", color=colors[idx], linestyle="--")

    # Customize the plot
    plt.xlabel("Generations")
    plt.ylabel("Allele Frequency")
    plt.title("Change in Allele Frequencies Over Generations for Different Initial Conditions")
    plt.legend()
    plt.show()

# Define initial conditions as a list of tuples (A1_initial, A2_initial)
initial_conditions = [
    (0.1, 0.9),  # Condition 1: A1 = 0.1, A2 = 0.9
    (0.3, 0.7),  # Condition 1: A1 = 0.3, A2 = 0.7
    (0.5, 0.5),  # Condition 2: A1 = 0.5, A2 = 0.5
    (0.7, 0.3),  # Condition 3: A1 = 0.7, A2 = 0.3
    (0.9, 0.1),  # Condition 1: A1 = 0.9, A2 = 0.1
]

# Define fitness values
fitness_A1 = 1.0  # Fitness of A1
fitness_A2 = 0.8  # Fitness of A2

# Number of generations to simulate
generations = 50

# Run the simulation and plot the results
simulate_generations(initial_conditions, fitness_A1, fitness_A2, generations)

# %% [markdown]
# ## Modelling  the spread of an allele for antibiotic resistance in a bacterial population
# 
# Let us imagine a haploid bacterium with two alleles, R and S, at a locus. Allele R confers resistance to an antibiotic. When a patient is treated with the antibiotic, 0% of individuals with the R genotype die, but 50% of indiviudals with the S genotype do.
# 
# Let's imagine we start with only 10% R bacteria and 90% S bacteria. How many generations will it take for the R bacteria take over in a patient being treated with this antibiotic?
# 
# To make the problem simpler at first we'll also assume:
#   * There are no mutations.
#   * The population size is infinte.
#   * There is no migration into or out of the population
#   * Generations are non-overlapping (we can model them as discrete 'turns').
# 
# **Approach:** We'll set our initial frequencies and fitnesses, and empty lists to hold results, then use a for loop to run our function for selection in a haploid once per generation. As we do so we'll take care to update the current allele freqency -- f(R) -- in each generation. Since there are only two alleles f(S) + f(R) = 1.0, and we will be able to calculate f(S) as 1.0 - f(R). Each generation we'll append the result to a list, which we will graph at the end.

# %%
def simulate_selection_on_haploid(f_R, W_R, W_S, generations=100):
    """
    Simulate the change in allele frequencies in a haploid population over time.

    Parameters:
    f_R (float): Initial frequency of the resistant allele (R).
    W_R (float): Relative fitness of the resistant allele (R).
    W_S (float): Relative fitness of the susceptible allele (S).
    generations (int): Number of generations to simulate.

    Returns:
    tuple: (generations_list, resistant_frequencies, susceptible_frequencies)
    """
    # Initialize variables
    curr_allele_frequency = f_R
    generations_list = list(range(generations))
    resistant_frequencies = []
    susceptible_frequencies = []

    # Simulate allele frequency changes across generations
    for _ in range(generations):
        f_S = 1.0 - curr_allele_frequency
        curr_allele_frequency = (curr_allele_frequency * W_R) / (curr_allele_frequency * W_R + f_S * W_S)
        resistant_frequencies.append(curr_allele_frequency)
        susceptible_frequencies.append(f_S)

    return generations_list, resistant_frequencies, susceptible_frequencies

# Define multiple initial conditions
initial_conditions = [0.1, 0.3, 0.5, 0.7, 0.9]  # Initial frequencies of R
W_R = 1.0  # Fitness of R
W_S = 0.8  # Fitness of S
generations = 100

# Plot settings
plt.figure(figsize=(20, 8))
colors = ["blue", "green", "red", "yellow", "white"]  # Colors for different initial conditions

# Simulate and plot for each initial condition
for idx, f_R in enumerate(initial_conditions):
    generations_list, resistant_frequencies, susceptible_frequencies = simulate_selection_on_haploid(
        f_R, W_R, W_S, generations
    )
    plt.plot(generations_list, resistant_frequencies, label=f'Resistant (R), Initial f_R = {f_R}',
             color=colors[idx], linestyle="-", linewidth=2)
    plt.plot(generations_list, susceptible_frequencies, label=f'Susceptible (S), Initial f_R = {f_R}',
             color=colors[idx], linestyle="--", linewidth=2)

# Customize the plot
plt.title(f"Allele Frequency Dynamics (W_R = {W_R}, W_S = {W_S})", fontsize=16, pad=20)
plt.xlabel('Generations', fontsize=14)
plt.ylabel('Allele Frequency', fontsize=14)
plt.legend(fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')  # Legend outside the plot
plt.tight_layout()  # Adjust layout to prevent overlap
plt.show()

# %% [markdown]
# ### Explanation of the Results
# 
# When we simulate the spread of the resistant allele (R) in a bacterial population under antibiotic treatment, the following scenarios can occur based on the input parameters:
# 
# 1. **Scenario 1: Initial Conditions**
#    - The simulation starts with 10% resistant bacteria (f_R = 0.1) and 90% susceptible bacteria (f_S = 0.9).
#    - The resistant allele has a fitness advantage (W_R = 1.0) compared to the susceptible allele (W_S = 0.5).
# 
# 2. **Scenario 2: Selection Pressure in Favor of R**
#    - The resistant allele frequency increases generation by generation because R confers an advantage in an environment with antibiotics.
#    - The susceptible allele (S) becomes less common due to reduced survival and reproduction rates.
# 
# 3. **Scenario 3: Fixation of R**
#    - Over time, R becomes dominant in the population, approaching a frequency of 1.0 (fixation), while S is almost entirely eliminated.
# 
# 4. **Scenario 4: No Change in Frequencies**
#    - If W_R = W_S, there would be no selection. Both alleles would maintain their initial frequencies over time.
# 
# 5. **Scenario 5: Selection Against R**
#    - If W_R < W_S, the resistant allele would decrease in frequency, eventually being outcompeted by the susceptible allele.
# 
# ### Key Insights:
# - The strength of selection (difference in W_R and W_S) directly impacts how quickly R becomes fixed in the population.
# - The initial frequency of the resistant allele (f_R) determines how long it takes for R to dominate.
# - If the fitness values are equal, no evolutionary change occurs.
# 
# These results highlight the powerful effect of selection pressure in shaping allele frequencies within populations.
# 

# %% [markdown]
# ## Selection in a diploid organism
# 
# When Darwin formulated the theory of evolution by natural selection, the source of variation was a tough problem. Several solutions proposed at the time -- like Darwin's notion (following Lamark) that the environment modified an organisms genetics and those changes could be inherited -- were tested and abandoned.
# 
# The **modern synthesis** in the 1930s combined Darwin's theory of evolution by natural selection with **mutation** as the source of biological variation, and inheritance following the **genetic** rules of inheritance derived from Gregor Mendel's work on inheritance in pea plants.
# 
# As this next section requires a bit more theory, let us talk briefly about the key concepts we'll need to extend our model to diploid organisms. Diploid organisms, as we noted above, by definition have two copies of each chromosome. Because all genetic information is encoded on a physical chromosome, this implies that diploid organisms have two copies of each allele (I am ignoring for the moment some important exceptions like sex chromosomes). During sexual reproduction, one copy comes from each parent.
# 
# This means that -- unlike in haploids -- we cannot represent the **genotype** of a particular organism with just R or S. It must instead either be RR, RS, or SS. Organisms that have two copies of the same allele -- like individuals that are RR or SS -- are called *homozygous* because they have the same (*homo*) allele on each chromosome. RS organisms are called *heterozygous* because they have different (*hetero*) alleles on the two chromosomes.
# 
# ## Moving from allele frequencies to genotype frequencies using Hardy-Weinberg Equilibrium
# 
# Under certain simplifying assumptions, called **Hardy-Weinberg Equilibrium** it is possible to predict the frequency of each genotype given the frequency of the alleles that compose it. If the conditions of Hardy-Weinberg Equilibrium hold for a set of two alleles at the same locus in a diploid organism, and if $p$ represents the frequency of the allele $A_1$ and $q$ represents the frequency of allele $A_2$, then the frequencies of the genotypes will be:
# 
# $$ freq(A_1A_1) = p^2$$
# 
# $$ freq(A_1A_2) = 2pq$$
# 
# $$ freq(A_2A_2) = q^2$$
# 
# Because the above list comprises all possible genotypes, it follows that the genotype frequencies above sum to 1.0:
# 
# $$p^2 + 2pq + q^2 = 1.0$$
# 
# ## Modelling how genotypes translate into phenotypes
# 
# Critically, selection does not act on the genotype directly. Instead, that genotype -- in interaction with its environmental context -- is translated into a phenotype or set of realized characters . The way in which different genotypes translate into phenotypes is called the **mode of inheritance**.
# 
# For our purposes of modelling selection in diploids, we care most about how the fitness values $W_{pp}$, $W_{pq}$ and $W_{qq}$ differ under different modes of inheritance.
# 
# **Dominance** -- The simplest type of inhertiance - the one first discovered by *Gregor Mendel* in pea plants - is called Complete Mendelian Dominance. In this mode of inheritance, one allele is dominant over another. In heterozygotes, the dominant allele determines the phenotype. In our above example if R were dominant to S, then RS heterozygotes would show the phenotype -- the same level of antibiotic resistance -- as RR individuals - even though they only have one R allele. Because RS and RR individuals both show the same phenotype, we would expect them to have the same fitness. So if R showed complete Mendelian Domiance over S, we would expect $W_{pp}$ == $W_{pq}$.
# 
# Other modes of inheritance are described below:
# 
# **Recessive** -- heterozygotes show the same phenotype as homozygotes for the other allele. So if R is recessive to S, then RS heterozygotes would be identical to SS homozygotes, and so be susceptible to the antibiotic. Therefore we expect $W_{pq}$ == $W_{qq}$ (the opposite of Dominance)
# 
# **Incomplete Dominance** -- in incomplete dominance, RS heterozygotes show an intermediate phenotype between RR and SS. If we observe $W_{pp}$ > $W_{pq}$ > $W_{qq}$ we can rule out a purely dominant or recessive trait (since in that case $W_{pq}$ == $W_{pp}$ or $W_{pq}$ == $W_{qq}$), and might instead suspect incomplete dominance or codominance (below)
# 
# **Codominance** -- if two alleles are co-dominant, then heterozygotes show a mosaic or patchwork of both phenotypes. For example, for an allele governing color codominance might result in a patchy pattern.
# 
# This is far from an exclusive list of patterns of inheritance, but it should at least help to illustrate how Mendelian inheritance works. It is worth noting at this point that the above is a useful simplification - in reality more complex patterns of inheritance are possible, and an organisms environment (not just its genetics) can influence how many genotypes translate into phenotypes.
# 
# ### Summary
# 
# Diploids have two copies of each genetic variant or **allele**. The set of alleles is called the **genotype**. These genotypes often map to phenotypes -- the actual observed traits of the organism -- following a pattern of inheritance (dominant, recessive, incomplete dominance or codominance). **Natural selection acts on phenotypes not genotypes**, so if we want to understand fitness we have to also understand how genotypes translate into phenotypes.
# 
# ## Designing the simulation.
# 
# We now have the basic background we need to develop a more sophisticated simulation of the evolution of our antibiotic resistance allele if the organism is diploid.
# 
# **Goal**: determine the allele frequency of an allele in a population over the generations. We will assume that criteria for Hardy-Weingberg equilibrium are met *except* for the presence of selection.
# 
# **Procedure**:
# 
# 1. Write a function to simulate the change in allele frequency in one generation,
# 2. Run that function over the desired number of generations to generate a list of  allele frequencies over time,
# 3. Graph that result.
# 
# 
# ### Writing a function for the allele frequency change in one generation ###
# 
# To review, we previously calculated how selection would change the frequency of allele R in haploids by multiplying the frequency of the R genotype by its fitness and then dividing by the same for each allele.  Now we will do the same again- the only trick is that now there are 3 possible genotypes (RR, RS, or SS) and 3 possible fitnesses $(W_{rr}, W_{rs}, or W_{ss})$ instead of just 2 possible genotypes (R or S) and 2 possible fitnesses $W_r$  and  $W_s$.
# 
# All this is summarized in the equation for Hardy-Weinberg equilibrium, modified to account for selection:
# 
# $$
# p'   =  \frac{p^2W_{pp} + (1/2)2pqW_{pq}}{p^2W_{pp} + 2pqW_{pq} + q^2W_{qq}}
# $$
# 
# (Note that  the RS genotype contributes only half as many R alleles as the RR genotype to the next generation, so the numerator multiplies the term for the number of surviving RS genotypes by 1/2)
# 
# This equation can be written as a python function as follows:

# %%
def HW_with_selection(p, q, Wpp, Wpq, Wqq):
    """
    Calculate the new frequency of allele A1 under Hardy-Weinberg equilibrium with selection.

    Parameters:
    p (float): Frequency of allele A1.
    q (float): Frequency of allele A2.
    Wpp (float): Fitness of A1A1 homozygotes.
    Wpq (float): Fitness of A1A2 heterozygotes.
    Wqq (float): Fitness of A2A2 homozygotes.

    Returns:
    float: New frequency of allele A1.
    """
    numerator = p**2 * Wpp + p * q * Wpq
    denominator = p**2 * Wpp + 2 * p * q * Wpq + q**2 * Wqq
    return numerator / denominator

# %% [markdown]
# **Run the simulation over 100 generations**. Now that we can get allele frequencies in the next generation from the previous one, we can run the simulation over many generations in a for loop:

# %%
def simulate_HW_with_selection_over_time(p, q, Wpp, Wpq, Wqq, generations=100):
    """
    Simulate the change in allele frequencies under Hardy-Weinberg equilibrium with selection over generations.

    Parameters:
    p (float): Initial frequency of allele A1.
    q (float): Initial frequency of allele A2.
    Wpp (float): Fitness of A1A1 homozygotes.
    Wpq (float): Fitness of A1A2 heterozygotes.
    Wqq (float): Fitness of A2A2 homozygotes.
    generations (int): Number of generations to simulate.

    Returns:
    tuple: (generations_list, A1_frequencies, A2_frequencies)
    """
    generations_list = list(range(generations))
    A1_frequencies = [p]
    A2_frequencies = [q]

    for _ in range(1, generations):
        new_p = HW_with_selection(A1_frequencies[-1], A2_frequencies[-1], Wpp, Wpq, Wqq)
        new_q = 1.0 - new_p
        A1_frequencies.append(new_p)
        A2_frequencies.append(new_q)

    return generations_list, A1_frequencies, A2_frequencies

# Define multiple initial conditions
initial_conditions = [
    (0.1, 0.9),  # Condition 1: p = 0.1, q = 0.9
    (0.3, 0.7),  # Condition 1: p = 0.3, q = 0.7
    (0.5, 0.5),  # Condition 2: p = 0.5, q = 0.5
    (0.7, 0.3),  # Condition 3: p = 0.7, q = 0.3
    (0.9, 0.1),  # Condition 1: p = 0.9, q = 0.1
]

# Define different sets of fitness values
fitness_sets = [
    {"Wpp": 0.6, "Wpq": 1.0, "Wqq": 0.5},  # Set 1: Heterozygote advantage
    {"Wpp": 1.0, "Wpq": 0.8, "Wqq": 0.6},  # Set 2: A1A1 has highest fitness
    {"Wpp": 0.5, "Wpq": 0.5, "Wqq": 1.0},  # Set 3: A2A2 has highest fitness
]

# Number of generations to simulate
generations = 100

# Plot settings
plt.figure(figsize=(18, 12))  # Larger figure for multiple subplots
colors = ["blue", "green", "red", "yellow", "white"]  # Colors for different initial conditions

# Create a subplot for each set of fitness values
for idx, fitness in enumerate(fitness_sets):
    plt.subplot(2, 2, idx + 1)  # Arrange subplots in a 2x2 grid
    Wpp, Wpq, Wqq = fitness["Wpp"], fitness["Wpq"], fitness["Wqq"]

    # Simulate and plot for each initial condition
    for cond_idx, (p, q) in enumerate(initial_conditions):
        generations_list, A1_frequencies, A2_frequencies = simulate_HW_with_selection_over_time(
            p, q, Wpp, Wpq, Wqq, generations
        )
        plt.plot(generations_list, A1_frequencies, label=f'A1 (Initial p = {p})',
                 color=colors[cond_idx], linestyle="-", linewidth=2)
        plt.plot(generations_list, A2_frequencies, label=f'A2 (Initial q = {q})',
                 color=colors[cond_idx], linestyle="--", linewidth=2)

    # Customize each subplot
    plt.title(f"Fitness: Wpp = {Wpp}, Wpq = {Wpq}, Wqq = {Wqq}", fontsize=14, pad=10)
    plt.xlabel('Generations', fontsize=12)
    plt.ylabel('Allele Frequency', fontsize=12)
    plt.legend(fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')  # Legend outside the plot
    plt.grid(True, linestyle='--', alpha=0.6)  # Add a grid for better readability

# Adjust layout and display the plot
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Results and Explanation of Diploid Selection Simulation
# 
# ### Key Observations:
# 1. **Initial Conditions**:
#    - We start with an initial allele frequency of $p$ (A1) = 0.2 and $q$ (A2) = 0.8.
#    - Fitness values for the three genotypes are set as:
#      - $W_{pp}$ (A1A1 homozygotes): 0.6
#      - $W_{pq}$ (A1A2 heterozygotes): 1.0
#      - $W_{qq}$ (A2A2 homozygotes): 0.5
# 
# 2. **Selection Effects**:
#    - Over time, the allele frequency of A1 ($p$) changes due to differential fitness values.
#    - Heterozygotes (A1A2) have the highest fitness ($W_{pq} = 1.0$), indicating a scenario of **overdominance**, where the heterozygote has a selective advantage over both homozygotes.
# 
# 3. **Equilibrium**:
#    - As generations progress, $p$ and $q$ approach an equilibrium state where the changes in allele frequency stabilize.
#    - This equilibrium represents the balance point where the effects of selection maintain a stable allele frequency.
# 
# 4. **Impact of Fitness Values**:
#    - If $W_{pp} > W_{pq}$, allele A1 would increase in frequency faster, eventually leading to fixation (p = 1.0).
#    - Conversely, if $W_{qq} > W_{pq}$, allele A2 would dominate.
# 
# ### Summary:
# The results illustrate how allele frequencies change over generations in diploid populations under selection. The mode of inheritance and fitness differences between genotypes determine the rate and direction of allele frequency changes. This highlights the importance of genotype-phenotype mapping and selection in evolutionary dynamics.
# 

# %% [markdown]
# ## Influence of Ploidy on the Evolution of Antibiotic Resistance
# 
# **Ploidy**: The number of chromosomes in the nucleus of a cell.
# 
# ### Key Question:
# How does the ploidy of an organism influence the rate at which antibiotic resistance alleles spread in a population?
# 
# ### Simulation Overview:
# We simulated the evolution of an antibiotic resistance allele in both haploid and diploid organisms under the same selection pressures.
# 
# #### Parameters:
# 1. **Diploid Organisms**:
#    - Fitness of A1A1 homozygotes ($W_{pp}$): 1.0
#    - Fitness of A1A2 heterozygotes ($W_{pq}$): 1.0
#    - Fitness of A2A2 homozygotes ($W_{qq}$): 0.50
#    - Initial allele frequencies: $p = 0.1$ (A1), $q = 0.9$ (A2)
#    - Generations simulated: 200
# 
# 2. **Haploid Organisms**:
#    - Fitness of A1 ($W_p$): 1.0
#    - Fitness of A2 ($W_q$): 0.50
#    - Initial allele frequencies: $p = 0.1$ (A1), $q = 0.9$ (A2)
#    - Generations simulated: 200
# 
# #### Results:
# 1. **Comparison of Frequencies**:
#    - The frequency of allele $q$ decreases over time in both haploids and diploids, but the rate of decrease differs.
#    - Diploids exhibit a slower decline in $q$ frequency compared to haploids.
# 
# 2. **Differences in Evolutionary Dynamics**:
#    - A key factor contributing to this difference is that diploids have heterozygotes (A1A2) that can partially mask the effects of selection, slowing the fixation of beneficial alleles.
# 
# ### Conclusion:
# Ploidy significantly influences the evolutionary trajectory of antibiotic resistance alleles. In diploids, heterozygotes act as a buffer, slowing the impact of selection. This results in different dynamics compared to haploids, where allele frequencies shift more rapidly due to direct selection on phenotypes.
# 
# Below is the Python code used to simulate and visualize these dynamics.

# %%
# Simulate diploid selection
generations = 200
Wpp = 1.0  # Fitness of A1A1 homozygotes
Wpq = 1.0  # Fitness of A1A2 heterozygotes
Wqq = 0.50  # Fitness of A2A2 homozygotes
p = 0.1  # Initial frequency of A1
q = 0.9  # Initial frequency of A2

xs, ps_diploid, qs_diploid = simulate_HW_with_selection_over_time(p, q, Wpp, Wpq, Wqq, generations)

# Simulate haploid selection
def simulate_selection_on_haploid(p: float, q: float, Wp: float, Wq: float, generations: int):
    """
    Simulate allele frequency dynamics in haploids under selection.

    Parameters:
    p (float): Initial frequency of allele A1.
    q (float): Initial frequency of allele A2.
    Wp (float): Fitness of allele A1.
    Wq (float): Fitness of allele A2.
    generations (int): Number of generations to simulate.

    Returns:
    tuple: (generations_list, A1_frequencies, A2_frequencies)
    """
    generations_list = []
    A1_frequencies = []
    A2_frequencies = []

    for generation in range(generations):
        p_prime = (p * Wp) / (p * Wp + q * Wq)
        q_prime = 1.0 - p_prime
        generations_list.append(generation)
        A1_frequencies.append(p_prime)
        A2_frequencies.append(q_prime)
        p, q = p_prime, q_prime

    return generations_list, A1_frequencies, A2_frequencies


Wp = 1.0  # Fitness of allele A1
Wq = 0.50  # Fitness of allele A2
p = 0.1  # Initial frequency of A1
q = 0.9  # Initial frequency of A2

xs, ps_haploid, qs_haploid = simulate_selection_on_haploid(p, q, Wp, Wq, generations)

# Calculate the difference in q frequency between haploids and diploids
freq_diff = [qs_haploid[i] - qs_diploid[i] for i in range(len(xs))]

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(xs, qs_diploid, '-r', label="Diploid")
plt.plot(xs, qs_haploid, '-b', label="Haploid")
plt.title("Frequency Dynamics of Antibiotic Resistance Allele (q)")
plt.xlabel("Generations")
plt.ylabel("Frequency of allele q")
plt.legend()
plt.show()

# %% [markdown]
# ## Visualizing the Impact of Ploidy on Allele Frequency Dynamics
# 
# ### Key Question:
# How does ploidy influence the rate of change in allele frequency, and how can we quantify the difference?
# 
# ### Visualization:
# We calculated the difference in the frequency of the antibiotic-sensitive allele $q$ between haploids and diploids for each generation. This difference ($f(q)_{\text{haploid}} - f(q)_{\text{diploid}}$) reflects how much faster or slower the allele $q$ decreases in haploids relative to diploids.
# 
# ### Insights:
# 1. A positive value of the difference indicates that the frequency of $q$ is higher in haploids than in diploids at that generation.
# 2. The magnitude of the difference highlights the buffering effect of diploidy due to heterozygotes (A1A2), which slows down the elimination of allele $q$.
# 
# ### Conclusion:
# This visualization underscores the importance of ploidy in shaping evolutionary dynamics. Diploid populations evolve differently from haploid populations due to the additional layer of genetic complexity introduced by heterozygosity.
# 
# Below is the Python code used to generate the visualization.

# %%
# Visualize the difference due to ploidy
plt.figure(figsize=(10, 6))
plt.plot(xs, freq_diff, '-y', label="Difference due to ploidy", linewidth=5)
plt.xlabel("Generations")
plt.ylabel("Difference in q frequency\n(f(q) for haploid - f(q) for diploid)")
plt.title("Impact of Ploidy on Allele Frequency Dynamics")
plt.legend()
plt.show()

# %% [markdown]
# ## Simulating Genetic Drift
# 
# ### What is Genetic Drift?
# Genetic drift is a stochastic process that leads to random changes in allele frequencies, especially in small populations. It arises due to sampling effects during reproduction.
# 
# ### Why Does Population Size Matter?
# The Hardy-Weinberg equilibrium assumes infinite population sizes, where $f(A_1A_1)$ equals $f(A_1)f(A_1)$. However, in small populations, allele frequencies fluctuate due to chance, making predictions based on Hardy-Weinberg less precise.
# 
# ### Simulation Overview:
# We simulate random mating in finite populations with no selection, mutation, or migration. By varying the population size, we can observe how genetic drift becomes more pronounced in smaller populations.
# 
# ### Key Insights:
# 1. **Small populations:** Allele frequencies fluctuate wildly and may fixate (reach 0 or 1) in fewer generations.
# 2. **Large populations:** Drift is less pronounced, and allele frequencies remain closer to their initial values over time.
# 
# ### Simulation Procedure:
# 1. Generate offspring genotypes through random mating.
# 2. Compute allele frequencies for each generation.
# 3. Visualize allele frequency dynamics for varying population sizes.
# 
# ### Simulation Results:
# The Python code below simulates genetic drift for different population sizes and generates plots to visualize how allele frequencies change over time.
# 

# %%
# Function to generate a random genotype
def random_genotype(f_A1):
    sperm_allele = 'A1' if random() <= f_A1 else 'A2'
    egg_allele = 'A1' if random() <= f_A1 else 'A2'
    return sperm_allele + egg_allele

# Function to simulate random mating in a finite population
def simulate_random_mating(pop_size, f_A1):
    genotypes = defaultdict(int)
    for _ in range(pop_size):
        genotype = random_genotype(f_A1)
        genotypes[genotype] += 1

    # Normalize genotype frequencies
    for g in genotypes.keys():
        genotypes[g] /= float(pop_size)
    return genotypes

# Function to simulate genetic drift over generations
def simulate_genetic_drift(n_generations, f_A1=0.5, pop_size=100):
    allele_freqs = []
    for _ in range(n_generations):
        genotypes = simulate_random_mating(pop_size, f_A1)
        f_A1 = (
            genotypes['A1A1'] +
            0.5 * (genotypes['A1A2'] + genotypes['A2A1'])
        )
        allele_freqs.append(f_A1)
    return list(range(n_generations)), allele_freqs

# %% [markdown]
# Now let's run the simulation for a variety of population sizes:

# %%
# Run simulations for different population sizes
f_A1 = 0.5
pop_sizes = [5, 50, 500, 5000]
replicates = 20

plt.figure(figsize=(15, 20))
for i, pop_size in enumerate(pop_sizes):
    print(f"Simulating pop_size: {pop_size}")
    # Create a subplot for each population size
    plt.subplot(len(pop_sizes), 1, i + 1)
    plt.ylabel(f"f(A1)\n n = {pop_size}")
    for _ in range(replicates):
        generations = pop_size * 6
        xs, ys = simulate_genetic_drift(generations, f_A1, pop_size)
        plt.plot(xs, ys, '-', linewidth=0.5)
    plt.ylim(0.0, 1.0)
    plt.xlim(0, generations)

# Add a global x-axis label and adjust subplot spacing
plt.xlabel("Generations")
plt.subplots_adjust(hspace=0.2)
plt.show()

# %% [markdown]
# ## Genetic Drift
# 
# Let's look at the figures produced up above. Each graph from top to bottom shows the frequency of allele $A_1$ changing over 500 generations in 20 replicate populations. (Replicate here just means that these are 20 *separate* populations - like flies breeding in different test tubes or lizards on different islands - that don't influence each other at all). The *only difference* between these graphs is the population size. There are several really important things to notice here...
# 
# Yet you can see that the frequency of allele $A_1$ differs quite a bit between these replicates. You may also see that random drift in the frequency of allele $A_1$ is stronger in.
# 
# ### Genetic drift can cause variation
# 
# There are several things to notice here. First, each line takes a somewhat different path. That means that given the same exact parameters (same pop. size, same allele frequency, same selective regime, etc) evolution has taken a slightly different course in each population. Importantly, we've shown that this can happen *without* mutation and *without* migration. Why then is the allele frequency changing over time?
# 
# ### Genetic drift can occur due to sampling effects during sexual reproduction
# 
# We haven't yet introduced any mutation or migration at all into the population so we've shown that genetic drift *does not* depend on mutation or migration. This is a very important point, and one that is often mixed up by folks learning about genetic drift for the first time.
# 
# Yet despite the absence of mutation and migration, each replicate population differs from the others. Those differences compound over generations, causing populations to 'drift' randomly in allele frequency.
# 
# In our simulation the true cause is what's called a sampling effect: while the probability of getting an $A_1A_1$ phenotype when crossing two $A_1A_2$ heterozygotes is 0.25, if we only have 4 offspring it will be fairly common for 0,1,2,3 or even 4 of them to be $A_1A_1$. So the *observed* count of $A_1A_1$ can be anywhere from 0-100%, even though the *expected* frequency is exactly 1. The same thing happens at a population level:  _when the number of observations is anything less than infinite the observed counts for each genotype won't perfectly match the expected counts of each genotype._ But it is the observed genotypes (i.e. the ones that actually occur in the organisms) that will be passed on to the next generation. Therefore this process changes allele frequencies.

# %% [markdown]
# ## Combining Genetic Drift and Selection
# 
# ### Key Concept:
# Natural selection and genetic drift are two forces driving evolution. Genetic drift is a stochastic process that affects allele frequencies in small populations, while selection favors certain alleles due to their fitness advantage. This simulation combines both processes to study how allele frequencies evolve under these interacting forces.
# 
# ### Simulation Parameters:
# - **Fitness Values**:
#   - $W_{pp}$: Fitness of genotype $A_1A_1$.
#   - $W_{pq}$: Fitness of genotype $A_1A_2$.
#   - $W_{qq}$: Fitness of genotype $A_2A_2$.
# - **Initial Allele Frequencies**:
#   - $p$: Frequency of allele $A_1$.
#   - $q$: Frequency of allele $A_2$.
# - **Population Size**: The number of individuals simulated each generation.
# - **Number of Generations**: Total number of generations over which drift and selection are applied.
# 
# ### Key Results:
# The simulation demonstrates how selection and drift interact:
# 1. **Selection Effects**: If one allele has higher fitness, its frequency increases over generations.
# 2. **Drift Effects**: Fluctuations in allele frequencies occur, especially in smaller populations.
# 3. **Combined Effects**: Selection generally dominates in larger populations, while drift may obscure selection effects in smaller populations.
# 

# %%
def simulate_genetic_drift_with_selection(n_generations, f_A1=0.5, pop_size=100,
                                          Wpp=1.0, Wpq=1.0, Wqq=1.0):
    """
    Simulate genetic drift with selection over generations.

    Parameters:
    n_generations (int): Number of generations to simulate.
    f_A1 (float): Initial frequency of allele A1.
    pop_size (int): Population size.
    Wpp (float): Fitness of genotype A1A1.
    Wpq (float): Fitness of genotype A1A2.
    Wqq (float): Fitness of genotype A2A2.

    Returns:
    tuple: (generations, allele_frequencies)
    """
    generations = range(n_generations)
    allele_freqs = [f_A1]  # To store f_A1 over generations

    for _ in generations[1:]:
        # Simulate random mating
        genotypes = simulate_random_mating(pop_size, allele_freqs[-1])

        # Apply selection: Calculate new allele frequency based on fitness values
        numerator = (genotypes['A1A1'] * Wpp +
                     0.5 * genotypes['A1A2'] * Wpq +
                     0.5 * genotypes['A2A1'] * Wqq)
        denominator = (genotypes['A1A1'] * Wpp +
                       genotypes['A1A2'] * Wpq +
                       genotypes['A2A1'] * Wpq +
                       genotypes['A2A2'] * Wqq)
        f_A1 = numerator / denominator
        allele_freqs.append(f_A1)  # Store updated frequency

    return list(generations), allele_freqs

# Define multiple initial conditions
initial_conditions = [0.1, 0.3, 0.5, 0.7, 0.9]  # Initial frequencies of A1

# Define different sets of fitness values
fitness_sets = [
    {"Wpp": 0.1, "Wpq": 0.1, "Wqq": 0.0},  # Set 1: A2A2 has lowest fitness
    {"Wpp": 1.0, "Wpq": 0.8, "Wqq": 0.6},  # Set 2: A1A1 has highest fitness
    {"Wpp": 0.5, "Wpq": 1.0, "Wqq": 0.5},  # Set 3: Heterozygote advantage
]

# Number of generations to simulate
n_generations = 100

# Population size
pop_size = 1000

# Plot settings
plt.figure(figsize=(18, 12))  # Larger figure for multiple subplots
colors = ["blue", "green", "red", "yellow", "white"]  # Colors for different initial conditions

# Create a subplot for each set of fitness values
for idx, fitness in enumerate(fitness_sets):
    plt.subplot(2, 2, idx + 1)  # Arrange subplots in a 2x2 grid
    Wpp, Wpq, Wqq = fitness["Wpp"], fitness["Wpq"], fitness["Wqq"]

    # Simulate and plot for each initial condition
    for cond_idx, f_A1 in enumerate(initial_conditions):
        generations, allele_freqs = simulate_genetic_drift_with_selection(
            n_generations, f_A1, pop_size, Wpp, Wpq, Wqq
        )
        plt.plot(generations, allele_freqs, label=f'A1 (Initial f_A1 = {f_A1})',
                 color=colors[cond_idx], linestyle="-", linewidth=2)
        plt.plot(generations, [1.0 - p for p in allele_freqs], label=f'A2 (Initial f_A1 = {f_A1})',
                 color=colors[cond_idx], linestyle="--", linewidth=2)

    # Customize each subplot
    plt.title(f"Fitness: Wpp = {Wpp}, Wpq = {Wpq}, Wqq = {Wqq}", fontsize=14, pad=10)
    plt.xlabel('Generations', fontsize=12)
    plt.ylabel('Allele Frequency', fontsize=12)
    plt.legend(fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')  # Legend outside the plot
    plt.grid(True, linestyle='--', alpha=0.6)  # Add a grid for better readability

# Adjust layout and display the plot
plt.tight_layout()
plt.show()

# %% [markdown]
# ### Key Insights from the Output
# * Drift and Selection Interplay: Even with selection favoring $A_1$ (higher fitness $W_{pp}$ and $W_{pq}$), genetic drift causes fluctuations in allele frequencies. The larger the population, the smoother the trajectories of $p$ and $q$ due to reduced drift effects.
# 
# * Directional Selection: As $W_{qq} = 0.0$, allele $A_2$ (q) is strongly selected against. This drives $A_1$ (p) towards fixation over generations.
# The rate at which $A_1$ approaches fixation depends on the strength of selection (difference in fitness values) and population size.
# 
# * Stochastic Variability: With a moderate population size ($n=200$), drift introduces variability in $p$ and $q$, which would not be present in deterministic selection models. This is evident in the slight nonlinearity of the trajectories.
# 
# * Fixation and Extinction: Over time, one allele ($A_1$) becomes fixed (frequency = 1), while the other ($A_2$) is lost (frequency = 0). This aligns with the fitness landscape, where $A_2A_2$ has zero fitness.

# %% [markdown]
# ## Maintenance of a Harmful Allele
# 
# Some alleles produce harmful phenotypes when present in two copies (homozygous recessive), yet they persist in populations for extended periods. Why might this happen?
# 
# Two key mechanisms can explain this:
# 1. **Heterozygote Advantage**: The heterozygous genotype has a fitness advantage over both homozygotes, maintaining both alleles in the population.
# 2. **Mutation/Selection Balance**: A harmful allele can persist if the rate of mutation introducing the allele balances the selection pressure removing it.
# 
# #### Simulating Heterozygote Advantage
# Heterozygote advantage (or overdominance) occurs when the fitness of the heterozygote exceeds that of either homozygote. By definition, heterozygote advantage does not occur in cases of complete Mendelian dominance, where heterozygotes and homozygous dominant individuals have the same phenotype, and thus selection cannot distinguish between them.
# 
# Below, we simulate a scenario where the harmful allele persists due to heterozygote advantage.
# 
# ### Parameters:
# - **Fitness Costs**:
#   - $s$: Fitness reduction of the $A_2A_2$ genotype.
#   - $t$: Fitness reduction of the $A_1A_1$ genotype.
#   - $W_{pp} = 1 - t$: Fitness of $A_1A_1$.
#   - $W_{pq} = 1$: Fitness of $A_1A_2$.
#   - $W_{qq} = 1 - s$: Fitness of $A_2A_2$.
# - **Generations**: The simulation runs for 50 generations.
# - **Range of Initial $p$ Values**: A range of starting frequencies for allele $A_1$ to observe how equilibrium frequencies are reached.
# 
# ### Results:
# The simulation demonstrates that under heterozygote advantage:
# - Both alleles ($A_1$ and $A_2$) persist in the population.
# - The equilibrium frequency depends on the relative fitness costs ($s$ and $t$).
# 

# %%
def HW_with_selection(p, q, Wpp, Wpq, Wqq):
    """
    Calculate the new frequency of allele A1 under Hardy-Weinberg equilibrium with selection.

    Parameters:
    p (float): Frequency of allele A1.
    q (float): Frequency of allele A2.
    Wpp (float): Fitness of A1A1 homozygotes.
    Wpq (float): Fitness of A1A2 heterozygotes.
    Wqq (float): Fitness of A2A2 homozygotes.

    Returns:
    float: New frequency of allele A1.
    """
    numerator = p**2 * Wpp + p * q * Wpq
    denominator = p**2 * Wpp + 2 * p * q * Wpq + q**2 * Wqq
    return numerator / denominator

# Define simulation parameters
generations = 50  # Number of generations

# Define multiple fitness parameter sets
fitness_sets = [
    {"s": 0.3, "t": 0.7},  # Set 1: s = 0.3, t = 0.7
    {"s": 0.1, "t": 0.5},  # Set 2: s = 0.1, t = 0.5
    {"s": 0.5, "t": 0.9},  # Set 3: s = 0.5, t = 0.9
]

# Set up the plot
plt.figure(figsize=(18, 12))  # Larger figure for multiple subplots

# Create a subplot for each set of fitness parameters
for idx, fitness in enumerate(fitness_sets):
    plt.subplot(2, 2, idx + 1)  # Arrange subplots in a 2x2 grid
    s, t = fitness["s"], fitness["t"]
    Wpp = 1.0 - t  # Fitness of A1A1
    Wpq = 1.0      # Fitness of A1A2
    Wqq = 1.0 - s  # Fitness of A2A2

    # Simulate for a range of initial allele frequencies (p values)
    for p in np.linspace(0.0, 1.0, 50):  # Range of initial frequencies for A1
        q = 1.0 - p  # Frequency of A2
        starting_p = p  # Store the initial value of p for color mapping

        # Store results for plotting
        xs = []  # Generations
        qs = []  # Frequencies of A2

        # Simulate over generations
        for g in range(generations):
            new_p = HW_with_selection(p, q, Wpp, Wpq, Wqq)  # Update p with selection
            xs.append(g)       # Append generation
            q = 1.0 - new_p    # Update q
            qs.append(q)       # Append new q
            p = new_p          # Update p for next generation

        # Plot frequency of A2 (q) over generations
        plt.plot(xs, qs, '-', color=plt.cm.Spectral(starting_p),
                 linewidth=0.5)

    # Add final plot elements
    plt.ylabel('Allele Frequency of Allele A2', fontsize=12)
    plt.xlabel('Generations', fontsize=12)
    plt.title(f"Fitness: s = {s}, t = {t}", fontsize=14, pad=10)
    plt.grid(True, linestyle='-', alpha=0.6)  # Add a grid for better readability

# Adjust layout and display the plot
plt.tight_layout()
plt.show()

# %% [markdown]
# ### Key Insights from the Simulation
# * Persistence of Both Alleles: Both $A_1$ and $A_2$ persist in the population due to heterozygote advantage. The exact equilibrium frequencies depend on the relative fitness of the genotypes.
# * Equilibrium Frequency: For $s = 0.6$ and $t = 0.4$, the harmful allele ($A_2$) stabilizes at an intermediate frequency. This equilibrium occurs where the fitness advantage of heterozygotes balances the fitness disadvantages of the homozygotes.
# * Initial Conditions: Regardless of the initial frequency of $A_1$ ($p_0$), the population converges to the same equilibrium frequencies for $A_1$ and $A_2$, showing the stable nature of heterozygote advantage.
# * Harmful Allele Maintenance: Even though $A_2A_2$ individuals suffer a significant fitness cost, the allele persists due to the higher fitness of $A_1A_2$ individuals.

# %% [markdown]
# ## Heterozygote Disadvantage
# 
# In cases of heterozygote disadvantage (or underdominance), the fitness of heterozygotes ($A_1A_2$) is lower than either homozygote ($A_1A_1$ or $A_2A_2$). This creates a unique dynamic where the equilibrium frequency of the harmful allele ($A_2$) depends strongly on its **initial frequency**.
# 
# Here’s why:
# 1. If $A_1A_2$ heterozygotes have low fitness and the population starts with all $A_1$, it is challenging for $A_2$ to increase in frequency. This is because the transition involves passing through a range of allele frequencies dominated by unfit heterozygotes before reaching a high enough frequency of $A_2$ to form $A_2A_2$ homozygotes (which may have higher fitness).
# 2. Consequently, under heterozygote disadvantage, allele $A_2$ can persist or be eliminated depending on:
#    - Its **initial frequency** ($f(A_2)$ at generation 0).
#    - The **relative fitness costs** of $A_1A_2$ heterozygotes versus $A_2A_2$ homozygotes.
# 
# The plot below demonstrates how the starting frequency and fitness parameters shape the persistence or elimination of the $A_2$ allele.
# 
# ### Parameters:
# - **Fitness Costs**:
#   - $h$: Fitness cost of $A_1A_2$ heterozygotes.
#   - $s$: Fitness cost of $A_2A_2$ homozygotes.
#   - $W_{pp} = 1.0$: Fitness of $A_1A_1$.
#   - $W_{pq} = 1.0 - h$: Fitness of $A_1A_2$.
#   - $W_{qq} = 1.0 - s$: Fitness of $A_2A_2$.
# - **Generations**: The simulation runs for 300 generations.
# - **Range of Initial $p$ Values**: A range of starting frequencies for allele $A_1$ to examine its impact on equilibrium outcomes.
# 
# ### Results:
# The simulation reveals how $A_2$ behaves under heterozygote disadvantage:
# - Allele $A_2$ persists or declines depending on its **initial frequency** and the fitness disadvantages of heterozygotes.

# %%
# Define simulation parameters
generations = 300  # Number of generations

# Define multiple fitness parameter sets
fitness_sets = [
    {"h": 0.06, "s": 0.03},  # Set 1
    {"h": 0.1, "s": 0.05},   # Set 2
    {"h": 0.2, "s": 0.1},    # Set 3
]

# Set up the plot
plt.figure(figsize=(18, 12))  # Larger figure for multiple subplots

# Create a subplot for each set of fitness parameters
for idx, fitness in enumerate(fitness_sets):
    plt.subplot(2, 2, idx + 1)  # Arrange subplots in a 2x2 grid
    h, s = fitness["h"], fitness["s"]
    Wpp = 1.0       # Fitness of A1A1
    Wpq = 1.0 - h   # Fitness of A1A2
    Wqq = 1.0 - s   # Fitness of A2A2

    # Simulate for a range of initial allele frequencies (p values)
    for p in np.linspace(0.0, 1.0, 30):  # Range of initial frequencies for A1
        q = 1.0 - p  # Frequency of A2
        starting_p = p  # Store the initial value of p for color mapping

        # Store results for plotting
        xs = []  # Generations
        qs = []  # Frequencies of A2

        # Simulate over generations
        for g in range(generations):
            new_p = HW_with_selection(p, q, Wpp, Wpq, Wqq)  # Update p with selection
            xs.append(g)       # Append generation
            q = 1.0 - new_p    # Update q
            qs.append(q)       # Append new q
            p = new_p          # Update p for next generation

        # Plot frequency of A2 (q) over generations
        plt.plot(xs, qs, '-', color=plt.cm.Spectral(starting_p), linewidth=0.5)

    # Add final plot elements
    plt.ylabel('Allele Frequency of A2', fontsize=12)
    plt.xlabel('Generations', fontsize=12)
    plt.title(f"Fitness: h = {h}, s = {s}", fontsize=14, pad=10)
    plt.grid(True, linestyle='-', alpha=0.6)  # Add a grid for better readability

# Adjust layout and display the plot
plt.tight_layout()
plt.show()

# %% [markdown]
# ### Key Insights from the Simulation
# * Dependence on Starting Frequencies: The harmful allele $A_2$ can persist or disappear depending on its initial frequency. If $A_2$ starts at a low frequency, it struggles to increase due to the unfit heterozygote ($A_1A_2$).
# If $A_2$ starts at a high enough frequency, it may overcome the heterozygote barrier and stabilize at an equilibrium frequency.
# * Threshold Frequency: There exists a threshold starting frequency above which $A_2$ persists and below which it is eliminated.
# * Fitness Costs Matter: The higher the fitness cost of heterozygotes relative to $A_2A_2$ homozygotes, the higher the threshold frequency for $A_2$ to persist.
# * Generational Dynamics: The transition dynamics are slow when $A_2$ starts near the threshold frequency, as the population spends significant time in heterozygote-dominated states.

# %% [markdown]
# ## Mutation-Selection Balance
# 
# Mutation-selection balance describes a scenario where harmful alleles persist in a population due to a balance between natural selection (which removes them) and mutation (which introduces them).
# 
# #### Key Parameters:
# - **Fitness Costs**:
#   - $W_{pp} = 1.0$: Fitness of $A_1A_1$.
#   - $W_{pq} = 1.0$: Fitness of $A_1A_2$.
#   - $W_{qq} = 1.0 - s$: Fitness of $A_2A_2$ ($s$ is the fitness cost of the $A_2A_2$ genotype).
# - **Mutation Rate**:
#   - `mutation_rate_p_q`: The probability of allele $A_1$ mutating to $A_2$ in one generation. (We assume no back mutation from $A_2$ to $A_1$.)
# - **Generations**: The simulation runs for 300 generations.
# - **Range of Starting $p$ Values**: The initial frequency of $A_1$ is varied across simulations.
# 
# ### Results:
# This simulation demonstrates how mutation and selection act in opposing directions:
# 1. **Selection** reduces the frequency of the harmful allele $A_2$ (by reducing the fitness of $A_2A_2$ homozygotes).
# 2. **Mutation** introduces $A_2$ alleles into the population, opposing selection.
# 3. Over time, the population stabilizes at a balance point where the rate of mutation equals the rate of selection against $A_2$.
# 
# ### Observations:
# - When the mutation rate is small, $A_2$ can be nearly eliminated by selection.
# - At higher mutation rates, $A_2$ persists at higher frequencies, even with a significant fitness cost.
# 
# The plot below illustrates how allele frequencies change over generations for different initial values of $A_1$ ($p$).
# 

# %%
# Define simulation parameters
generations = 300  # Number of generations
mutation_rate_p_q = 0.001  # Mutation rate from A1 to A2

# Define multiple fitness parameter sets
fitness_sets = [
    {"s": 0.05},  # Set 1
    {"s": 0.1},   # Set 2
    {"s": 0.2},   # Set 3
]

# Set up the plot
plt.figure(figsize=(18, 12))  # Larger figure for multiple subplots

# Create a subplot for each set of fitness parameters
for idx, fitness in enumerate(fitness_sets):
    plt.subplot(2, 2, idx + 1)  # Arrange subplots in a 2x2 grid
    s = fitness["s"]
    Wpp = 1.0       # Fitness of A1A1
    Wpq = 1.0       # Fitness of A1A2
    Wqq = 1.0 - s   # Fitness of A2A2

    # Simulate for a range of initial allele frequencies (p values)
    for p in np.linspace(0.0, 1.0, 30):  # Range of initial frequencies for A1
        q = 1.0 - p  # Frequency of A2
        starting_p = p  # Store the initial value of p for color mapping

        # Store results for plotting
        xs = []  # Generations
        qs = []  # Frequencies of A2

        # Simulate over generations
        for g in range(generations):
            new_p = HW_with_selection(p, q, Wpp, Wpq, Wqq)  # Update p with selection
            new_p -= new_p * mutation_rate_p_q  # Apply mutation from A1 to A2
            xs.append(g)       # Append generation
            q = 1.0 - new_p    # Update q
            qs.append(q)       # Append new q
            p = new_p          # Update p for next generation

        # Plot frequency of A2 (q) over generations
        plt.plot(xs, qs, '-', color=plt.cm.Spectral(starting_p), linewidth=0.5)

    # Add final plot elements
    plt.ylabel('Allele Frequency of A2', fontsize=12)
    plt.xlabel('Generations', fontsize=12)
    plt.title(f"Fitness: s = {s}", fontsize=14, pad=10)
    plt.grid(True, linestyle='-', alpha=0.6)  # Add a grid for better readability

# Adjust layout and display the plot
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Key Insights from the Simulation
# * Balance Point: The harmful allele $A_2$ stabilizes at a frequency where the rate of mutation introducing $A_2$ equals the rate of selection removing it.
# * Effect of Mutation Rate: Low mutation rates result in very low equilibrium frequencies of $A_2$. Higher mutation rates allow $A_2$ to persist at higher frequencies, even if it is strongly selected against.
# * Starting Frequencies: Regardless of the initial frequency of $A_2$, the population converges to the same equilibrium frequency dictated by the mutation-selection balance.
# * Generational Dynamics: The time to reach equilibrium depends on both the initial frequency of $A_2$ and the strength of selection and mutation.
# 

# %% [markdown]
# ### Testing Nearly Neutral Theory
# 
# Nearly neutral theory suggests that alleles subject to weak selection ($s \ll \frac{1}{N_e}$, where $N_e$ is the effective population size) behave as if they are effectively neutral. This means that the fixation rates of these alleles are primarily driven by genetic drift, rather than selection.
# 
# #### Key Concepts:
# - **Effectively Neutral Alleles**:
#   - When the selection coefficient ($s$) is much smaller than $\frac{1}{N_e}$, selection is too weak to influence allele frequencies significantly.
#   - These alleles are subject to drift, and their fixation or loss is random.
#   
# - **Testing the Prediction**:
#   - The simulation varies population sizes ($N_e$) and tracks the fixation rates of an allele with a weak selection coefficient ($s = \frac{1}{10}$).
#   - It categorizes populations where $s < \frac{1}{N_e}$ as effectively neutral and those with $s \geq \frac{1}{N_e}$ as subject to selection.
#   
# - **Fixation Rate**:
#   - The proportion of simulations where the allele becomes fixed (reaches a frequency of 1.0).
# 
# #### Observations:
# 1. For small populations (where $N_e$ is small), alleles with $s \ll \frac{1}{N_e}$ behave neutrally, with fixation rates close to the expectation for neutral alleles.
# 2. As population size increases, selection becomes stronger relative to drift, and fixation rates deviate from neutral predictions.
# 3. Fixation rates in neutral vs. selected cases are plotted against population size, illustrating the transition from drift-dominated to selection-dominated regimes.
# 
# #### Simulation Results:
# The plot below shows:
# - **Blue Line**: Fixation rates for populations where $s \ll \frac{1}{N_e}$ (neutral).
# - **Red Line**: Fixation rates for populations where $s \geq \frac{1}{N_e}$ (selected).
# 
# These results demonstrate how the nearly neutral theory operates under varying population sizes and selection intensities.

# %%
def simulate_genetic_drift_with_selection(n_generations, f_A1=0.5, pop_size=100, Wpp=1.0, Wpq=1.0, Wqq=0.99):
    """Simulates allele frequency change over generations under weak selection."""
    allele_freqs = []
    for generation in range(n_generations):
        genotypes = simulate_random_mating(pop_size, f_A1)
        numerator = (
            genotypes["A1A1"] * Wpp
            + 0.5 * genotypes["A1A2"] * Wpq
            + 0.5 * genotypes["A2A1"] * Wpq
        )
        denominator = (
            genotypes["A1A1"] * Wpp
            + genotypes["A1A2"] * Wpq
            + genotypes["A2A2"] * max(Wqq, 1e-5)  # Prevent Wqq from being too small
        )
        
        # Prevent division by zero
        f_A1 = numerator / max(denominator, 1e-9)
        allele_freqs.append(f_A1)
    return allele_freqs

# %%
# Simulation parameters
plt.figure(figsize=(10, 6))
neutral_xs, neutral_ys = [], []
selected_xs, selected_ys = [], []
s = 1.0 / 10.0  # Selection coefficient
population_sizes = np.logspace(0.3, 3, num=20, dtype=int)  # Log-spaced population sizes

for pop_size in population_sizes:
    print(f"Simulating for population size: {pop_size}")
    losses, fixation = 0, 0

    for _ in range(1000):  # 1000 replicate simulations
        generations = 200
        Wpp, Wpq, Wqq = 1.0, 1.0, 1.0 - s
        f_A1s = simulate_genetic_drift_with_selection(generations, f_A1=0.5, pop_size=pop_size, Wpp=Wpp, Wpq=Wpq, Wqq=Wqq)

        if f_A1s[-1] == 0.0:
            losses += 1
        elif f_A1s[-1] == 1.0:
            fixation += 1

    fixation_freq = fixation / max(losses + fixation, 1)  # Avoid division by zero

    if s < 1 / pop_size:
        neutral_xs.append(pop_size)
        neutral_ys.append(fixation_freq)
    else:
        selected_xs.append(pop_size)
        selected_ys.append(fixation_freq)

# Plot results
plt.plot(neutral_xs, neutral_ys, "bo-", label="Neutral (s ≪ 1/Nₑ)", linewidth=2)
plt.plot(selected_xs, selected_ys, "ro-", label="Selected (s ≥ 1/Nₑ)", linewidth=2)

# Add theoretical neutral expectation line
plt.axhline(y=0.5, color="gray", linestyle="dashed", label="Neutral Expectation f(A₁) = 0.5")

# Improve readability
plt.xscale("log")  # Log scale for population size
plt.ylabel("Fixation Probability of A1")
plt.xlabel("Effective Population Size (Nₑ)")
plt.legend()
plt.title("Fixation Probability vs. Population Size")
plt.show()