## The effect of single versus successive warm summer on intertidal communities

This experiment set out to test how temperature shapes the composition of intertidal barnacle bed communities, asking the question: how do single vs. successive warm summers affect this community? We used a passive warming manipulation of black and white settlement tiles to generate differences in substratum temperatures and tracked communities over two years. After the first year, the treatment of half of each group (warm and cool) were swapped to manipulate thermal stress through time.

We expected that:

1.  Barnacle bed communities that are exposed to hotter temperatures during summer, even for a single year, will have lower organism abundances and diversity (species richness, Shannon diversity) than those that are exposed to ambient/cooler conditions during the same period.

2.  There will be an interactive effect between the temperature treatments of the first and second summer on the same response metrics. Previously 'cool' communities, since they have more established, larger barnacle beds with a more diverse array of microhabitats and thermal refugia, will be less perturbed by warming than previously 'warm' communities that have less structurally complex biogenic habitat. 

# Data types (in raw_data and clean_data folders)
- Temperature: hourly temperature collected for experimental tiles and adjacent bedrock
- Biological survey data: communities were surveyed monthly during summer and every two months during winter to identify and count species on tiles (or record percent cover, as in the case of algae)
- Epifaunal community data: tiles were destructively sampled the fall and winter of year two to look at cryptic diversity within treatments.
- Tide and daylight hours data: these are pulled from Canadian government sources for analyses of temperature.

For a full data dictionary, please see the project_summary.pdf file within the metadata folder.

# Analyses conducted (in scripts folders)
1) Temperature differences between treatments, in terms of maximum daily temperature and mean temperature -- 03_temperature.R 
2) Barnacle abundance (recruits and adults) differences between treatments -- 04_barnacles.R
3) Grazer (littorine snail and limpet) abundance differences between treatments -- 05_grazers.R
4) Algal cover differences between treatments -- 06_algae.R
5) Alpha diversity differences between treatments -- 07_diversity.R
6) Community composition and beta diversity between treatments -- 07_diversity.R
7) Trajectory of community structure through time -- 08_cta.R
