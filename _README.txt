Bastiaen Boekelo, March 2022
Short description of Script content


Part 1: First results -> Scripts "00" to "06":
- Convert original input in easier to handle formats (e.g. same naming conventions, or multiband to singleband)
- Create different vegetation indices, with different spatial operations (weighing, removing shadows and shines)
- Perform zonal statistics, write as csv (mostly written in Data/2_Intermediate, 
- Inspect some of the results (mostly written in Data/3_Output)

These operations can be considered the base and 'pre-work' for the other scripts.


Part 2: Creating spatial variables -> Scripts "07" to "15"


Part 3: Key scripts creating key graphs -> "17" to "18"


17_Combine_and_clean_data.R:
	- Creates 1 .csv file combining multiple variables created earlier in the process
	- Shows what calculation / conversions have been performed to derive the nutrient variables 

18_Create_corplots.R:
	- Creating key graphs
