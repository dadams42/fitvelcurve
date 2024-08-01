# fitvelcurve
Overview: This work is the result of a one-week group project from the 2024 Code/Astro workshop hosted at the CIERA Center, Northwestern University, on software engineering for astrophysics.

Title: Open-source package for derived dark matter density profiles from rotation curves.

Description: We fit an NFW profile [1] using the nested sampling package dynesty [2] to synthetic (mock) rotation curve data and perform an end-to-end test to validate against the underlying profile for the halo. 

Acknowledgments: We thank the team at the Code/Astro workshop for their invaluable help in making this possible, including Sarah Blunt, Jason Wang, and Suchitra Narayanan, who oversaw this project. 

Sources: For information on dynesty, see: https://dynesty.readthedocs.io/en/stable/ 

References:

[1]: Navarro, Julio F., Carlos S. Frenk, and Simon DM White. "A universal density profile from hierarchical clustering." The Astrophysical Journal 490.2 (1997): 493.

[2]: Joshua S Speagle, "dynesty: a dynamic nested sampling package for estimating Bayesian posteriors and evidences." Monthly Notices of the Royal Astronomical Society, Volume 493, Issue 3, April 2020, Pages 3132–3158, https://doi.org/10.1093/mnras/staa278
