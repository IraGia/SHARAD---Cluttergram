# SHARAD-Cluttergram
A Matlab function that generates cluttergrams for SHARAD based on x,y,z topography data and SHARAD's orbit. 


Ground-Penetrating Radar (GPR) is an ideal tool for imaging the subsurface of planetary bodies due to its high resolution, continuous measurement capability, and ability to operate in dry, low-conductivity environments. As a result, both satellite and in-situ radars are becoming mainstream geophysical approaches for imaging terrestrial planets, with four active radar sounders (SHARAD, MARSIS, Tianwen-1, SELENE) [1] and three in-situ GPR systems (Chang'E-4, Tianwen-1, Perseverance) [2]. Data from satellite orbiters have been successfully used for mapping subglacial water on Mars and detecting Lunar and Martian caves [3]. Caves such as lava tubes, piping caves, and sub-ice volcanic caves are unique geological structures of great importance in planetary science [3]. Unlike exposed surfaces, caves are shielded environments protected from physicochemical decay and isolated from strong fluxes of ultraviolet cosmic and solar ionizing radiation [3]. For these reasons, caves are ideal places for searching for evidence of current and past life and for potential human habitation and planetary bases [3].

Detecting and mapping planetary subsurface structures requires processing and interpreting a large volume of radar data in a semi-automatic manner. Despite the plethora of available data from radar orbiters, interpretation remains challenging because radargrams are often corrupted with surface reflections that mask signatures from subsurface targets [4]. Therefore, the first and most important processing step is generating a cluttergram that simulates the surface reflections expected in the area of interest. These reflections are then filtered out, enhancing and isolating the radar signatures of subsurface targets [4]. In this function we use the approahc described in [5], with the addition of incorporating the directivity pattern of a dipole in the calculations. 



[1] Xu, Y., Cummer, S. A., Farrell, W. M., “Application of an orbital radar sounder model to detecting Martian polar subsurface features,” Journal of Geophysical Research, vol. 111, pp. E06S17.

[2] Hamran, S.-E., Paige, D. A., Amundsen, H. E. F., et al., "Radar Imager for Mars’ Subsurface Experiment—RIMFAX," Space Science Reviews, vol. 216, 128, 2020.

[3] Sam, L., Bhardwaj, A., Singh, S., Martin-Torres, F. J., et al., “Small Lava Caves as Possible Exploratory Targets on Mars: Analogies Drawn from UAV Imaging of an Icelandic Lava Field,” Remote Sensing, vol. 12, pp. 1970, 2020.

[4] Ilyushin, Y. A., Orosei, R., Witasse, O., Sánchez-Cano, B., "CLUSIM: A synthetic aperture radar clutter simulator for planetary exploration," Radio Science, vol. 52, pp. 1200–1213, 2017.

[5] Choudhary, P., Holt, J. W., Kempf, S. D., "Surface Clutter and Echo Location Analysis for the Interpretation of SHARAD Data From Mars," IEEE Geoscience and Remote Sensing Letters, vol. 13, pp. 1285–1289, 2016.
