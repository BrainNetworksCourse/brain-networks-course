In this exercise you will become acquainted with how to visualize images and atlases using the FSLeyes viewer.  

To being, start the virtual machine by running "vagrant up" in the appropriate directory.  Once the VM is running, start a terminal window by clicking the bottom left icon and select Accessories->LXTerminal.  Start FSLeyes by typing "fsleyes" in the terminal.

In Fsleyes:

- Load template image
- File->Add standard
    - Select “MNI152_T1_2mm_brain.nii.gz”
- Open atlas panel
    - Settings->OrthoView1->Atlas Panel
- Find Broca’s area
     - Click “atlas search” in the atlas tools panel
     - Type “Broca” in the search window, and then click on “Juelich Histological Atlas” which should be highlighted
    - Click the check box next to “GM Broca’s Area BA44 L” to present a probabilistic map of Brodmann’s area 44
    - Select the Broca’s image from the overlay list (to the left), and adjust the min and max at the top to range from 20 to 100 and the colormap to Cool
    - Go to voxel location [71, 69, 51] using the Location panel - what is the probability that this voxel is in BA44?
- Identify white matter tracts
    - Click “Atlas information” in the atlas tools panel, and click the check box next to JHU White Matter Tractography Atlas
    - Scroll down to the listing for the White matter atlas to the right of the atlas listing, and click “show/hide” to show the atlas overlaid on the brain
    - Surf around in the region of BA44 - do you see any white matter tracts that appear to terminate in the vicinity of BA44?
- Load activation image
    - File->add from file
        - Select data/neurovault/nv_304.nii.gz
    - Change color slider from Grayscale to Red-Yellow
    - Set min and max values to 3 and 10 respectively
    - Does the activation area overlap with Broca’s area?  You can click on the nv304 image in the Overlay list and then vary its opacity using the slider in the top toolbar.
    - Does the white matter tract nearest Broca’s area also connect to any other areas of activation? To find this, first turn off the JHU-tracts overlay by clicking the eye next to it in the overlays list.  Then use the atlas search tool to find the specific tract from the JHU tractography atlas, and click its box to create an overlay (setting it to 50-100 range and coloring it green to stand out).  Then turn on the Harvard-Oxford Cortical Structural atlas to identify the posterior region closest to the tract. Record the MNI152 coordinates of the peak activation in this posterior area.
- Using neurosynth
    - Go to http://neurosynth.org/locations/ and enter those coordinates
    - First examine the maps section, which shows a map of other voxels that are correlated with this voxel in resting fMRI data.  Does the map include Broca’s area as well?
    - Then examine the Associations tab - what are the terms that are most strongly associated with this image?
    - Based on these data, what do you think the task is that generated these activation data?
- White matter anatomy exercises:
    - Find each of the following white matter tracts, and describe which regions they appear to connect:
        - Forceps major
        - Inferior fronto-occipital fasciculus
        - Uncinate fasciculus
        - Forceps minor
        - Inferior longitudinal fasciculus
