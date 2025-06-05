# Fiber-to-Fiber Spacing Algorithm
A Python tool for interactive analysis of fiber bundle images, including rotation correction, geometric fitting, and measurement of spacing and alignment.

![Alt text](FullObjBundle.bmp)

project/

├── rotated/                 # Processed outputs (centroid & ellipse data)

├── plots/                   # Saved histograms and visualization images

├── Centroid.py              # Main script

└── README.md                # You're here!

## Overview

The script performs the following main functions:
1. Loads fiber bundle images (BMP format)
2. Allows interactive point selection for scaling and alignment
3. Detects fiber positions in a hexagonal pattern
4. Fits ellipses to each fiber core
5. Analyzes fiber sizes and grades the bundle quality
6. Generates visualizations and output files

## Dependencies

- Python 3.9
- Required libraries:
  - numpy
  - Operating System (os)
  - Pillow (PIL)
  - math
  - matplotlib
  - scipy 

## Defined Functions

### Key Functions
Interactive Functions
two_clicks(): Handles two right-clicks for point selection

one_clicks(): Handles single right-click for point selection

two_clicks_y_scaling(): Handles y-axis scaling point selection

### Geometric Functions
RotationM(): Creates 2D rotation matrix

GetRotatedMajorMinor(): Rotates ellipse axes

drawEllipse(): Draws ellipse given parameters

### Analysis Functions
EllipseFitter(): Fits ellipse to fiber core image

Centroid calculation using center_of_mass

Histogram analysis with peak detection

## Algorithm
Core Detection
Preprocessing

- Image rotation

- Dual-axis scaling

Grid Generation

- row_spacing = y_scale * diameter * (√3/2)

Centroid Refinement

- Initial center_of_mass

- Radial thresholding

- Secondary centroid

- Boundary Fitting

Ellipse parameters:

- angle, major, minor = EllipseFitter(thimg)

## 1. How to Run

1.Place your fiber bundle images (BMP format) in a directory (e.g., C:\Users\user\Desktop\Masters\Fiber_Spacing\)

2.  - Update the path variable in the script to point to your image directory.
    - In the script look for the line with equation, rough_scale = float(x2 - x1) / (20 * 239), change the value 20 to the number of fibres in the middle row of you fiber bundle minus 1, change the 239 value to the diameter of a single fibre.
    - halfwidth = int(rough_scale * 119.5), change 119.5 to the raius (diameter/2) of your fibre.
    - rsq = (119.5 * rough_scale)**2, change 119.5 to the raius (diameter/2) of your fibre.
    - y_rough_scale = float(y_y2 - y_y1) / (18 * (239*np.sqrt(3)/2)), change 18 to the number of rows in your fiber bundle minus 1.
    - threshold_area = np.pi * (((y_rough_scale + rough_scale) / 2) * 50)**2, change 50 to your fiber core radius.

3.Go the Centroid.py file directory in the python terminal/cmd and run the file using #python Centroid.py

4.Run the script by pressing enter on your computer:

![Alt text](Images/Screenshot(16).png)

A folder named rototed will be created, along with '_.dat'and '-centroid.txt' files.

## 2. Interactive steps
### Step 1: X-Axis Calibration and Image Rotation
A window will pop-up with an image prompting you to select two points on the image. You can use the zoom icon indicated in the below image using a blue arrow to zoom in on a single fibre in order to acurately select a point, then you will left-click the left facing arrow to go back to the full image and zoom in again on the next fibre and select the next point.
Right click the centre of the first fibre in the middle row, then right click the centre of the last fibre on the middle row of the fibre bundle. The image will close by itself after you select the required number of points.

![Alt text](\Images/Screenshot(15)1.png)
![Alt text](\Images/Screenshot(11).png)

The coordinates of the two points that you right-clicked on will be printed out in the terminal. These two points will be used in the code to rotate/align the image such that the the fibre bundle is not tilted.

### Step 2: Y-Axis Calibration
![Alt text](\Images/Screenshot(15)2.png)

Now a second pop-up with the rotated fibre bundle will appear. You will also zoom on the target fibres as explained previously. You have to right click the first fiber on the last row of the fiber bundle (or right click the centre of one of the fibres on the last row of the fibre bundle that will be colinear with the centre of the first fibre on the first row.), then you right click the centre of the first fibre on the first row. The image will close by itself after you select the required number of points.

![Alt text](\Images/Screenshot(15)3.png)
![Alt text](\Images/Screenshot(12).png)

The coordinates of the two points that you right-clicked on will be printed out in the terminal. Then the values of the x and y -scaling will the displayed on the terminal.

![Alt text](\Images/Screenshot(15)4.png)

### Step 3: Bundle Center
Now a third pop-up with the rotated fibre bundle will show. Zoom into the middle fibre/centre point of the fibre bundle. You have to right click the centre of the middle fibre on the middle row of the fibre bundle (or simply the centre point of the entire fibre bundle if the rows and/or number of fibers in the middle row are even). Once the pont is selected the image will close.

![Alt text](\Images/Screenshot(15)5.png)
![Alt text](\Images/Screenshot(13).png)

The centre coordinates of the middle fibre that you right-clicked will be printed on the terminal. 

### Step 4: Parameters
Next, you will be required to entre the diameter of a single fibre in micronmetres, followed by the number of rows of the fibre bundle, followed by the number of fibres in the middle row of the fiber bundle, you will have to press enter after inputing each of the values to proceed to the next step.

![Alt text](\Images/Screenshot(15).png)

### 2.1 Generate Grid Algorithm

An algorithm that calculates the centres of each of the fibres in the fibre bundle will use the coordinates from the right-click events, to generate a fibre grid. Then a second algorithm will generate a fibre grid using the centre_of_mass function to get the actual centres of the fibres on the image.

A histogram of the Number of pixels versus the light intensity will be generated for a subimg of each of the fibres and each histogram will be saved to the plots folder which is in the same folder as the main python script.  

![Alt text](\plots/histogram_75.png)

The value of light intensity between two of the highest peaks (D_p) will be used to threshold the subimg and decide whether to fit an ellipse/octagon or not.

![Alt text](\Images/screenshot.png)

Then the graph of the generated fibre grid based on the scaling factors is shown. You have to close the image for the next one to show.

![Alt text](\Images/Figure_2.png)

The graph of the generated fibre grid based of the center_of_mass function is shown. You have to close the image for the next one to show.

![Alt text](\Images/Figure_3.png)

An image of the generated fibre grids based of the center_of_mass function plotted on the image of the fibres is shown.

![Alt text](\Images/Figure_4.png)

An image of the generated fibre grids based on the scaling factors and generated fibre grids based of the center_of_mass function plotted on the image of the fibres is shown. You have to close the image for the next one to show.

![Alt text](\Images/Figure_5.png)

An image of the centroids of the fibres plotted on the image of the fibre bundle is shown. You have to close the image for the next one to show.

![Alt text](\Images/Figure_1.png)

Then the following image which can be optional, is the image of the fitted ellipses (change to octagons if necessary) plotted on the edges of the fibers on the fibre bundle image. This is only necessary if you are not sure that the centroids calculated/plotted for the fibres are accurate. You close the image when you done viewing it.

![Alt text](\Images/Figure_6.png)

### 2.2 Grading
Once you close the image of the fitted ellipses/octagons the holes/fibres of the image will be graded accoeding to the following scheme:

Fiber cores are graded based on minor axis size.

Grading scale:

| Grade | Size Range (µm) | Quality |
|----------|----------|----------|
| Grade 0  | < 120 µm  | Poor/Failed  |
| Grade 1  | 120–130 µm  | Marginal  |
| Grade 2  | 130–140 µm  | Good  |
| Grade 3  | 140 µm  | Excellent  |

![Alt text](\Images/Screenshot(14).png)

### 2.3 Output Files

project/

├── rotated/               # Processed images and data

│   ├── *_centroid.txt     # Coordinate files

│   └── *_params.txt       # Ellipse parameters

├── plots/                 # Analysis visuals

│   ├── histogram_*.png    # Intensity profiles

│   └── fitting_*.png      # Geometric fits

└── Centroid.py            # Main analysis script

## Troubleshooting
### Common Issues:

Incorrect scaling: Ensure calibration points are on exact fiber centers

Missing fibers: Verify input matches actual bundle geometry

Poor fitting: Check image contrast and focus
  - radius_threshold = (fiber_diameter * y_rough_scale) / 2 * 0.9
  - all_peaks, _ = find_peaks(hist_smooth, height=np.max(hist_smooth)*0.01, distance=100), peak lines missing
  - D_p = 550  # Default value
  - thimg = (thresholded_subimg > (0.40 * D_p))

### Debug Tips:

Review terminal output for coordinates

Examine histogram plots for clear bimodal distribution
