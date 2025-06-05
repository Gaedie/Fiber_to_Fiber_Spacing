import numpy as np
import os
from PIL import Image
import math
import matplotlib.pyplot as plt
from scipy import optimize, ndimage
from scipy.ndimage import center_of_mass
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d  # For smoothing
import cv2
from skimage.measure import label, regionprops

# Global variables to store click coordinates
x1, y1, x2, y2 = -1, -1, -1, -1
def two_clicks(event):
    """Handles two right-clicks to select points."""
    global x1, y1, x2, y2
    if event.button == 3:  # Right-click detected
        if x1 == -1:  # First click
            x1, y1 = event.xdata, event.ydata
            print(f"First click at: ({x1}, {y1})")
        else:  # Second click
            x2, y2 = event.xdata, event.ydata
            print(f"Second click at: ({x2}, {y2})")
            plt.close()

def one_clicks(event):
    """Handles one right-click to select a point."""
    global x1, y1
    if event.button == 3:  # Right-click detected
        x1, y1 = event.xdata, event.ydata
        print(f"Clicked point: ({x1}, {y1})")
        plt.close()

def RotationM(theta):
    """Creates a 2D rotation matrix."""
    return np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])

def GetRotatedMajorMinor(major, minor, angle, woundAngle):
    """Rotates the major and minor axes of an ellipse."""
    th = (woundAngle - angle) * np.pi / 180
    R = RotationM(th)
    Rinv = RotationM(-th)
    I = np.sqrt(np.dot(Rinv, np.dot([[major**2, 0], [0, minor**2]], R)))
    return I[0][0], I[1][1]

def GetRotatedBoundingBox(width_oct, height_oct, angle, woundAngle):
    """
    Rotates the bounding box dimensions (width and height) from OctagonFitter by a given reference angle.
    Parameters:
        width_oct (float): The width of the octagon's bounding box.
        height_oct (float): The height of the octagon's bounding box.
        angle (float): The measured orientation (in degrees) of the shape from OctagonFitter.
        woundAngle (float): The reference (or “wound”) angle you wish to align to.   
    Returns:
        rotated_width (float): The rotated width.
        rotated_height (float): The rotated height.
    """
    # Calculate the rotation angle in radians (difference between desired and measured angles)
    th = (woundAngle - angle) * np.pi / 180.0
    # Build the 2x2 rotation matrices (assumes you have a function RotationM defined)
    R = RotationM(th)
    Rinv = RotationM(-th)
    # Apply the rotation transformation to the squared dimensions.
    # Here we treat the bounding box dimensions similarly to the ellipse axes.
    I = np.sqrt(np.dot(Rinv, np.dot([[width_oct**2, 0], [0, height_oct**2]], R)))
    return I[0][0], I[1][1]


def drawOctagon(angle, side_length, xCenter, yCenter, record=False):
    """
    Draws a regular octagon based on the computed angle and side length.
    Parameters:
        angle (float): Orientation angle (in radians) for the octagon.
        side_length (float): The side length of the octagon (as computed by OctagonFitter).
        xCenter (float): x-coordinate of the octagon's center.
        yCenter (float): y-coordinate of the octagon's center.
        record (bool): If True, the function will record vertices instead of drawing.
    Note:
        For drawing, we convert the side length to the circumradius R.
        For a regular octagon:  R = side_length / (2 * sin(pi/8))
    """
    if side_length == 0:  # If there is no shape, stop!
        return
    
    # Convert the side length to the circumradius.
    R = side_length / (2 * np.sin(np.pi / 8))
    
    # Round the center position to whole numbers (pixels)
    xc = int(np.round(xCenter))
    yc = int(np.round(yCenter))
    
    # Define the number of sides for the octagon
    num_sides = 8
    
    # Lists to hold the coordinates of the octagon's vertices
    xCoords = []
    yCoords = []
    
    # Calculate the angle increment between each vertex (in radians)
    angle_increment = 2 * np.pi / num_sides
    
    # Loop to calculate the coordinates of each vertex
    for i in range(num_sides):
        # Calculate the angle for each vertex
        theta = angle + i * angle_increment
        
        # Calculate the x and y coordinates of the vertex using the circumradius R
        x = xc + R * np.cos(theta)
        y = yc + R * np.sin(theta)
        
        # Add the calculated vertex to the list
        xCoords.append(x)
        yCoords.append(y)
    
    if record:
        # If recording, save or return the coordinates as needed.
        return xCoords, yCoords
    else:
        # Draw the octagon by connecting the vertices.
        for i in range(num_sides):
            # Get the current vertex and the next (with wrap-around)
            x1, y1 = xCoords[i], yCoords[i]
            x2, y2 = xCoords[(i + 1) % num_sides], yCoords[(i + 1) % num_sides]
            # Draw a line between vertices
            plt.plot([x1, x2], [y1, y2], 'r-', linewidth=2)


def OctagonFitter(arr, usePrint=False):
    # Retrieves the number of rows (height) and columns (width) from the image array.
    height, width = arr.shape
    
    # Initialize sums for moment calculations
    xsum = 0.0
    ysum = 0.0
    x2sum = 0.0
    y2sum = 0.0
    xysum = 0.0
    bitCount = 0
    # Initializes variables to track the bounding box.
    min_x, max_x = width, 0
    min_y, max_y = height, 0
    
    # Loop through every pixel in the image
    for y in range(height):
        for x in range(width):
            if arr[y, x] != 0:
                # For every nonzero pixel, it adds the x and y coordinates (and their squares and product) to compute the first and second moments.
                xsum += x
                ysum += y
                x2sum += x * x
                y2sum += y * y
                xysum += x * y
                bitCount += 1
                
                # Track the bounding box
                if x < min_x: min_x = x
                if x > max_x: max_x = x
                if y < min_y: min_y = y
                if y > max_y: max_y = y
    
    # If no shape is found (no nonzero pixels), the function returns default zeros to avoid further computation errors.
    if bitCount == 0:
        return 0, 0, 0, 0, 0, 0
    
    # Computes the centroid (mean position) of the shape.
    xCenter = xsum / bitCount
    yCenter = ysum / bitCount
    
    # Computes the width and height of the bounding box.
    width_oct = max_x - min_x
    height_oct = max_y - min_y
    
    # Refined approximation of octagonal side length based on ideal octagonal proportions
    side_length = (width_oct + height_oct) / (2 + 2 * np.sqrt(2))
    
    # Compute angle using second-order moments
    u20 = x2sum / bitCount - xCenter ** 2
    u02 = y2sum / bitCount - yCenter ** 2
    u11 = xysum / bitCount - xCenter * yCenter
    # Uses the standard formula for orientation from second-order moments.
    theta = 0.5 * np.arctan2(2 * u11, u20 - u02)
    angle = np.degrees(theta) % 180
    
    if usePrint:
        print(f"Angle: {angle:.2f} degrees")
        print(f"Side Length: {side_length:.2f}")
        print(f"Bounding Box Width: {width_oct}, Height: {height_oct}")
        print(f"Center: ({xCenter:.2f}, {yCenter:.2f})")
    # Returns a tuple containing the orientation, estimated octagon side length, bounding box dimensions, and the centroid coordinates.
    return angle, side_length, width_oct, height_oct, xCenter, yCenter

def drawEllipse(angle, major, minor, xCenter, yCenter, maxY):
    """Draws an ellipse given its parameters."""
    if major == 0 or minor == 0:
        return [], []
    
    xc = int(np.round(xCenter))
    yc = int(np.round(yCenter))
    
    txmin = np.zeros(maxY)
    txmax = np.zeros(maxY)

    sint = np.sin(angle)
    cost = np.cos(angle)
    rmajor2 = 1.0 / (major / 2)**2
    rminor2 = 1.0 / (minor / 2)**2

    g11 = rmajor2 * (cost)**2 + rminor2 * (sint)**2
    g12 = (rmajor2 - rminor2) * sint * cost
    g22 = rmajor2 * sint**2 + rminor2 * cost**2
    k1 = -g12 / g11
    k2 = (g12**2 - g11 * g22) / g11**2
    k3 = 1.0 / g11
    
    ymax = int(np.floor(np.sqrt(np.abs(k3 / k2))))
    if ymax > maxY:
        ymax = maxY
    if ymax < 1:
        ymax = 1
    ymin = -ymax
    
    xCoordinates = []
    yCoordinates = []
    
    for y in range(ymax):
        j2 = np.sqrt(k2 * (y**2) + k3)
        j1 = k1 * y
        txmin[y] = int(np.round(j1 - j2))
        txmax[y] = int(np.round(j1 + j2))
    
    for y in range(ymin, ymax):
        x = txmax[-y] if y < 0 else -txmin[y]
        xCoordinates.append(xc + x)
        yCoordinates.append(yc + y)
    
    for y in range(ymax, ymin, -1):
        x = txmin[-y] if y < 0 else -txmax[y]
        xCoordinates.append(xc + x)
        yCoordinates.append(yc + y)
    
    return xCoordinates, yCoordinates

def EllipseFitter(arr, centroid2, usePrint=False):
    width=arr.shape[0]
    height=arr.shape[1]
    HALFPI = np.pi/2
    nCoordinates = 0 #int # Initialized by makeRoi()
    bitCount = 0 #int
    sqrtPi = np.sqrt(np.pi)
    # computeSums
    xsum = 0.0
    ysum = 0.0
    x2sum = 0.0
    y2sum = 0.0
    xysum = 0.0
    #int bitcountOfLine
    #double   xe, ye
    #int xSumOfLine
    for y in range(height):
        bitcountOfLine = 0
        xSumOfLine = 0
        #offset = y*width # int
        for x in range(width):
            if arr[x,y] != 0:
                bitcountOfLine+=1
                xSumOfLine += x
                x2sum += x * x
        
        xsum += xSumOfLine
        ysum += bitcountOfLine * y
        ye = y
        xe = xSumOfLine
        xysum += xe*ye
        y2sum += ye*ye*bitcountOfLine
        bitCount += bitcountOfLine
    # getMoments
    #double   x1, y1, x2, y2, xy
    if bitCount != 0:
        x2sum += 0.08333333 * bitCount
        y2sum += 0.08333333 * bitCount
        n = bitCount
        x1 = xsum/n
        y1 = ysum / n
        x2 = x2sum / n
        y2 = y2sum / n
        xy = xysum / n
        xm = x1
        ym = y1
        u20 = x2 - (x1 * x1)
        u02 = y2 - (y1 * y1)
        u11 = xy - x1 * y1
    m4 = 4.0 * np.abs(u02 * u20 - u11 * u11)
    if m4 < 0.000001:
        m4 = 0.000001
    a11 = u02 / m4
    a12 = u11 / m4
    a22 = u20 / m4
    xoffset = xm
    yoffset = ym
    tmp = a11 - a22
    if tmp == 0.0:
        tmp = 0.000001
    theta = 0.5 * np.arctan(2.0 * a12 / tmp)
    if theta < 0.0:
        theta += HALFPI
    if a12 > 0.0:
        theta += HALFPI
    elif a12 == 0.0:
        if a22 > a11:
            theta = 0.0
            tmp = a22
            a22 = a11
            a11 = tmp
        elif a11 != a22:
            theta = HALFPI
    tmp = np.sin(theta)
    if tmp == 0.0:
        tmp = 0.000001
    z = a12 * np.cos(theta) / tmp
    major = np.sqrt (1.0 / np.abs(a22 + z))
    minor = np.sqrt (1.0 / np.abs(a11 - z))
    scale = np.sqrt (bitCount / (np.pi * major * minor)) #equalize areas
    major = major*scale*2.0
    minor = minor*scale*2.0
    angle = 180.0 * theta / np.pi
    if angle == 180.0:
        angle = 0.0
    if major < minor:
        tmp = major
        major = minor
        minor = tmp
    yCenter, xCenter = [centroid2[0], centroid2[1]]
    # xCenter = left + xoffset + 0.5
    # yCenter = top + yoffset + 0.5
    if usePrint:
        print(angle)
        print(major,minor)
        print(xCenter,yCenter)
    return angle,major,minor,xCenter,yCenter
    
# Global variables to store click coordinates for y-scaling
y_x1, y_y1, y_x2, y_y2 = -1, -1, -1, -1
def two_clicks_y_scaling(event):
    """Handles two right-clicks to select points for y-axis scaling."""
    global y_x1, y_y1, y_x2, y_y2
    if event.button == 3:  # Right-click detected
        if y_x1 == -1:  # First click
            y_x1, y_y1 = event.xdata, event.ydata
            print(f"First click for y-scaling at: ({y_x1}, {y_y1})")
        else:  # Second click
            y_x2, y_y2 = event.xdata, event.ydata
            print(f"Second click for y-scaling at: ({y_x2}, {y_y2})")
            plt.close()

# Main function
path = r"C:\Users\user\Desktop\Masters\Fiber_Spacing\\"
ld = os.listdir(path)
ld.sort()

if not os.path.exists(path + 'rotated'):
    os.makedirs(path + 'rotated')

rough_scale = -1
y_rough_scale = -1  # Initialize y_rough_scale

for filename in ld:
    if filename[-4:] in '.bmp':
        fl = open(path + 'rotated/' + filename[0:-4] + '_.dat', 'w')
        flc = open(path+'rotated/'+filename[0:-4]+'_centroid.txt', 'w')
        img1 = Image.open(path + filename)

        # Create combined R+G+B image
        img_array = np.array(img1)
        combined_rgb = img_array[:,:,0] + img_array[:,:,1] + img_array[:,:,2]
            
        # Step 1: Calculate rough_scale (x-axis scaling)
        print("\033[91m# Click two points on the following image, one on the center of the first fiber of the first row and the second one on the center of the last fiber of the first row.\033[0m")
        fig, ax = plt.subplots()
        ax.imshow(np.array(img1), cmap='gray', interpolation='none', origin="lower")
        ax.set_title("Right-click twice to select points for x-axis scaling")
        fig.canvas.mpl_connect('button_release_event', two_clicks)
        plt.show()

        if (x1 == x2) and (y1 == y2):
            img1.save(path + 'rotated/' + filename)
            print(filename + '  :   grade  0 ')
        else:
            theta = math.degrees(math.atan((y2 - y1) / (x2 - x1)))
            img2 = img1.rotate(theta)
            img2.save(path + 'rotated/' + filename)
            npimg = np.array(img2.getdata()).reshape(img2.size[1], img2.size[0], 3)
            npimg = npimg[:, :, 0] + npimg[:, :, 1] + npimg[:, :, 2]
            fiber_diameter = float(input('Enter the diameter of a single fibre in micrometer:'))
            if rough_scale == -1:
                rough_scale = float(x2 - x1) / (27 * fiber_diameter)  # Calculate rough_scale for x-axis
                halfwidth = int(rough_scale * fiber_diameter/2)
                mask = np.zeros((2 * halfwidth + 1, 2 * halfwidth + 1))
                rsq = (fiber_diameter/2 * rough_scale)**2

                for ii in np.arange(2 * halfwidth + 1):
                    for jj in np.arange(2 * halfwidth + 1):
                        if ((ii - halfwidth)**2 + (jj - halfwidth)**2) < rsq:
                            mask[ii, jj] = 1

            # Step 2: Calculate y_rough_scale (y-axis scaling)
            print("\033[91m# Click two points on the following image, one on the center of the first fiber of the first column and the second one on the center of the last fiber of the same column.\033[0m")
            fig, ax = plt.subplots()
            ax.imshow(np.array(img2), cmap='gray', interpolation='none', origin="lower")
            ax.set_title("Right-click twice to select points for y-axis scaling")
            fig.canvas.mpl_connect('button_release_event', two_clicks_y_scaling)
            plt.show()

            if (y_x1 == y_x2) and (y_y1 == y_y2):
                print("Error: The same point was clicked twice for y-axis scaling.")
            else:
                y_rough_scale = float(y_y2 - y_y1) / (27 * (fiber_diameter*np.sqrt(3)/2))  # Calculate y_rough_scale for y-axis
                print(f"rough_scale (x-axis): {rough_scale}")
                print(f"y_rough_scale (y-axis): {y_rough_scale}")

            # Step 3: Select the center point of the fiber bundle
            print("\033[91m# Click one point which should be the centre point of the fibre bundle\033[0m")
            fig, ax = plt.subplots()
            ax.imshow(np.array(img2), cmap='gray', interpolation='none', origin="lower")
            ax.set_title("Right-click centre point of the fibre bundle")
            fig.canvas.mpl_connect('button_release_event', one_clicks)
            plt.show()

            print(f"Final coordinates: x1={x1}, y1={y1}")

            
            num_rows = int(input('Enter an odd number of fibre rows: '))
            max_fibers_in_middle = int(input("Enter number of fibers in middle row: "))
            row_spacing = y_rough_scale * fiber_diameter * (np.sqrt(3) / 2)  # Use y_rough_scale for row spacing
            fiber_centers = []
            fitted_centers = []
            holesizearr = np.zeros((max_fibers_in_middle, num_rows))
            
            for row in range(num_rows):
                row_index = row - (num_rows // 2)
                num_fibers = max_fibers_in_middle - abs(row_index)
                y = y1 + row_index * row_spacing
                x_start = x1 - (num_fibers - 1) * fiber_diameter * rough_scale / 2
                for i in range(num_fibers):
                    x = x_start + i * fiber_diameter * rough_scale
                    fiber_centers.append((x, y))
            # Calculate and plot centroids for each fiber position
            plt.figure(figsize=(12, 12))
            plt.imshow(npimg, cmap='gray', interpolation='none', origin="lower")
                        
            for idx, (x, y) in enumerate(fiber_centers):
                x_int = int(np.round(x))
                y_int = int(np.round(y))
                
                # Extract subimage around the expected fiber position
                subimg = npimg[y_int - halfwidth - 10:y_int + halfwidth + 10 , 
                               x_int - halfwidth-10:x_int + halfwidth + 10]
                centroid1 = ndimage.center_of_mass(subimg)
                centroid1_x = x_int - halfwidth - 10 + centroid1[1]
                centroid1_y = y_int - halfwidth - 10 + centroid1[0]
                
                # Create a grid of distances from the initial centroid
                y_indices, x_indices = np.indices(subimg.shape)
                distances = np.sqrt((x_indices - centroid1[1])**2 + (y_indices - centroid1[0])**2)
                
                # Determine the radius to threshold (you can adjust this factor)
                radius_threshold = (fiber_diameter * y_rough_scale) / 2 * 0.8  # 80% of expected radius
                
                # Create a mask for pixels within the radial distance
                mask = distances <= radius_threshold
                
                # Apply the mask to the subimg (set outside pixels to zero)
                thresholded_subimg = subimg.copy()
                thresholded_subimg[~mask] = 0
                
                # Calculate second centroid using only the thresholded region
                centroid2 = ndimage.center_of_mass(thresholded_subimg)
                centroid2_x = x_int - halfwidth - 10 + centroid2[1]
                centroid2_y = y_int - halfwidth - 10 + centroid2[0]
                # Plot the centroid
                plt.plot(centroid1_x, centroid1_y, 'r+', markersize=10)
                plt.plot(centroid2_x, centroid2_y, 'r+', markersize=10)
                plt.text(centroid2_x, centroid2_y, str(idx), color='yellow', fontsize=8)  
            plt.title('Detected Fiber Centroids')
            plt.show()
            plt.figure()
            plt.imshow(subimg)
            plt.show()
            # Optional: Save the centroid coordinates to a file
            with open(path + 'rotated/' + filename[0:-4] + '_centroids.txt', 'w') as f:
                for idx, (x, y) in enumerate(fiber_centers):
                    f.write(f"{idx}, {x}, {y}\n")
            ellipse_file = open(path + 'rotated/' + filename[0:-4] + '_ellipse_params.txt', "w")

            save_dir = os.path.join(path, "plots")
            os.makedirs(save_dir, exist_ok=True)
            # Inside your fiber centers loop
            fibre_core = int(input('Enter the fibre core diameter in microns: '))
            for idx, (x, y) in enumerate(fiber_centers):
                x_int = int(np.round(x))
                y_int = int(np.round(y))
                subimg = npimg[y_int - halfwidth -10:y_int + halfwidth + 10, x_int - halfwidth -10:x_int + halfwidth + 10]
                centroid1 = ndimage.center_of_mass(subimg)
                centroid1_x = x_int - halfwidth - 10 + centroid1[1]
                centroid1_y = y_int - halfwidth - 10 + centroid1[0]
                
                # Create a grid of distances from the initial centroid
                y_indices, x_indices = np.indices(subimg.shape)
                distances = np.sqrt((x_indices - centroid1[1])**2 + (y_indices - centroid1[0])**2)
                
                # Determine the radius to threshold (you can adjust this factor)
                radius_threshold = (fiber_diameter * y_rough_scale) / 2 * 0.8  # 80% of expected radius
                
                # Create a mask for pixels within the radial distance
                mask = distances <= radius_threshold
                
                # Apply the mask to the subimg (set outside pixels to zero)
                thresholded_subimg = subimg.copy()
                thresholded_subimg[~mask] = 0
                
                # Calculate second centroid using only the thresholded region
                centroid2 = ndimage.center_of_mass(thresholded_subimg)
                centroid2_x = x_int - halfwidth - 10 + centroid2[1]
                centroid2_y = y_int - halfwidth - 10 + centroid2[0]
                # Histogram analysis for each fiber
                bins = np.arange(np.max(thresholded_subimg)+1)
                hist, bin_edges = np.histogram(thresholded_subimg, bins)
                hist1 = hist[0]
                bin_centers = (bins[:-1] + bins[1:]) / 2
                
                # Smooth histogram and find peaks
                hist_smooth = gaussian_filter1d(hist, sigma=3)
                all_peaks, _ = find_peaks(hist_smooth, height=np.max(hist_smooth)*0.01, distance=100) #0.01
                D_p = 550  # 550 Default value
                
                # Plot and save each histogram
                plt.figure()
                if len(all_peaks) >= 2:
                    midpoint = np.median(bin_centers[all_peaks])
                    left_peaks = [p for p in all_peaks if bin_centers[p] < midpoint]
                    right_peaks = [p for p in all_peaks if bin_centers[p] >= midpoint]
                    
                    if left_peaks and right_peaks:
                        left_peak = max(left_peaks, key=lambda p: hist_smooth[p])
                        right_peak = max(right_peaks, key=lambda p: hist_smooth[p])
                        D_p = abs(bin_centers[left_peak] - bin_centers[right_peak])
                    print(f"Histogram {idx} - Peaks: {bin_centers[left_peak]:.2f}, {bin_centers[right_peak]:.2f}, D_p: {D_p:.2f}")
                    plt.scatter(bin_centers[all_peaks], hist_smooth[all_peaks], color='red', label='All Peaks')
                    plt.axvline(bin_centers[left_peak], color='green', linestyle="--", label="Left Parabola Peak")
                    plt.axvline(bin_centers[right_peak], color='orange', linestyle="--", label="Right Parabola Peak")
                plt.title(f'Histogram {idx}')
                plt.xlabel('Pixel Intensity')
                plt.ylabel('Number of Pixels')
                plt.legend()
                plt.plot(bin_centers, hist_smooth, label='Smoothed Histogram', color='blue')
                plt.savefig(os.path.join(save_dir, f"histogram_{idx}.png"))
                plt.close()
                
                # Thresholding and ellipse fitting
                thimg = (thresholded_subimg > (0.40 * D_p)) #* mask 0.40
                threshold_area = np.pi * (((y_rough_scale + rough_scale) / 2) * fibre_core)**2 #100
                minor = 0  # Initialize minor axis variable
                if np.sum(thimg) > threshold_area:
                    # Fit an ellipse to the thresholded sub-image.
                    # params = OctagonFitter(thimg)
                    # angle, side_length, width_oct, height_oct, xCenter, yCenter = params
                    angle, major, minor, xCenter, yCenter = EllipseFitter(thimg, centroid2)
                    xf = xCenter + (x_int - halfwidth)
                    yf = yCenter + (y_int - halfwidth)
                    fitted_centers.append((xf, yf))
                    ellipse_file.write(
                        f"{idx}: Major={major / rough_scale:.2f} px, "
                        f"Minor={minor / y_rough_scale:.2f} px, "
                        f"Angle={angle:.2f} deg, "
                        f"Center=({xf / rough_scale:.2f} um, {yf / y_rough_scale:.2f} um)\n"
                    )
                    print(f"Grid point {idx}: Ellipse fitted. Angle={angle:.2f} deg, "
                      f"Major={major:.2f} px, Minor={minor:.2f} px, "
                      f"Center=({xf:.2f}, {yf:.2f}) px")
                    # Store the minor axis value for grading
                    holesizearr[i, row] = minor / y_rough_scale #minor
                else:
                    fl.write(f"{i + 1} , {row + 1} , Hole smaller than 100 um.\n")
                    holesizearr[i, row] = 0  # Mark as invalid
            ellipse_file.close()
            # Grading should be after all fibers are processed
            valid_sizes = holesizearr[holesizearr > 0]
            if len(valid_sizes) > 0:
                min_size = np.min(valid_sizes)
                if min_size < 120.0:
                    grade = 0
                elif min_size < 130.0:
                    grade = 1
                elif min_size < 140.0:
                    grade = 2
                else:
                    grade = 3
            else:
                grade = 0  # No valid fibers found
            print(f"Final Grade for {filename}: {grade}")

            # After processing all fibers, display the result image with ellipses
        
            # Plot the hexagonal fiber bundle grid
            fiber_centers = np.array(fiber_centers)
            fitted_centers = np.array(fitted_centers)
        
            plt.figure(figsize=(8, 8))
            plt.scatter(fiber_centers[:, 0], fiber_centers[:, 1], c='cyan', s=100)
            for idx, (x, y) in enumerate(fiber_centers):
                plt.text(x, y, str(idx), fontsize=8, color='blue', ha='center', va='center')
            plt.gca().set_aspect('equal')
            plt.xlabel("X (pixel)")
            plt.ylabel("Y (pixel)")
            plt.title("Hexagonal Fiber Bundle Grid")
            plt.show()
        
            # Plot fitted centers overlaid on the original image
            plt.figure(figsize=(12, 12))
            if len(fitted_centers) > 0:
                    plt.scatter(fitted_centers[:, 0], fitted_centers[:, 1], c='red', s=70, label='Fitted Centers')
            for idx, (x, y) in enumerate(fiber_centers):
                    plt.text(x, y, str(idx), fontsize=8, color='blue', ha='center', va='center')
            plt.gca().set_aspect('equal')
            plt.xlabel("X (pixel)")
            plt.ylabel("Y (pixel)")
            plt.title("Fitted Hexagonal Fiber Bundle - Close the image to continue")
            plt.show()
        
            if not os.path.exists(path + 'rotated'):
                os.makedirs(path + 'rotated')
            for filename in ld:
                if filename[-4:] in '.bmp':
                        # Display the image with fiber centers overlaid
                        img3 = Image.open(path + 'rotated/' + filename)
                        img_width, img_height = img3.size
        
                        plt.figure(figsize=(12, 12))
                        plt.imshow(img3, cmap='gray', interpolation='none', origin="lower")
                        plt.scatter(fiber_centers[:, 0], fiber_centers[:, 1], c='red', s=30, label='Adjusted')
                        for idx, (x, y) in enumerate(fiber_centers):
                            plt.text(x, y, str(idx), fontsize=8, color='white', ha='center', va='center')
                        plt.gca().set_aspect('equal')
                        plt.xlabel("X (pixel)")
                        plt.ylabel("Y (pixel)")
                        plt.title("Hexagonal Fiber Bundle Grid - Close the image to continue")
                        plt.show()
        
                        # Overlay generated grid and fitted centers
                        plt.figure(figsize=(8, 8))
                        plt.imshow(npimg, cmap='gray', interpolation='none', origin="lower")
                        if len(fitted_centers) > 0:
                            plt.scatter(fitted_centers[:, 0], fitted_centers[:, 1], c='red', s=70, label='Fitted Centers')
                        for idx, (x, y) in enumerate(fiber_centers):
                            plt.text(x, y, str(idx), fontsize=8, color='blue', ha='center', va='center')
                        plt.scatter(fiber_centers[:, 0], fiber_centers[:, 1], c='cyan', s=50, label='Generated Grid')
                        plt.gca().set_aspect('equal')
                        plt.xlabel("X (pixel)")
                        plt.ylabel("Y (pixel)")
                        plt.title("Fiber Bundle Grid & Fitted Ellipse Centers")
                        plt.legend()
                        plt.show()
        
                        # Plot all detected ellipses on the same figure
                        plt.figure(figsize=(8, 8))
                        plt.imshow(npimg, cmap='gray', interpolation='none', origin="lower")
        
                        for idx, (xf, yf) in enumerate(fitted_centers):
                            # Fit the octagon to the thresholded image
                            # angle, side_length, width_oct, height_oct, xCenter, yCenter = OctagonFitter(thimg)
                            
                            # # Get the coordinates of the octagon vertices
                            # xCoords, yCoords = drawOctagon(np.radians(angle), side_length, xf, yf, record=True)
                            angle, major, minor, xCenter, yCenter = EllipseFitter(thimg, centroid2)
                            xCoords, yCoords = drawEllipse(np.radians(angle), major, minor, xf, yf, thimg.shape[0])
                            plt.plot(xCoords + [xCoords[0]], yCoords + [yCoords[0]], 'r-', linewidth=2)
                            plt.scatter([xf], [yf], color='blue', s=30)
        
                        plt.xlabel("X (pixel)")
                        plt.ylabel("Y (pixel)")
                        plt.title("Detected Ellipses on Fiber Cores")
                        plt.legend(["Ellipses", "Centers"])
                        plt.show()
