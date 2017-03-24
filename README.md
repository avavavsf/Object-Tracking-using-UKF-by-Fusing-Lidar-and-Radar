# Unscented Kalman Filter Project Starter Code
[//]: # (Image References)
[image1]: ./output_images/track1.png
[image2]: ./output_images/track1_radarNIS.png
[image3]: ./output_images/track1_laserNIS.png
[image4]: ./output_images/track2.png
[image5]: ./output_images/track2_radarNIS.png
[image6]: ./output_images/track2_laserNIS.png
[image7]: ./output_images/ctrv.png

## Workflow
The workflow is the same as the [EKF fusion project](https://github.com/ymshao/Object-Tracking-using-EKF-by-Fusing-Lidar-and-Radar/blob/master/README.md), but there are two differences:
  * UKF used instead of EKF to fusion Radar and Laser data.
  * ctrv model used as the process model with 5 parametes in the state vector.
![alt text][image7]
## Input data format
Input data is located in the "data" directory. Here is the data format description.
```
#Input file format:
#L(for laser) meas_px meas_py timestamp gt_px gt_py gt_vx gt_vy
#R(for radar) meas_rho meas_phi meas_rho_dot timestamp gt_px gt_py gt_vx gt_vy
-----------------------------
Example
-----------------------------
R	8.60363	0.0290616	-2.99903	1477010443399637	8.6	0.25	-3.00029	0
L	8.45	0.25	1477010443349642	8.45	0.25	-3.00027	0 

```

## Results
We have generated two simulated tracking data to test our implementation, here is the tracking results. We noticed that the UKF fusion resutls is more smooth than [EKF](https://github.com/ymshao/Object-Tracking-using-EKF-by-Fusing-Lidar-and-Radar/blob/master/README.md) at the circle part, and we archive a smaller RMSE values. 

We also found that it is difficult to archive the best RMSE and the best filter consistancy (indicated by the NIS figures), by only adjusting the process noise parametes.In this project, we scrifice the NIS to archive a little bit to archive best RMSE values.
  * Track 1.       px, py, vx, and vy RMSE: [0.0740189, 0.0837208, 0.583002, 0.586191] meters
![alt text][image1]
  * Track 1 Radar NIS
![alt text][image2]
   * Track 1 Laser NIS
![alt text][image3]
  
  * Track 2.       px, py, vx, and vy RMSE: [0.197399, 0.189804, 0.548703, 0.543699] meters
![alt text][image4]
  * Track 2 Radar NIS
![alt text][image5]
  * Track 2 Laser NIS
![alt text][image6]

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`


