# Baseline correction via Adaptive Smoothness Penalized Least Squares
The **asPLS** (*Adaptive Smoothness Penalized Least Squares*) algorithm is a method to remove spurious baseline trends that might be present in timeseries data. It is especially useful in the context of peak detection, where spurious trends can completely break the identification and analysis.
The retionale behind this method is to fit and remove the baseline, but not the peaks present in the timeseries, thereby preserving the relevant information.<br>
This repository proposes a Matlab implementation of asPLS based in the [original paper](https://www.tandfonline.com/doi/abs/10.1080/00387010.2020.1730908) by Zhang et al. [1]

![asPLS_removal (1)](https://github.com/CNelias/asPLS/assets/34754896/951d28eb-40d1-4989-9a1c-ec9168dda1cd)

## Usage
Download **aspls.m**, and place it in the folder where you intend to run it. To perform baseline correction, call the ```aspls``` function. 

If X is you data, running ```[detrended_signal, baseline] = aspls(X);``` returns the input series without baseline (detrended_signal) and the identified baseline trend(baseline). This however might not always produce optimal results, since it uses the default value for **the most important parameter for the analysis**: ```lambda```.
```lambda``` controls the smoothness of the fit, the higher it is, the straighter the estimated baseline is going to be. If ```lambda``` is choosen to low, the estimated baseline might start to fit the peaks in the timeseries, which is usually undesirable. The fitted baseline is not extremely sensitive to small changes in the value of ```lambda```, it is therefore common practice to vary ```lambda``` as power of 10, i.e. ```lambda = 10e5, 10e6 or 10e7...```. Other parameters of the aspls have less influence on the final estimated baseline, so we refer the reader to the original paper [1] for more details. <br>

Here are all the parameters of ```aspls```:
```
Inputs:
   - X : array, The y-values of the measured data, with t data points. 
       Must not contain missing data (NaN) or Inf.
   - lambda : float, The smoothing parameter. Larger values will create smoother baselines.
       Default is 1e5.
   - order : int, The order of the differential matrix. Must be greater than 0. Default is 2
       (second order differential matrix). Typical values are 2 or 1.
   - k : int, coefficient of hardness for the decay of the weights as
       a function of distance to baseline. A higher k means a harder
       decay which will push the baseline towards lower value of the
       input data. 0.5 is empirically recommended (Default).
   - itermax : int, The max number of fit iterations. Default is 20.
   - epsilon : float, The exit criteria. Default is 1e-4.
 
Returns
   - detrended_signal : array, The input signal - the calculated baseline.
   - baseline : array, The calculated baseline.
```






[1] Zhang, F., Tang, X., Tong, A., Wang, B., Wang, J., Lv, Y., ... & Wang, J. (2020). Baseline correction for infrared spectra using adaptive smoothness parameter penalized least squares method. Spectroscopy Letters, 53(3), 222-233.
