# TCP/NTCP Radiobiological model
Python code used to model the variation in TCP and NTCP with dose fraction specific dose variations.
There is huge flexibility in variation of the delivered dose per fraction.
Doses can be generated within the model (as used heavily within the thesis) or can be supplied as a list of doses (which must be done for the NTCP model; a list of doses can easily be gneerated from the TCP model and then supplied for NTCP calcualtions).

## Example usage
This example generates a set of TCP reuslts for 1000 patients and then plots the individual results and the population mean.
There are 1000 patients with alpha/beta mean of 10 and standard deviation of 20% of the mean.
Nominal dose is 2Gy/#. There is no offset of dose from the mean, but delivered dose per fraction has a standard deviation of 1.5%.
A 5% drift in beam output per year is included and specific results are returned at 74Gy.
Usage of the NTCP model is similar.

```python
# Import the modules
import TCP_NTCP as model
import matplotlib.pyplot as plt


# Generate a set TCP of results
results = model.completeTCPcalc(n=1000,                     # number of patients in population to model
                                alphabeta_use=10,           # mean alpha/beta
                                alphabeta_sd_use=20,        # SD of alpha/beta (%)
                                d=2,                        # nominal dose (Gy/fraction)
                                d_shift=0,                  # initial dose difference (%)
                                d_sd=1.5,                   # standard deviation of delivered dose (%)
                                d_trend=5/365,              # dose drift (%/day)
                                max_d=100,                  # maximum dose for which TCP is calcualted (Gy)
                                dose_of_interest=74,        # results of TCP at this dose are returned seperately for simpler analysis.
                                n0 = 74)                    # Specified N0 value (can be determined by fitting to a defined population)


# Produce a plot the results

# plot the first 100 individual curves
for i in range(100):
    plt.plot(results['nom_doses'],results['TCPs'][i],c='C1',alpha=0.3,lw=1)
    
# plot population mean
plt.plot(results['nom_doses'],results['TCP_pop'],ls='-',c='C0',marker='o',ms=3)
plt.show()

# Generate set of NTCP results

ntcp_results = model.complete_NTCP_calc(d_data=[64,74],             # dose values on NTCP curve
                                       ntcp_data=[0.36,0.47],       # NTCP values corresponding to the above doses
                                       irrad_perc = 100,            # scaling factor
                                       frac_doses=all_d_array,      # this will be a list of dose values to compute the NTCP for
                                       max_dose=100,                # maximum calculation dose
                                       ntcp_params={'td50_1':(31.4,0),
                                                    'v': (0.72,0), ## v=1 for parotid (e.g. low for rectum (0.08,10))
                                                    'm':(0.53,0),   # m paremter in NTCP equation (may have ot be found through optimisation
                                                    'n':(1,0)},     # n parameter in NTCP equation
                                       fit_vals=False)['patient_ntcps'] #[0][29] ## this will reuturn the NTCP after 30#
```
