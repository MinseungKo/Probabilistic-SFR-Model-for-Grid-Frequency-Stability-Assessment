# Probabilistic System Frequency Response (SFR) Model for Grid Frequency Stability Assessment
MATLAB Code for the paper "Voltage-Coupled Probabilistic Frequency Security Assessment of Inverter-Dominated Grid."

## Code Overview
This code aims to enhance the system frequency response (SFR) model, which is widely used for assessing power system frequency.
Three improvements are made as follows:

### 1) Unified SFR model with IBRs
In addition to synchronous generators and fast-frequency response resources, grid-following and grid-forming inverters are added with realistic control settings.

### 2) Load bus voltage generation method
Based on the generator bus voltage, we estimate the load bus voltage to incorporate the voltage-dependent behaviors of the load and DER.
The conventional estimation method uses the fixed voltage sensitivity of the generator to the load bus, which represents the load bus voltage as a linear combination of the generator bus voltages.
Our method enhances this approach by utilizing the SVD results on V-Q sensitivity. This enables us to obtain a more accurate load bus voltage estimation.

We further consider a case where we do not have measurements for the target case (instead, we assume that we have measurements for the other reference case.)
In this case, we first estimate the generator bus voltage for the target case based on SVD & orthogonal alignment.
Then, we conduct the same process with the previous case and can also get accurate load voltage estimation results for the target case.

### 3) Stochastic SFR model
We modify the conventional deterministic SFR to incorporate stochastic elements, thereby accounting for various types of power system uncertainties. Solutions can be obtained by solving stochastic differential equations.
Specifically, uncertainties implemented in this code are as follows:
  i) Steady-state uncertainties:
  Load amount, IBR output, BTM DER distribution & output
  ii) Dynamic uncertainties:
  Load model parameters, Voltage recovery, Voltage jump, Completely unknown uncertainties

## File Description
- Fin_samescene_SFR.m: main file to run probabilistic SFR when we have voltage measurements for the target case
- Fin_diffscene_SFR.m: main file to run probabilistic SFR when we have voltage measurements for the reference case
- loadbus_K_B.xlsx: PSS/E simulation output on IEEE 39-bus system. A (K) MW load increase at Bus (B) is simulated.
- volt_sensitivity.csv: Voltage-to-Voltage sensitivity between load/generator bus to generator bus. This is obtained by injecting reactive power perturbation signals at generator buses. Simulations are conducted by using PSS/E.
- voltH.csv: Voltage-to-Reactive power sensitivity. We use this information for our proposed estimation of the load bus voltages.
