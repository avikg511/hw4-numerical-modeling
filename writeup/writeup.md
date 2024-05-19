# SIOC 224 - Homework 4
## Avik Ghosh

### Code:
Here's the GitHub link for it!
\url{https://github.com/avikg511/hw4-numerical-modeling.git}

## Question 1 - Diffusion Equation

### 1a - von Neumann stability analysis for Matsuno Time Stepping (Centered Difference space)
Let's consider the equation $$\frac{\partial C}{\partial t} = \kappa \frac{\partial^2 C}{\partial x^2}$$

First, let's discretize in space. This would give us the semidiscrete equation, so $$\left( \frac{\partial C}{\partial t} \right)_k = \kappa \left( \frac{C_{k+1} - 2C_k + C_{k-1}}{(\Delta x)^2}\right)$$

where $C_k$ refers to the kth space step. Writing Implicit Euler Scheme, we need to have $C_k^{n+1} = C_k^n + \Delta t * \left(\frac{\partial C}{\partial t}\right)_k^{n+1}$, where $C^{n+1}$ refers to the $(n+1)^{th}$ timestep, and similarly for $C^n$. But note that we're doing the predictor corrector step, so we use Forward Euler to determine the $(n+1)$ derivative. Here, we get a set of test values using $$(\vec{C^{n+1}})^* =  \vec{C}^n + \Delta t * \vec f(\vec C^n)$$ 
where $\vec f$ refers to the derivative acting on the vector $\vec C$. Then we can use $(\vec{C^{n+1}})^*$ to calculate the derivative at time $(n+1)$. Without the vector notation, we have $$(C_k^{n+1})^* = C_k^n + \Delta t \cdot (C_{k+1}^n - 2C_k^n + C_{k-1}^n)$$

Therefore, let's write $$C^{n+1} = C_k^n + \Delta t * \vec f((\vec{C^{n+1}})^*)$$ or writing it with the index k, we can see that $$C_k^{n+1} = C_k^n + \Delta t * \frac{\kappa}{(\Delta x)^2} \cdot \left( (C_{k+1}^{n+1})^*  - 2(C_{k}^{n+1})^* + (C_{k}^{n+1} )^*\right)$$ 

Plugging in for $\sigma = \frac{\kappa \Delta t}{(\Delta x)^2}$, we can write it as $$C_k^{n+1} = C_k^n + \sigma \cdot \left( (C_{k+1}^{n+1})^*  - 2(C_{k}^{n+1})^* + (C_{k}^{n+1} )^*\right)$$ 

Let's also plug in for the test step, so we can write the entire equation as $$\begin{aligned} C_k^{n+1} = C_k^n + \sigma \cdot ([C_{k+1}^n + \sigma \cdot (C_{k+2}^n - 2C_{k+1}^n + C_{k}^n)] \\ + [C_{k}^n + \sigma \cdot (C_{k+1}^n - 2C_{k}^n + C_{k-1}^n)]  + [C_{k-1}^n +  \sigma \cdot (C_{k}^n - 2C_{k-1}^n + C_{k-2}^n)]) \end{aligned}$$

Rewriting it, we can say $$\begin{aligned} C_k^{n+1} = C_k^n + \sigma (C^n_{k+1} -2 C_k^{n} + C^n_{k-1}) + \\ \sigma^2 ((C_{k+2}^n - 2C_{k+1}^n + C_{k}^n) -2(C_{k+1}^n - 2C_{k}^n + C_{k-1}^n) + (C_{k}^n - 2C_{k-1}^n + C_{k-2}^n)) \end{aligned}$$ 

Then we expand and write

$$C_{k}^{n+1}=C_{k}^n + \sigma(C^n_{k+1} -2 C_k^{n} + C^n_{k-1}) + \sigma^2 (C_{k+2}^n - 4C_{k+1}^n + 6C_{k}^n -4C_{k-1}^n + C_{k-2}^n)$$

We plug in for $C_n^k = A \rho^n e^{jmk\Delta x}$. We're going to constrain A to be nonzero so we divide out the A terms on both sides. Therefore, we have:


$$\begin{aligned} \rho^{n+1}e^{jmk\Delta x} = \\
\rho^{n}e^{jmk\Delta x} + \sigma (\rho^{n}e^{jm(k+1)\Delta x}-2\rho^{n}e^{jmk\Delta x} + \rho^{n}e^{jm(k-1)\Delta x}) + \\
\sigma^2 (\rho^{n}e^{jm(k+2)\Delta x} - 4\rho^{n}e^{jm(k+1)\Delta x} + 6\rho^{n}e^{jmk\Delta x} - 4\rho^{n}e^{jm(k-1)\Delta x} + \rho^{n}e^{jm(k-2)\Delta x}) \end{aligned}$$

Dividing away the $\rho^n e^{jmk\Delta x}$ from both sides, we're left with:

$$\rho = 1 + \sigma (e^{jm\Delta x} -2 + e^{-jm\Delta x}) + \sigma^2(e^{2jm\Delta x} -4e^{jm\Delta x} + 6(1) - 4e^{-jm\Delta x} + e^{-2jm\Delta x})$$

The first term is 1, the second was given in lecture as $4\sigma sin^2(\frac{m\Delta x}{2})$. The last term has the coefficients from Pascal's Triangle that are from a binomial raised to the fourth power, and we know the middle term is $6x^2y^2 = 6$, so x and y must be complex conjugates of one another with magnitude 1. Also, the fourth order term has an exponent of 2, so it looks like 

$$\begin{aligned} \left(\exp\left(\frac{jm\Delta x}{2}\right) - \exp\left(\frac{-jm\Delta x}{2}\right)\right)^4 = \left(2j \cdot sin\left(\frac{m \Delta x}{2}\right) \right)^4 = \\[8pt]
\left( 2^4 j^4 sin^4\left(\frac{m \Delta x}{2}\right)\right) = 16sin^4\left(\frac{m \Delta x}{2}\right) \end{aligned}$$

Plugging it all in, we see

$$\rho = 1 - 4\sigma sin^2\left(\frac{m \Delta x}{2}\right) + 16 \sigma^2 sin^4\left(\frac{m \Delta x}{2}\right)$$

We're again considering the most dangerous mode, and considering that $m = \frac{\pi}{\Delta x}$ causes a discrete time cosine to flip with each timestep, we pick that and note that the sines go to 1. Therefore we have $\rho = 1 - 4\sigma + 16 \sigma^2$. We note that since we're modeling $C^n_k = A\rho^ne^{jmk\Delta x}$, we need $\rho <1$ to keep our system bounded. Therefore we solve for $\rho = 1-4\sigma +16\sigma^2 < 1$. (I don't think a plot is entirely necessary here so I just solved algebraicly). 

$$ \begin{aligned} 1 - 4\sigma + 16\sigma^2 < 1  \\[6pt]
16 \sigma^2 - 4\sigma < 0 \\[6pt]
16 \sigma^2 < 4 \sigma \\[6pt]
\sigma^2 < \frac{1}{4} \sigma \\[6pt]
\sigma < \frac{1}{4} \end{aligned}$$

where in the last step, we used the fact that $\sigma = \frac{\kappa \Delta t}{\Delta x ^2} > 0$. Therefore, we now know that $$ \begin{aligned} \sigma < \frac{1}{4} \implies \\ \rho < 1 \implies \end{aligned}$$ our system is bounded and therefore stable.

### 1b - Nondimensionalizing + FE

We have $$\frac{\partial C}{\partial \tilde t} = \kappa \frac{\partial^2 C}{\partial \tilde x ^2}$$ 

where $\tilde x, \tilde t$ represent the dimensional quantities of position and time. Let $\tilde t = \frac{t}{u_0}$, where $[u_o] = \rm seconds$. (Not sure how standard this is but [] refers to units here). Similarly, define $\tilde x = \frac{\kappa}{u_0}$. Then, plugging in for the nondimensionalized quantities, we get $$\frac{\partial C}{\partial t} u_0 =  \kappa \frac{\partial^2 C}{\partial x^2} \frac{u_0}{\kappa}$$

Simplifying by dividing out the $\kappa, u_0$, we get $$\frac{\partial C}{\partial t} = \frac{\partial^2 C}{\partial x^2}$$

(nondimensionalized).

Here, $\sigma$ is actually a different value. We still divide by the $\Delta x^2$ for the semidiscrete equation, and multiply by $\Delta t$, but now we don't really need the $\kappa$, so $\sigma = \frac{\Delta t}{\Delta x^2}$.

### Plot of RMS Error vs Time Forward Euler
![](images/image.png)
Plot of RMS error vs Time (x axis). The dt range is 0 to 0.0014 for time, which corresponds to 0 to about 0.56 for sigma. Duration at 0.3

The minimum is at about 0.0004-0.0005, which corresponds to $\sigma = \frac{0.0004}{0.05^2} = 0.16$ to about $\sigma \approx 0.2$, which is close to the $\frac{1}{6}$ we were expecting. This plot starts to show the divergence at values above 0.5, but as mentioned in the question, the duration is too small. 

I stepped through the times with a resolution of $ $5 \times 10^{-7}$ to get the grid this fine. For the next graph, I increased duration to about 0.5 seconds.

![](images/image-1.png)
Plot of RMS error vs Time (x). There is a strong peak at about 0.0013 units time, and that corresponds to about sigma 0.56 (using $\frac{dt}{dx^2}$)

Here it's a lot clearer that the 0.5 bound is accurate. One thing that was a little concerning was how much noise was in this data, and I'm curious as to why, and if that's really just floating point math again.

![](images/image-2.png)
Main method that controlled FE code. Sigma is a temporary set of values that correspond to a proper range for sigma. Scaling that gives the tVals array for the range of time

![](images/image-3.png)
Class outline (header file). Here I'm actually using a real matrix class (boost::numeric::ublas) to make my life simpler, instead of doing my own indexing system. In the end, they also store it in a single line of memory so there isn't any performance tradeoffs, which is nice. data contains all the concentration values at all time, indexed by time (row) and space (collumn)

![](images/image-4.png)
Constructor that includes initial conditions

![](images/image-5.png)
Forward Euler Code

![](images/image-6.png)
RMS error

### 1c - Matsuno RMS Plot
![](images/image-7.png)
RMS Error plot vs sigma (x), semilog. Error begins to diverge right around 0.25, demonstrating the derivation was correct.

![](images/image-8.png)
RMS Error plot again, but this time against time. I plotted the $\sigma$ first just because it's easier to read it than convert in my head from timestep sizes to sigma. Here, the divergence happens around t = 0.00063 or so. Given that sigma is still dt / dx^2, and dx = 0.05, which is about 0.252 for sigma

![](images/image-12.png)
Matsuno Code to solve diffusion equation. Index calculation method was separated because it's the same code multiple times.

![](images/image-11.png)
Matsuno Index Calculation code used in Models::runMatsuno above

### 1d - Implicit Euler Plot
![](images/image-9.png)
RMS Error plot vs time. Error is always stable, and increasing time resolution will decrease error.

![](images/image-13.png)
![](images/image-14.png)
Code to solve the Implicit Euler method.