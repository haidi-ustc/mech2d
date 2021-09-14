# use the stress strain approach to calculate the elastic tensor
m2d init -n 9 -a stress -m 0.015 -b
m2d run -a stress input.yaml

# use the energy approach to calculate the elastic tensor
m2d init -n 9 -a energy -m 0.025 -b 
m2d run -a energy input.yaml

