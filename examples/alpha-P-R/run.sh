# use the stress strain approach to calculate the elastic tensor
m2d init -n 9 -a stress -m 0.015 -b
m2d run -a stress input.yaml

# use the energy approach to calculate the elastic tensor
m2d init -n 9 -a energy -m 0.02 -b
m2d run -a energy input.yaml

# calculate the stress strain curve along different direction
m2d init -a stress -r 0.0  0.1 -n 11 -p ssc -d 'xx' 'yy' 'xy' 'bi' -b
m2d run -a stress -p ssc  input.yaml
