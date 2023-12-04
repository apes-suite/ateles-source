-- Sample configuration for solve_euler_riemann

-- Isentropic expansion coefficient
isen_coef = 1.4

-- State left in the order of density, velocity, pressure
velocityL = 0.0
densityL = 1.0
pressureL = 1e5

left = { 1.0, 0.0, 1e5 }

speedOfSoundL = math.sqrt(1e5*1.4)
speedOfSoundR = speedOfSoundL - (300*0.4/2)
pressureRatioRL = (speedOfSoundR / speedOfSoundL)^(2*1.4/(1.4 -1))
print(pressureRatioRL)

densityR = densityL * pressureRatioRL^(1/1.4) 

velocityR = velocityL - ((2 * speedOfSoundL/(1.4-1))* (pressureRatioRL^(0.4/(2*1.4))-1))

pressureR = pressureL*(speedOfSoundR / speedOfSoundL)^(2*1.4/(1.4 -1)) 

-- State right in the order of density, velocity, pressure
right = { densityR, velocityR, pressureR }

-- Write resulting time slices into files
create_files = true

-- Shift the origin by the given value.
-- (the original discontinuity is located at the given x-value)
origin_offset = 0.4

-- file_prefix = 'exact'

-- Auxilary function to read data from an ascii file
function readcsv_col(file, col)
  local fp = assert(io.open (file))
  local csv = {}
  for line in fp:lines() do
    local row = {}
    for value in line:gmatch("[^ ]*") do
      if value ~= '' then
        row[#row+1] = value
      end
    end
    csv[#csv+1] = row[col]
  end

  -- csv contains now given column from the csv file
  return csv
end

-- Points in time, where to get sample the solution:
--t = {0.0, 0.1, 0.2, 0.3, 0.4}
t = {8e-04}
--for i=0,4 do
--  table.insert(t, i*0.1)
--end

-- Points in space, where to probe the solution in every time slice:
-- x = {-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5}
-- x = readcsv_col('xdat.ascii', 1)
x = {}
for i=0,1000 do
  table.insert(x, i*0.00064)
end
