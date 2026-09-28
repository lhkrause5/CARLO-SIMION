--[[
CARLO FINAL production run — 3-pitch / 0.0025 mm/gu / 100k / 2-group MEMORY-SAFE
========================================
SIMION 8.2 workbench user program.

PURPOSE
-------
FINAL HIGH-RESOLUTION PRODUCTION SWEEP.

This file is for the final converged geometry:
  carlo_six_grid_final_3p_0025_2adj.pa0
with 3-pitch half-span and 0.0025 mm/gu resolution.

It runs the full nine-case production voltage sweep with the accepted
1000 K / exact-4-degree O+ phase-space ensemble.  This is intended to
serve as the final high-resolution, laterally converged paper baseline.

IMPORTANT PHYSICAL SCOPE
------------------------
The adopted electrical design maintains the continuous collector near
spacecraft common through the virtual-ground input of the current amplifier.
Therefore the collector is a FIXED 0-V conductor in the final geometry.

The final PA uses only two adjustable electrostatic groups:
  adjustable 1 = G3 + G4 retarding pair = +V_R
  adjustable 2 = G6 suppressor          = -15 V

Fixed at 0 V:
  G1, G2, G5, collector, computational upstream reference plane.

The program assumes:
  carlo_six_grid_final_3p_0025_2adj.pa0

with workbench coordinates:
  x,y = 0 .. 3.12 mm
  z   = 0 .. 5.85 mm
  central lattice origin = (1.56, 1.56) mm
  G1 entrance            = 1.00 mm
  G6 exit                = 3.80 mm
  collector front        = 5.80 mm

The PA instance must be unrotated/unscaled relative to these coordinates.

PARTICLE DEFINITION
-------------------
Load the supplied deterministic placeholder particle file:
  carlo_placeholder_100000.fly2

The FLY2 file only creates the required number of particles. This Lua
program replaces every particle's mass, charge, position, and velocity in
segment.initialize().

The Monte Carlo sampling is reproducible because SIMION's RNG is reseeded
to the same value at the start of every voltage run. Therefore each bias
state sees the same phase-space ensemble.

CURRENT NORMALIZATION
---------------------
Launch positions are uniform over one primitive hexagonal-lattice cell.
The raw collector estimator therefore contains the cell open fraction.

Raw explicit-cell current:
    I_raw = q n A_face < v_n * C >

where C=1 for a collector hit and 0 otherwise.

To preserve the accepted single-head amplitude convention, also report:
    I_adopted = I_raw * tau_stack / phi_hex

with tau_stack = 0.40 and phi_hex = 0.410855070...

This is the normalization proposed in Phase 3A. BOTH raw and adopted
currents are written so the physics remains traceable.

Do not multiply by phi_hex again or use phi_hex^6.

PRODUCTION STATISTICS
---------------------
This version is intended for 100,000 launched particles per voltage case.
With nine cases, one Fly'm executes 900,000 trajectories total.
Trajectory retention remains disabled for production statistics.

PRODUCTION VOLTAGE CASES
------------------------
One Fly'm executes the full production set:
  control: all electrodes equipotential at 0 V
  VR = 0.00 V   (G6 remains -15 V)
  VR = 3.00 V
  VR = 4.00 V
  VR = 4.50 V
  VR = 4.75 V
  VR = 5.00 V
  VR = 6.00 V
  VR = 7.00 V

For biased cases:
  G1=0, G2=0, G3=+VR, G4=+VR, G5=0, G6=-15 V,
  collector=0 V (virtual ground), upstream reference plane=0 V.

OUTPUT
------
  carlo_FINAL_3p_0025_1000K_4deg_100k.csv
--]]

simion.workbench_program()

-- User-adjustable model controls.
adjustable minimum_recommended_particles = 100000
adjustable random_seed = 9142026
adjustable lattice_center_x_mm = 1.56
adjustable lattice_center_y_mm = 1.56
adjustable launch_z_mm = 0.50
adjustable collector_witness_z_mm = 5.795
adjustable production_trajectory_image_control = 1

-- Accepted physical/model constants.
local q_C            = 1.602176634e-19
local density_m3     = 1.0e11
local A_face_m2      = 6.25e-4
local tau_stack      = 0.40
local phi_hex        = 0.410855070047639

local pitch_mm       = 0.52
local hex_row_mm     = math.sqrt(3) * pitch_mm / 2

local ion_mass_amu   = 16.0
local ion_charge_e   = 1.0

local mean_vx_mm_us  = 0.530149200455
local mean_vy_mm_us  = 0.0
local mean_vz_mm_us  = 7.581486781975
local sigma_v_mm_us  = 0.720870247558

local accepted_modelB_fieldfree_nA = 4.385028158
local suppressor_V = -15.0
local upstream_reference_V = 0.0

local cases = {
  {label="equipotential_control", equipotential=true,  VR=0.00},
  {label="VR_0p00",              equipotential=false, VR=0.00},
  {label="VR_3p00",              equipotential=false, VR=3.00},
  {label="VR_4p00",              equipotential=false, VR=4.00},
  {label="VR_4p50",              equipotential=false, VR=4.50},
  {label="VR_4p75",              equipotential=false, VR=4.75},
  {label="VR_5p00",              equipotential=false, VR=5.00},
  {label="VR_6p00",              equipotential=false, VR=6.00},
  {label="VR_7p00",              equipotential=false, VR=7.00},
}
local current_case = cases[1]

local diag_names = {"after_G1","after_G2","after_G3","after_G4","after_G5","after_G6"}
local diag_z = {1.30,1.85,2.40,2.95,3.50,4.80}

local outfile = nil
local n_initialized = 0
local n_positive_flux = 0
local weight_mps = {}
local last_z = {}
local counted_collector = {}
local counted_diag = {}
local sum_w_all = 0.0
local sum_w_col = 0.0
local sum_w2_col = 0.0
local collector_hits = 0
local diag_count = {}
local diag_sum_w = {}
local first_col_i = {}
local first_col_j = {}
local changed_column = {}
local samecol_hits = 0
local neighbor_hits = 0
local sum_w_samecol = 0.0
local sum_w_neighbor = 0.0

local function urand()
  local u = simion.rand()
  if u <= 1e-15 then u = 1e-15 end
  return u
end

local function gaussian_pair()
  local u1 = urand()
  local u2 = urand()
  local r = math.sqrt(-2.0 * math.log(u1))
  local a = 2.0 * math.pi * u2
  return r * math.cos(a), r * math.sin(a)
end

local function sample_launch_position()
  local u = simion.rand() - 0.5
  local v = simion.rand() - 0.5
  local x = lattice_center_x_mm + u*pitch_mm + v*(0.5*pitch_mm)
  local y = lattice_center_y_mm + v*hex_row_mm
  return x, y
end

local function nearest_hex_indices(x_mm, y_mm)
  local xr = x_mm - lattice_center_x_mm
  local yr = y_mm - lattice_center_y_mm
  local j0 = math.floor(yr / hex_row_mm + 0.5)
  local best_d2 = math.huge
  local best_i, best_j = 0, 0

  for j = j0-2, j0+2 do
    local offset = 0.0
    if (j % 2) ~= 0 then offset = 0.5*pitch_mm end
    local i0 = math.floor((xr - offset)/pitch_mm + 0.5)
    for i = i0-2, i0+2 do
      local xc = i*pitch_mm + offset
      local yc = j*hex_row_mm
      local dx = xr - xc
      local dy = yr - yc
      local d2 = dx*dx + dy*dy
      if d2 < best_d2 then
        best_d2 = d2
        best_i, best_j = i, j
      end
    end
  end
  return best_i, best_j
end

local function crossed_forward(z0, z1, zplane)
  return z0 < zplane and z1 >= zplane and ion_vz_mm > 0
end

local function reset_run_state()
  n_initialized = 0
  n_positive_flux = 0
  weight_mps = {}
  last_z = {}
  counted_collector = {}
  counted_diag = {}
  sum_w_all = 0.0
  sum_w_col = 0.0
  sum_w2_col = 0.0
  collector_hits = 0
  diag_count = {}
  diag_sum_w = {}
  first_col_i = {}
  first_col_j = {}
  changed_column = {}
  samecol_hits = 0
  neighbor_hits = 0
  sum_w_samecol = 0.0
  sum_w_neighbor = 0.0

  for k=1,#diag_names do
    diag_count[k] = 0
    diag_sum_w[k] = 0.0
    counted_diag[k] = {}
  end
end

local function case_electrode_voltages()
  if current_case.equipotential then
    return 0,0,0,0,0,0,0,0
  else
    return 0,0,current_case.VR,current_case.VR,0,suppressor_V,
           0.0,upstream_reference_V
  end
end

-- MEMORY-SAFE FAST ADJUST FOR THE FINAL TWO-GROUP PA
--
-- adjustable electrode 1 = G3 + G4 retarding pair
-- adjustable electrode 2 = G6 suppressor
--
-- G1, G2, G5, collector, and upstream reference are fixed 0-V conductors
-- embedded in the PA geometry and are not independently fast-adjusted.
--
-- pa:fast_adjust() streams the basis solutions from disk and avoids the
-- large extra contiguous-memory allocation of the ordinary workbench path.
function segment.init_p_values()
  local vr = 0.0
  local vs = 0.0

  if not current_case.equipotential then
    vr = current_case.VR
    vs = suppressor_V
  end

  simion.wb.instances[1].pa:fast_adjust{
    [1]=vr,
    [2]=vs
  }
end

function segment.initialize_run()
  reset_run_state()
  simion.seed(random_seed)
  sim_trajectory_image_control = production_trajectory_image_control

  if sim_ions_count and sim_ions_count < minimum_recommended_particles then
    print(string.format(
      "WARNING: only %d particles are defined; >= %d recommended.",
      sim_ions_count, minimum_recommended_particles))
  end

  local v1,v2,v3,v4,v5,v6,v7,v8 = case_electrode_voltages()
  print("------------------------------------------------------------")
  print(string.format("CASE: %s", current_case.label))
  print(string.format(
    "Voltages: G1=%g G2=%g G3=%g G4=%g G5=%g G6=%g collector(model)=%g upstream=%g V",
    v1,v2,v3,v4,v5,v6,v7,v8))
  print("Collector is fixed at 0 V by the adopted virtual-ground amplifier design.")
end

function segment.initialize()
  n_initialized = n_initialized + 1
  ion_mass = ion_mass_amu
  ion_charge = ion_charge_e

  local x0,y0 = sample_launch_position()
  ion_px_mm = x0
  ion_py_mm = y0
  ion_pz_mm = launch_z_mm

  local g1,g2 = gaussian_pair()
  local g3,_  = gaussian_pair()

  ion_vx_mm = mean_vx_mm_us + sigma_v_mm_us*g1
  ion_vy_mm = mean_vy_mm_us + sigma_v_mm_us*g2
  ion_vz_mm = mean_vz_mm_us + sigma_v_mm_us*g3

  local w = math.max(ion_vz_mm * 1000.0, 0.0)
  weight_mps[ion_number] = w
  sum_w_all = sum_w_all + w

  if w > 0 then
    n_positive_flux = n_positive_flux + 1
  else
    ion_splat = 1
  end

  last_z[ion_number] = ion_pz_mm
  counted_collector[ion_number] = false
  changed_column[ion_number] = false
end

function segment.other_actions()
  local z0 = last_z[ion_number]
  local z1 = ion_pz_mm
  local w = weight_mps[ion_number] or 0.0

  if z0 ~= nil then
    for k=1,#diag_names do
      if not counted_diag[k][ion_number] and crossed_forward(z0,z1,diag_z[k]) then
        counted_diag[k][ion_number] = true
        diag_count[k] = diag_count[k] + 1
        diag_sum_w[k] = diag_sum_w[k] + w

        local ci,cj = nearest_hex_indices(ion_px_mm, ion_py_mm)
        if first_col_i[ion_number] == nil then
          first_col_i[ion_number] = ci
          first_col_j[ion_number] = cj
        elseif ci ~= first_col_i[ion_number] or cj ~= first_col_j[ion_number] then
          changed_column[ion_number] = true
        end
      end
    end

    if not counted_collector[ion_number]
       and crossed_forward(z0,z1,collector_witness_z_mm) then

      counted_collector[ion_number] = true
      collector_hits = collector_hits + 1
      sum_w_col = sum_w_col + w
      sum_w2_col = sum_w2_col + w*w

      if changed_column[ion_number] then
        neighbor_hits = neighbor_hits + 1
        sum_w_neighbor = sum_w_neighbor + w
      else
        samecol_hits = samecol_hits + 1
        sum_w_samecol = sum_w_samecol + w
      end
    end
  end

  last_z[ion_number] = z1
end

local function write_header(f)
  f:write("# CARLO 0.0025 mm/gu memory-safe high-resolution production sweep\n")
  f:write("# physical_case=O+ 1000K exact_4deg total_bulk_speed_7600mps density_1e5cm-3\n")
  f:write("# geometry=six_registered_grids d0.35mm p0.52mm t0.05mm h0.50mm G6_to_collector2.00mm\n")
  f:write("# final_geometry=3pitch_halfspan 0.0025mm_per_gu two_adjustable_groups\n")
  f:write("# G6=-15V for biased cases\n")
  f:write("# collector_V=0 V (adopted virtual-ground current-amplifier design)\n")
  f:write("# adjustable_1=G3+G4 retarding pair; adjustable_2=G6 suppressor\n")
  f:write(string.format("# random_seed=%d\n", random_seed))
  f:write(string.format("# accepted_fieldfree_ModelB_reference_nA=%.12g\n",
    accepted_modelB_fieldfree_nA))
  f:write("# current_adopted_legacy = current_raw * tau_stack/phi_hex; raw current is also retained\n")

  local cols = {
    "case_index","case_label","VR_V","equipotential",
    "G1_V","G2_V","G3_V","G4_V","G5_V","G6_V","collector_V","upstream_V",
    "N_initialized","N_positive_flux","N_collected",
    "sum_vnplus_started_mps","sum_vnplus_collected_mps",
    "unweighted_collection_fraction","flux_weighted_collection_fraction",
    "I_raw_cell_nA","I_adopted_nA","I_adopted_MC_SE_nA",
    "N_samecolumn_collected","N_neighborcolumn_collected",
    "I_adopted_samecolumn_nA","I_adopted_neighborcolumn_nA"
  }

  for k=1,#diag_names do
    table.insert(cols, diag_names[k].."_count")
    table.insert(cols, diag_names[k].."_weighted_fraction")
  end

  f:write(table.concat(cols,",") .. "\n")
end

local function write_case_row(f, case_index)
  local v1,v2,v3,v4,v5,v6,v7,v8 = case_electrode_voltages()

  local N = math.max(n_initialized,1)
  local unweighted_frac = collector_hits / N
  local weighted_frac = (sum_w_all > 0) and (sum_w_col/sum_w_all) or 0.0

  local mean_x = sum_w_col / N
  local Iraw_A = q_C * density_m3 * A_face_m2 * mean_x
  local Iadopt_A = Iraw_A * tau_stack / phi_hex

  local var_x = 0.0
  if N > 1 then
    var_x = (sum_w2_col - N*mean_x*mean_x) / (N-1)
    if var_x < 0 and math.abs(var_x) < 1e-12 then var_x = 0 end
    if var_x < 0 then var_x = 0 end
  end

  local se_mean_x = math.sqrt(var_x / N)
  local Iadopt_se_A = q_C * density_m3 * A_face_m2 *
                      se_mean_x * tau_stack / phi_hex

  local Isame_A = q_C * density_m3 * A_face_m2 *
                  (sum_w_samecol/N) * tau_stack / phi_hex
  local Ineigh_A = q_C * density_m3 * A_face_m2 *
                   (sum_w_neighbor/N) * tau_stack / phi_hex

  local row = {
    tostring(case_index),
    current_case.label,
    string.format("%.12g",current_case.VR),
    current_case.equipotential and "1" or "0",
    string.format("%.12g",v1), string.format("%.12g",v2),
    string.format("%.12g",v3), string.format("%.12g",v4),
    string.format("%.12g",v5), string.format("%.12g",v6),
    string.format("%.12g",v7), string.format("%.12g",v8),
    tostring(n_initialized), tostring(n_positive_flux), tostring(collector_hits),
    string.format("%.15g",sum_w_all),
    string.format("%.15g",sum_w_col),
    string.format("%.12g",unweighted_frac),
    string.format("%.12g",weighted_frac),
    string.format("%.12g",Iraw_A*1e9),
    string.format("%.12g",Iadopt_A*1e9),
    string.format("%.12g",Iadopt_se_A*1e9),
    tostring(samecol_hits), tostring(neighbor_hits),
    string.format("%.12g",Isame_A*1e9),
    string.format("%.12g",Ineigh_A*1e9)
  }

  for k=1,#diag_names do
    local wf = (sum_w_all > 0) and (diag_sum_w[k]/sum_w_all) or 0.0
    table.insert(row,tostring(diag_count[k]))
    table.insert(row,string.format("%.12g",wf))
  end

  f:write(table.concat(row,",") .. "\n")
  f:flush()

  print(string.format(
    "%s: collected %d/%d, weighted T=%.6f, I_adopted=%.6f +/- %.6f nA",
    current_case.label, collector_hits, n_initialized, weighted_frac,
    Iadopt_A*1e9, Iadopt_se_A*1e9))
end

function segment.flym()
  outfile = assert(io.open("carlo_FINAL_3p_0025_1000K_4deg_100k.csv","w"))
  write_header(outfile)

  for i=1,#cases do
    current_case = cases[i]
    run()
    write_case_row(outfile,i)
  end

  outfile:close()
  outfile = nil

  print("============================================================")
  print("CARLO sweep complete.")
  print("Results: carlo_FINAL_3p_0025_1000K_4deg_100k.csv")
  print("Collector fixed at 0 V: adopted virtual-ground design.")
  print("============================================================")
end
