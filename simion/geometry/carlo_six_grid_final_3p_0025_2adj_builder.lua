--[[
CARLO FINAL high-resolution six-grid builder — 3-pitch / 0.0025 mm/gu / 2 adjustable groups
===========================================================

Purpose
-------
Builds a fast-adjustable .PA# geometry for the CURRENT single-head model:

  G1   G2   G3   G4   G5   G6        collector
  |    |    |    |    |    |            |
  0.05-mm molybdenum foils, five 0.50-mm clear gaps
  circular registered holes, d = 0.35 mm on hexagonal pitch p = 0.52 mm
  G6 exit -> collector front = 2.00 mm

Optimized electrical grouping:
  FIXED 0-V conductors (non-adjustable in the PA#):
    G1, G2, G5, continuous collector, upstream reference plane

  adjustable electrode 1:
    G3 + G4 retarding pair -> +V_R

  adjustable electrode 2:
    G6 suppressor -> -15 V

The collector is now an author-adopted virtual-ground design and is
therefore fixed at spacecraft common (0 V).

This grouping is electrostatically identical to the locked nominal design
but requires only TWO fast-adjust basis solutions instead of eight.

IMPORTANT
---------
1. The adopted collector design is virtual ground at spacecraft common.
   G1, G2, G5, the collector, and the upstream reference boundary are fixed
   0-V conductors and do not need independent fast-adjust basis arrays.
2. The plate OUTER dimensions and collector lateral dimensions are not
   source-defined.  Therefore this script intentionally builds a LOCAL
   interior hexagonal-lattice patch, not a full sensor-head housing.
3. The hex lattice continues analytically to the PA edges.  Open PA side
   boundaries are SIMION Neumann boundaries.  Launch particles only in the
   central region and check lateral-domain convergence before publication.
4. The collector "thickness" below is only a computational downstream
   boundary slab.  Only the front-face location (2.00 mm behind G6 exit)
   is a hardware input.
5. Dimensions are evaluated in physical mm and then rasterized onto the
   selected SIMION grid.  Tighten gu_mm for a mesh-convergence check.

Usage
-----
From the SIMION main screen:
  Run Lua Program -> carlo_six_grid_geometry_builder.lua

The script writes:
  carlo_six_grid_local.pa#

Then Refine that .PA# to generate the .PA0/.PA1/... fast-adjust arrays.

This script requires SIMION 8.1+ (simion.pas).
--]]

-- ============================================================
-- Author-locked physical geometry (mm)
-- ============================================================

local aperture_d_mm       = 0.35
local aperture_r_mm       = aperture_d_mm / 2
local pitch_mm            = 0.52
local foil_t_mm           = 0.05
local clear_gap_mm        = 0.50
local n_grids             = 6
local g6_to_collector_mm  = 2.00

-- Derived: G1 entrance face to G6 exit face.
local L_col_mm = n_grids * foil_t_mm + (n_grids - 1) * clear_gap_mm
assert(math.abs(L_col_mm - 2.80) < 1e-12)

-- ============================================================
-- Numerical model controls (NOT hardware dimensions)
-- ============================================================

-- FINAL HIGH-RESOLUTION CASE: 0.0025 mm/gu.
-- Mesh-convergence testing showed this resolution is required for the
-- primary biased-current result.
--   pitch = 208 gu exactly
--   clear gap = 200 gu exactly
--   aperture diameter = 140 gu exactly
--   foil thickness = 20 gu exactly
--   G6-to-collector gap = 800 gu exactly
local gu_mm = 0.0025

-- Empty region in front of G1, used for particle launch / entrance field.
-- Numerical-domain choice, not hardware.
local upstream_mm = 1.00

-- FINAL laterally converged local patch.
-- 3-pitch half-span -> total transverse width ~3.12 mm.
-- Lateral-domain testing showed 3 -> 4 pitches changes the primary
-- response by << 1%, so 3 pitches is adopted.
local lateral_half_span_pitches = 3
local xy_half_mm = lateral_half_span_pitches * pitch_mm

-- Computational collector boundary slab thickness only.
-- Actual collector thickness is not asserted.
local collector_boundary_thickness_mm = 0.05

local output_name = "carlo_six_grid_final_3p_0025_2adj.pa#"

-- Include a one-grid-plane FIXED 0-V upstream reference boundary.
-- It is computational, not hardware, and is not fast-adjustable.
local include_upstream_reference_plane = true

-- ============================================================
-- Helpers
-- ============================================================

local function round_nearest(x)
  return math.floor(x + 0.5)
end

local function make_odd(n)
  if n % 2 == 0 then return n + 1 end
  return n
end

local function nearest_hex_center_distance2(x_mm, y_mm)
  -- Infinite registered hexagonal lattice:
  -- row spacing = sqrt(3)/2 * p
  -- alternate rows offset by p/2.
  local row_dy = math.sqrt(3) * pitch_mm / 2
  local j0 = round_nearest(y_mm / row_dy)
  local best = math.huge

  -- Search nearby rows/columns; enough to find the nearest center.
  for j = j0 - 2, j0 + 2 do
    local offset = 0
    if (j % 2) ~= 0 then offset = pitch_mm / 2 end

    local i0 = round_nearest((x_mm - offset) / pitch_mm)
    for i = i0 - 2, i0 + 2 do
      local xc = i * pitch_mm + offset
      local yc = j * row_dy
      local dx = x_mm - xc
      local dy = y_mm - yc
      local d2 = dx*dx + dy*dy
      if d2 < best then best = d2 end
    end
  end

  return best
end

local function is_aperture(x_mm, y_mm)
  return nearest_hex_center_distance2(x_mm, y_mm) <= aperture_r_mm^2
end

local function grid_z_start(k)
  -- k = 1..6
  return upstream_mm + (k - 1) * (foil_t_mm + clear_gap_mm)
end

local g6_exit_mm = grid_z_start(6) + foil_t_mm
local collector_front_mm = g6_exit_mm + g6_to_collector_mm
local zmax_mm = collector_front_mm + collector_boundary_thickness_mm

-- SIMION array dimensions.
-- Center x,y on the central aperture at x=y=0.
local nx = make_odd(math.ceil((2 * xy_half_mm) / gu_mm) + 1)
local ny = nx
local nz = math.ceil(zmax_mm / gu_mm) + 1

local xcenter_i = (nx - 1) / 2
local ycenter_i = (ny - 1) / 2

local function x_mm_of(i) return (i - xcenter_i) * gu_mm end
local function y_mm_of(j) return (j - ycenter_i) * gu_mm end
local function z_mm_of(k) return k * gu_mm end

local function point_in_closed_interval(z, a, b)
  -- Physical-coordinate inclusion. Small epsilon avoids binary roundoff.
  local eps = 1e-12
  return z >= a - eps and z <= b + eps
end

-- ============================================================
-- Create PA#
-- ============================================================

assert(simion and simion.pas,
  "This builder requires SIMION 8.1+ with simion.pas support.")

print("============================================================")
print("FINAL MODEL: 3-pitch half-span, 0.0025 mm/gu, two adjustable voltage groups")
print("Fixed 0 V: G1, G2, G5, collector, upstream reference")
print("Adj 1: G3+G4 = +VR; Adj 2: G6 = -15 V")
print("Building CARLO six-grid local geometry")
print(string.format("Output: %s", output_name))
print(string.format("Mesh: %.6f mm/gu", gu_mm))
print(string.format("PA size: %d x %d x %d grid points", nx, ny, nz))
print(string.format("Nominal transverse model width: %.6f mm", (nx-1)*gu_mm))
print(string.format("Nominal z extent: %.6f mm", (nz-1)*gu_mm))
print(string.format("G1 entrance z: %.6f mm", upstream_mm))
print(string.format("G6 exit z: %.6f mm", g6_exit_mm))
print(string.format("Collector front z: %.6f mm", collector_front_mm))
print("============================================================")

local pa = simion.pas:open()
-- A newly opened PA starts with nz=1. Resize it first, then switch to
-- 3-D planar symmetry; SIMION rejects 3dplanar while nz is still 1.
pa:size(nx, ny, nz)
pa.symmetry = '3dplanar'
pa.dx_mm = gu_mm
pa.dy_mm = gu_mm
pa.dz_mm = gu_mm
-- No pa.name field exists in simion.pas. The filename is assigned by pa:save().

-- Precompute the 2-D aperture mask.
print("Computing registered hexagonal aperture mask...")
local open_xy = {}
for i = 0, nx-1 do
  open_xy[i] = {}
  local x = x_mm_of(i)
  for j = 0, ny-1 do
    local y = y_mm_of(j)
    open_xy[i][j] = is_aperture(x, y)
  end
end

-- Fixed 0-V computational upstream reference plane.
-- A 0-V electrode in a PA# is non-adjustable; no basis array is required.
if include_upstream_reference_plane then
  print("Writing fixed 0-V computational upstream reference plane...")
  local k = 0
  for i = 0, nx-1 do
    for j = 0, ny-1 do
      pa:point(i, j, k, 0.0, true)
    end
  end
end

-- Six perforated molybdenum grids.
-- Fixed 0-V grids: G1, G2, G5
-- Adjustable electrode 1: G3 + G4
-- Adjustable electrode 2: G6
for eg = 1, n_grids do
  local z0 = grid_z_start(eg)
  local z1 = z0 + foil_t_mm

  local marker_value
  local marker_label
  if eg == 3 or eg == 4 then
    marker_value = 1
    marker_label = "adjustable 1 (+VR)"
  elseif eg == 6 then
    marker_value = 2
    marker_label = "adjustable 2 (-15 V)"
  else
    marker_value = 0.0
    marker_label = "fixed 0 V"
  end

  print(string.format(
    "Writing G%d [%s]: z = %.6f .. %.6f mm",
    eg, marker_label, z0, z1))

  for k = 0, nz-1 do
    local z = z_mm_of(k)
    if point_in_closed_interval(z, z0, z1) then
      for i = 0, nx-1 do
        for j = 0, ny-1 do
          if not open_xy[i][j] then
            pa:point(i, j, k, marker_value, true)
          end
        end
      end
    end
  end
end

-- Continuous collector boundary, fixed at virtual ground (0 V).
print(string.format(
  "Writing fixed 0-V collector: front z = %.6f mm",
  collector_front_mm))

for k = 0, nz-1 do
  local z = z_mm_of(k)
  if z >= collector_front_mm - 1e-12 then
    for i = 0, nx-1 do
      for j = 0, ny-1 do
        pa:point(i, j, k, 0.0, true)
      end
    end
  end
end

-- Save as fast-adjust definition array.
print("Saving fast-adjust geometry...")
pa:save(output_name)

-- ============================================================
-- Geometry report / sanity checks
-- ============================================================

local row_dy_mm = math.sqrt(3) * pitch_mm / 2
local phi_hex = math.pi * aperture_d_mm^2 /
                (2 * math.sqrt(3) * pitch_mm^2)

print("")
print("Geometry complete.")
print(string.format("  aperture diameter       = %.6f mm", aperture_d_mm))
print(string.format("  hex pitch               = %.6f mm", pitch_mm))
print(string.format("  hex row spacing         = %.9f mm", row_dy_mm))
print(string.format("  foil thickness          = %.6f mm", foil_t_mm))
print(string.format("  clear inter-grid gap    = %.6f mm", clear_gap_mm))
print(string.format("  G1 entrance -> G6 exit  = %.6f mm", L_col_mm))
print(string.format("  G6 exit -> collector    = %.6f mm", g6_to_collector_mm))
print(string.format("  single-plane open frac  = %.9f", phi_hex))
print("")
print("Electrical grouping in PA#:")
print("  fixed 0 V: G1, G2, G5, collector, upstream reference")
print("  adjustable 1: G3 + G4 retarding pair")
print("  adjustable 2: G6 suppressor")
print("")
print("Nominal fast-adjust map:")
print("  adjustable 1 = +VR")
print("  adjustable 2 = -15 V")
print("  all fixed conductors remain at 0 V")
print("")
print("Refine should generate only PA1 and PA2 basis arrays (plus PA0).")
print("")
print("")
print("IMPORTANT MODEL-SCOPE NOTE:")
print("  This is the exact locked aperture/axial geometry in a LOCAL")
print("  interior lattice patch. Plate outer size, border, housing, and")
print("  collector lateral size are not source-defined and are not claimed.")
print("  Check central-region field stability by increasing")
print("  lateral_half_span_pitches and by tightening gu_mm.")
print("")
print("Next step in SIMION: Refine the saved .PA# file.")
print("============================================================")
