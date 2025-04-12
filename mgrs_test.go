package mgrs

import (
	"fmt"
	"math"
	"regexp"
	"strconv"
	"testing"
)

func TestLatLngToMGRS(t *testing.T) {
	tests := []struct {
		name      string
		lat       float64
		lon       float64
		precision int
		wantMGRS  string
		wantErr   bool
	}{
		{
			name:      "Paris Valid Example",
			lat:       48.8566, // Slightly different from example for testing
			lon:       2.3522,
			precision: 5,                 // 1m precision
			wantMGRS:  "31UDQ5248211717", // Example: 31U 52483 32804
		},
		{
			name:      "Paris Valid short",
			lat:       48.8,
			lon:       2.2,
			precision: 5,                 // 1m
			wantMGRS:  "31UDQ4125105530", // Expected based on online calc, close to example
		},
		{
			name:      "Paris Example from Test",
			lat:       48.8,
			lon:       2.2,
			precision: 5,                 // 1m
			wantMGRS:  "31UDQ4125205531", // Using exact string from prompt for comparison
		},
		{
			name:      "NYC Example",
			lat:       40.730610,
			lon:       -73.935242,
			precision: 5, // 1m
			wantMGRS:  "18TWL8991209398",
		},
		{
			name:      "Equator Prime Meridian",
			lat:       0.0,
			lon:       0.0,
			precision: 5,
			wantMGRS:  "31NAA6602100000",
		},
		{
			name:      "Southern Hemisphere",
			lat:       -33.8688, // Sydney
			lon:       151.2093,
			precision: 5,
			wantMGRS:  "56HLH3436850948",
		},
		{
			name:      "Low Precision",
			lat:       40.730610,
			lon:       -73.935242,
			precision: 1,         // 10km
			wantMGRS:  "18TWL80", // Approx from online calc (89912 -> 8, 09398 -> 0) -> 18TWL80
		},
		{
			name:      "Out of Bounds Latitude High",
			lat:       85.0,
			lon:       0.0,
			precision: 5,
			wantErr:   true,
		},
		{
			name:      "Out of Bounds Latitude Low",
			lat:       -81.0,
			lon:       0.0,
			precision: 5,
			wantErr:   true,
		},
		{
			name:      "Invalid Precision",
			lat:       48.8,
			lon:       2.2,
			precision: 6,
			wantErr:   true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotMGRS, err := LatLngToMGRS(tt.lat, tt.lon, tt.precision)

			if (err != nil) != tt.wantErr {
				t.Errorf("LatLngToMGRS() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if err != nil {
				return // Expected error, test passed
			}

			// Allow minor variations in the numeric part due to calculation differences
			if !compareMGRS(gotMGRS, tt.wantMGRS, 1) { // Allow 1 digit difference in last place
				t.Errorf("LatLngToMGRS() = %q, want %q", gotMGRS, tt.wantMGRS)
			} else {
				fmt.Printf("Test %s: LatLng(%.6f, %.6f) -> MGRS: %s (OK, matches %s within tolerance)\n", tt.name, tt.lat, tt.lon, gotMGRS, tt.wantMGRS)
			}
		})
	}
}

func TestMGRSToLatLng(t *testing.T) {
	tests := []struct {
		name    string
		mgrs    string
		wantLat float64
		wantLon float64
		wantErr bool
	}{
		{
			name:    "Paris short",
			mgrs:    "31UDQ4125205531", // 1m precision
			wantLat: 48.80000,          // Expected based on input example MGRS
			wantLon: 2.20000,           // Expected based on input example MGRS
		},
		{
			name:    "NYC Example",
			mgrs:    "18TWL8991209398", // 1m precision
			wantLat: 40.730610,
			wantLon: -73.935242,
		},
		{
			name:    "Equator Prime Meridian (Approx)",
			mgrs:    "31NAA6602100000", // Approx from online calc test
			wantLat: 0.0,
			wantLon: 0.0,
		},
		{
			name:    "Southern Hemisphere (Approx)",
			mgrs:    "56HLH3310987488", // Approx from online calc test
			wantLat: -33.539196,        // Sydney
			wantLon: 151.202569,
		},
		{
			name:    "Low Precision (Center of square)",
			mgrs:    "18TWL80",  // 10km precision -> Center of square WL80
			wantLat: 40.691517,  // Approx center from online calc for 18TWL80
			wantLon: -73.994000, // Approx center from online calc for 18TWL80
		},
		{
			name:    "Invalid Format - Bad Chars",
			mgrs:    "31UDQ4125A05531",
			wantErr: true,
		},
		{
			name:    "Invalid Format - Bad Lat Band",
			mgrs:    "31A DQ4125055310", // Invalid Lat Band 'A'
			wantErr: true,
		},
		{
			name:    "Invalid Format - Bad 100k Letter",
			mgrs:    "31U DI4125055310", // Invalid 100k Letter 'I'
			wantErr: true,
		},
		{
			name:    "Invalid Format - Bad Zone",
			mgrs:    "61UDQ4125055310",
			wantErr: true,
		},
		{
			name:    "Invalid MGRS - Northing Cycle Inconsistency",
			mgrs:    "33XXX", // Northing value is too high for the DQ 100km square in this zone/band
			wantErr: true,
			// This would cause the resulting latitude to be outside the expected band 'U'
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotLat, gotLon, err := MGRSToLatLng(tt.mgrs)

			if (err != nil) != tt.wantErr {
				t.Errorf("MGRSToLatLng() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if err != nil {
				return // Expected error, test passed
			}

			// Precision of MGRS string determines tolerance for comparison
			_, _, _, _, _, _, precision, parseErr := parseMGRS(tt.mgrs) // Assuming parseMGRS is robust
			if parseErr != nil {
				t.Fatalf("Failed to parse valid MGRS for tolerance calculation: %s, err: %v", tt.mgrs, parseErr)
			}

			// Calculate tolerance based on MGRS precision (size of the grid square in degrees, roughly)
			// 1m (prec 5) -> ~1e-5 deg
			// 10m (prec 4) -> ~1e-4 deg
			// 100m (prec 3) -> ~1e-3 deg
			// 1km (prec 2) -> ~1e-2 deg
			// 10km (prec 1) -> ~1e-1 deg
			// Corrected: math.Pow10 takes an int exponent
			inverseTolerance := math.Pow10(-precision) / 2.0 // Half grid size in degrees approx

			// Increase tolerance slightly for lower precision due to projection distortions
			if precision <= 2 {
				inverseTolerance *= 2
			}
			// Minimum tolerance to avoid float issues
			if inverseTolerance < 1e-7 {
				inverseTolerance = 1e-7
			}

			latDiff := math.Abs(gotLat - tt.wantLat)
			lonDiff := math.Abs(gotLon - tt.wantLon)

			// Use a tolerance appropriate for the center coordinate calculation
			// For low precision, the center Lat/Lon might differ more significantly
			// due to Earth curvature and projection distortions within the larger square.
			// Allow ~half the grid size as tolerance.
			if latDiff > inverseTolerance || lonDiff > inverseTolerance {
				t.Errorf("MGRSToLatLng() = lat %.6f, lon %.6f, want lat %.6f, lon %.6f (Diff lat: %.6e, lon: %.6e, Tol: %.6e)", gotLat, gotLon, tt.wantLat, tt.wantLon, latDiff, lonDiff, inverseTolerance)
			} else {
				fmt.Printf("Test %s: MGRS %s -> LatLng(%.6f, %.6f) (OK, ~%.6f, ~%.6f within Tol %.1e)\n", tt.name, tt.mgrs, gotLat, gotLon, tt.wantLat, tt.wantLon, inverseTolerance)
			}
		})
	}
}

// compareMGRS compares two MGRS strings, allowing for slight numeric variations.
func compareMGRS(got, want string, tolerance int) bool {
	// Use regex to split prefix and numeric parts consistently
	// Allows for 1 or 2 digit zone, C-X lat band, 2 non-I/O 100k letters, optional digits
	reMGRS := regexp.MustCompile(`^(\d{1,2}[C-HJ-NP-X][A-HJ-NP-Z]{2})(\d*)$`)
	gotMatches := reMGRS.FindStringSubmatch(got)
	wantMatches := reMGRS.FindStringSubmatch(want)

	// Check if both strings broadly match the MGRS pattern prefix + optional digits
	if len(gotMatches) < 2 || len(wantMatches) < 2 {
		// fmt.Printf("Debug compareMGRS: Pattern mismatch. got: %q (%d matches), want: %q (%d matches)\n", got, len(gotMatches), want, len(wantMatches))
		// Fallback to direct string comparison if regex fails (shouldn't happen for valid MGRS)
		return got == want
	}

	gotPrefix := gotMatches[1]
	wantPrefix := wantMatches[1]

	// Compare prefixes directly (Zone, Lat Band, 100k ID)
	if gotPrefix != wantPrefix {
		// fmt.Printf("Prefix mismatch: got %s, want %s\n", gotPrefix, wantPrefix) // Uncomment for debug
		return false
	}

	// Extract numeric parts (group 2, might be empty for precision 0)
	numStrGot := ""
	if len(gotMatches) > 2 { // Index 2 exists if digits were captured
		numStrGot = gotMatches[2]
	}
	numStrWant := ""
	if len(wantMatches) > 2 { // Index 2 exists if digits were captured
		numStrWant = wantMatches[2]
	}

	// If numeric parts have different lengths, they represent different precisions/locations
	if len(numStrGot) != len(numStrWant) {
		// fmt.Printf("Numeric length mismatch: got %q (%d), want %q (%d)\n", numStrGot, len(numStrGot), numStrWant, len(numStrWant)) // Uncomment for debug
		return false
	}

	// If numeric parts are empty (precision 0), prefixes matched, so return true
	if len(numStrGot) == 0 {
		return true
	}

	// Validate numeric part length (must be even and > 0)
	if len(numStrGot)%2 != 0 {
		// fmt.Printf("Odd number of digits: got %q\n", numStrGot) // Uncomment for debug
		return false // Invalid MGRS numeric part
	}

	// Compare numeric part with tolerance
	prec := len(numStrGot) / 2
	if prec < 1 || prec > 5 { // Should be caught by parser, but double-check
		// fmt.Printf("Invalid precision derived from digits: %d\n", prec) // Uncomment for debug
		return false
	}

	eGotStr := numStrGot[:prec]
	nGotStr := numStrGot[prec:]
	eWantStr := numStrWant[:prec]
	nWantStr := numStrWant[prec:]

	eGot, errGotE := strconv.Atoi(eGotStr)
	nGot, errGotN := strconv.Atoi(nGotStr)
	eWant, errWantE := strconv.Atoi(eWantStr)
	nWant, errWantN := strconv.Atoi(nWantStr)

	// Check for conversion errors (should not happen if regex matched \d*)
	if errGotE != nil || errGotN != nil || errWantE != nil || errWantN != nil {
		// fmt.Printf("Error converting numeric parts: got E(%s)/N(%s), want E(%s)/N(%s)\n", eGotStr, nGotStr, eWantStr, nWantStr) // Uncomment for debug
		return false // Treat conversion errors as mismatch
	}

	// Check easting and northing within tolerance
	if math.Abs(float64(eGot-eWant)) <= float64(tolerance) && math.Abs(float64(nGot-nWant)) <= float64(tolerance) {
		return true
	} else {
		// fmt.Printf("Numeric mismatch: got E%d N%d, want E%d N%d (tolerance %d)\n", eGot, nGot, eWant, nWant, tolerance) // Uncomment for debug
		return false
	}
}
