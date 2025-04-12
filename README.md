## MGRS

A Military Grid Reference System (MGRS) computation library for Go

## Usage

```Go
// Convert latitude/longitude to MGRS
	lat, lon := 40.7306, -73.9352 // New York City
	precision := 5                // 1-meter precision (precision can be 1-5, with 5 being the highest)

	mgrsString, _ := mgrs.LatLngToMGRS(lat, lon, precision)
	fmt.Printf("Coordinates (%.4f, %.4f) in MGRS: %s\n", lat, lon, mgrsString)

	// Convert MGRS back to latitude/longitude
	lat2, lon2, _ := mgrs.MGRSToLatLng(mgrsString)
	fmt.Printf("MGRS %s as lat/lng: (%.4f, %.4f)\n", mgrsString, lat2, lon2)

	// Using lower precision (1 - 10km)
	lowPrecision := 1
	mgrsStringLow, _ := mgrs.LatLngToMGRS(lat, lon, lowPrecision)
	fmt.Printf("Coordinates (%.4f, %.4f) in MGRS (10km precision): %s\n", lat, lon, mgrsStringLow)
```

## Authors

Becaude this project was mostly generated by Gemini 2.5 pro, the code is attributed to all developers who published public code before 2025.
