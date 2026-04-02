# tests/testthat/test_formatData.R

# helper to create a clean base dataset
get_base_data <- function() {
  data.frame(
    id = rep("A1", 5),
    date = as.POSIXct(c("2023-01-01 10:00:00", "2023-01-01 12:00:00",
                        "2023-01-01 14:00:00", "2023-01-01 16:00:00",
                        "2023-01-01 18:00:00"), tz = "UTC"),
    x = c(500000, 501000, 502000, 503000, 504000), # Metric Easting
    y = c(4000000, 4001000, 4002000, 4003000, 4004000), # Metric Northing
    lc = c("3", "2", "G", "A", "B"),
    smaj = c(10, 15, NA, NA, NA),
    smin = c(5, 7, NA, NA, NA),
    eor  = c(45, 90, NA, NA, NA), # degrees from north
    x.sd = c(NA, NA, 5, NA, NA),
    y.sd = c(NA, NA, 5, NA, NA),
    extra_var = c("temp1", "temp2", "temp3", "temp4", "temp5")
  )
}

test_that("Standard POSIXct + Mixed Errors format correctly", {
  dat1 <- get_base_data()
  res1 <- formatData(dat1)

  expect_equal(nrow(res1), 5)
  expect_true("extra_var" %in% names(res1))
  expect_s3_class(res1, "dataLangevin")
})

test_that("Character dates and custom coordinates are handled safely", {
  dat2 <- get_base_data()
  dat2$date <- as.character(dat2$date) # Strip POSIXct
  names(dat2)[3:4] <- c("longitude", "latitude") # Custom coords

  res2 <- formatData(dat2, coord = c("longitude", "latitude"))

  expect_s3_class(res2$date, "POSIXt")
  expect_true(all(c("x", "y") %in% names(res2)))
})

test_that("Custom error parameter names are successfully standardized", {
  dat3 <- get_base_data()
  names(dat3)[6:10] <- c("err_maj", "err_min", "err_angle", "err_x", "err_y")

  res3 <- formatData(dat3,
                     epar = c("err_maj", "err_min", "err_angle"),
                     sderr = c("err_x", "err_y"))

  expect_true(all(c("smaj", "smin", "eor", "x.sd", "y.sd") %in% names(res3)))
})

test_that("EMF integration accurately fills missing NA errors", {
  dat5 <- get_base_data()
  my_emf <- get_emf()

  # ensure it is NA before formatting
  expect_true(is.na(dat5$x.sd[4]))

  res5 <- formatData(dat5, emf = my_emf)

  # ensure it was filled correctly after formatting
  expect_false(is.na(res5$x.sd[4]))
  expect_equal(res5$x.sd[4], my_emf$emf.x[my_emf$lc == "A"])
})

test_that("sf spatial objects are safely parsed and stripped", {
  dat6 <- get_base_data()
  # Use crs = 32611 (UTM Zone 11N) which is a projected coordinate system, bypassing the lon/lat error
  dat_sf <- sf::st_as_sf(dat6, coords = c("x", "y"), crs = 32611)

  res6 <- formatData(dat_sf)

  expect_true(all(c("x", "y") %in% names(res6)))
  expect_false(inherits(res6, "sf"))
})

test_that("Radians vs Degrees are detected and warned appropriately", {
  # test safe degrees (should not warn)
  dat7 <- get_base_data()
  expect_silent(res7 <- formatData(dat7))

  # test radians (should trigger warning)
  dat7_rads <- get_base_data()
  dat7_rads$eor <- c(0.5, 1.2, NA, NA, NA) # < pi

  expect_warning(formatData(dat7_rads), "appear to have been provided in radians rather than degrees")
})

test_that("Explicit lc='G' preserves user-specified ellipses or SDs", {
  # Create a dataset where every row is explicitly set to "G" but with different errors
  dat8 <- data.frame(
    id = rep("A1", 3),
    date = as.POSIXct(c("2023-01-01 10:00:00", "2023-01-01 12:00:00", "2023-01-01 14:00:00"), tz = "UTC"),
    x = c(500000, 501000, 502000), # Metric bounds
    y = c(4000000, 4001000, 4002000),
    lc = rep("G", 3),
    smaj = c(NA, 10, NA), # Row 2 has an ellipse
    smin = c(NA, 5, NA),
    eor  = c(NA, 45, NA),
    x.sd = c(5, NA, NA),  # Row 1 has SDs
    y.sd = c(5, NA, NA)
  ) # Row 3 has no errors

  res8 <- formatData(dat8)

  # row 1: G with SDs
  expect_equal(res8$x.sd[1], 5)
  expect_true(is.na(res8$smaj[1]))

  # row 2: G with Ellipse (ensure eor converted to radians correctly)
  expect_equal(res8$smaj[2], 10)
  expect_equal(res8$eor[2], 45 * pi / 180)
  expect_true(is.na(res8$x.sd[2]))

  # row 3: G with no errors
  expect_true(is.na(res8$x.sd[3]))
  expect_true(is.na(res8$smaj[3]))

  # all rows should still be class "G"
  expect_true(all(res8$lc == "G"))
})

test_that("Argos classes preserve user-specified errors and selectively bypass EMF", {
  # Create a dataset with Argos classes, but pre-filled with different user errors
  dat9 <- data.frame(
    id = rep("A1", 3),
    date = as.POSIXct(c("2023-01-01 10:00:00", "2023-01-01 12:00:00", "2023-01-01 14:00:00"), tz = "UTC"),
    x = c(500000, 501000, 502000), # Metric bounds
    y = c(4000000, 4001000, 4002000),
    lc = c("3", "A", "B"),
    smaj = c(100, NA, NA), # Row 1 (LC 3) has an explicit ellipse
    smin = c(50, NA, NA),
    eor  = c(45, NA, NA),
    x.sd = c(NA, 75, NA),  # Row 2 (LC A) has explicit SDs
    y.sd = c(NA, 75, NA)   # Row 3 (LC B) has NO errors
  )

  my_emf <- get_emf()
  res9 <- formatData(dat9, emf = my_emf)

  # Row 1 (LC 3): Should KEEP its user-provided ellipse and NOT get EMF SDs
  expect_equal(res9$smaj[1], 100)
  expect_true(is.na(res9$x.sd[1]))

  # Row 2 (LC A): Should KEEP its user-provided SDs (75) and NOT get overwritten by EMF
  expect_equal(res9$x.sd[2], 75)
  expect_true(res9$x.sd[2] != my_emf$emf.x[my_emf$lc == "A"])

  # Row 3 (LC B): Had no errors, so it SHOULD be filled by the EMF table
  expect_false(is.na(res9$x.sd[3]))
  expect_equal(res9$x.sd[3], my_emf$emf.x[my_emf$lc == "B"])
})

test_that("Missing 'lc' column is assigned NA without guessing", {
  dat10 <- get_base_data()
  dat10$lc <- NULL # Completely omit the location class column

  res10 <- formatData(dat10)

  expect_true("lc" %in% names(res10))
  expect_true(all(is.na(res10$lc)))
})

test_that("Invalid 'lc' classes throw an error", {
  dat11 <- get_base_data()
  dat11$lc[1] <- "INVALID"

  expect_error(formatData(dat11), "Invalid location classes detected")
})

test_that("predTimes correctly pads data, sorts, injects NAs, and recalculates dt", {
  dat12 <- get_base_data()

  # insert predictions between the 1st/2nd and 3rd/4th observations, and one at the end
  pt <- data.frame(
    id = "A1",
    date = as.POSIXct(c("2023-01-01 11:00:00",
                        "2023-01-01 15:00:00",
                        "2023-01-01 20:00:00"), tz = "UTC")
  )

  res12 <- formatData(dat12, predTimes = pt)

  # check total rows (5 original + 3 predicted = 8)
  expect_equal(nrow(res12), 8)

  # check temporal sorting (10:00, 11:00, 12:00...)
  expect_equal(res12$date[2], pt$date[1])
  expect_equal(res12$date[5], pt$date[2])
  expect_equal(res12$date[8], pt$date[3])

  # check NA injection for coordinates and observation errors at the predicted times
  expect_true(all(is.na(res12$x[c(2, 5, 8)])))
  expect_true(all(is.na(res12$y[c(2, 5, 8)])))
  expect_true(all(is.na(res12$lc[c(2, 5, 8)])))
  expect_true(all(is.na(res12$smaj[c(2, 5, 8)])))

  # check that dt is recalculated correctly (time.unit = "hours" by default)
  # 10:00 to 11:00 is 1 hr. 11:00 to 12:00 is 1 hr.
  expect_equal(res12$dt[2], 1)
  expect_equal(res12$dt[3], 1)
})

test_that("predTimes catches input formatting errors", {
  dat13 <- get_base_data()

  # not a data frame
  expect_error(formatData(dat13, predTimes = list(id="A1", date=dat13$date[1])),
               "must be a data frame")

  # missing required columns
  pt_bad_cols <- data.frame(Animal="A1", Time=dat13$date[1])
  expect_error(formatData(dat13, predTimes = pt_bad_cols),
               "must contain")

  # id not present in the actual dataset
  pt_bad_id <- data.frame(id = "GhostAnimal", date = dat13$date[2])
  expect_error(formatData(dat13, predTimes = pt_bad_id),
               "must exist in")

  # mismatched date classes (data is POSIXct, predTime is Numeric)
  pt_bad_class <- data.frame(id = "A1", date = 1672570800)
  expect_error(formatData(dat13, predTimes = pt_bad_class),
               "must be of the same type")

  # dates occurring before the first real observation of a track
  pt_early <- data.frame(id = "A1", date = as.POSIXct("2023-01-01 08:00:00", tz="UTC"))
  expect_error(formatData(dat13, predTimes = pt_early),
               "prior to its first observation")
})

test_that("predTimes works with custom user-defined column names", {
  dat14 <- get_base_data()
  # User renames columns to non-standard names
  names(dat14)[1:2] <- c("AnimalID", "Timestamp")

  pt <- data.frame(
    AnimalID = "A1",
    Timestamp = as.POSIXct("2023-01-01 11:00:00", tz="UTC")
  )

  # Pass the custom column names into formatData alongside the matching predTimes
  res14 <- formatData(dat14, id = "AnimalID", date = "Timestamp", predTimes = pt)

  # Verify it processed, appended, and standardized the names correctly
  expect_equal(nrow(res14), 6)
  expect_true("id" %in% names(res14))
  expect_true("date" %in% names(res14))

  # Verify the 11:00 inserted row (Row 2) got its NAs
  expect_true(is.na(res14$x[2]))
})

test_that("predTimes works when dates are provided as character strings", {
  dat15 <- get_base_data()
  # Convert both data and predTimes dates to raw character strings
  dat15$date <- as.character(dat15$date)

  pt <- data.frame(
    id = "A1",
    date = "2023-01-01 11:00:00"
  )

  # Should not trigger the mismatched class error or the "prior to first obs" error
  res15 <- formatData(dat15, predTimes = pt)

  expect_equal(nrow(res15), 6)
  expect_s3_class(res15$date, "POSIXt") # Ensures formatData still properly converted them internally
})

test_that("NA coordinates in the first observation of any track throws an error", {
  dat16 <- data.frame(
    id = c(rep("A1", 2), rep("A2", 2)),
    date = as.POSIXct(c("2023-01-01 10:00:00", "2023-01-01 12:00:00",
                        "2023-01-01 10:00:00", "2023-01-01 12:00:00"), tz = "UTC"),
    x = c(500000, 501000, NA, 503000),
    y = c(4000000, 4001000, NA, 4003000),
    lc = rep("G", 4)
  )

  expect_error(formatData(dat16), "The first observation for each track in 'data' cannot have missing \\(NA\\) coordinates.")
})
