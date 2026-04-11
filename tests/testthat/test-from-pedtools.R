## Tests for from_pedtools() — the paramlink-free entry point.

test_that("from_pedtools returns the expected list structure", {
  skip_if_not_installed("pedtools")
  ped <- pedtools::linearPed(2)
  m   <- pedtools::marker(ped, afreq = c(A = 0.3, B = 0.4, C = 0.3), name = "M1")
  ped <- pedtools::setMarkers(ped, m)

  pbn <- from_pedtools(ped)
  expect_type(pbn, "list")
  expect_named(pbn, c("ped", "markerEvidence", "markerLinkage", "alelFreq", "mmodel"))
  expect_true(all(c("pedigree", "ped", "orig.ids", "markerdata", "available", "nMark") %in% names(pbn$ped)))
})

test_that("from_pedtools records the pedigree topology correctly", {
  skip_if_not_installed("pedtools")
  ped <- pedtools::linearPed(2)
  m   <- pedtools::marker(ped, afreq = c(A = 0.5, B = 0.5), name = "M1")
  ped <- pedtools::setMarkers(ped, m)

  pbn <- from_pedtools(ped)
  expect_equal(pbn$ped$orig.ids, as.character(labels(ped)))
  expect_equal(nrow(pbn$ped$pedigree), length(labels(ped)))
  expect_true(all(c("ID", "FID", "MID", "SEX", "AFF") %in% colnames(pbn$ped$pedigree)))
  # Founders should have FID = "0" and MID = "0"
  for (id in c("1", "2", "4")) {
    row <- which(pbn$ped$pedigree[, "ID"] == id)
    if (length(row) == 1) {
      expect_equal(unname(pbn$ped$pedigree[row, "FID"]), "0")
      expect_equal(unname(pbn$ped$pedigree[row, "MID"]), "0")
    }
  }
})

test_that("from_pedtools extracts allele frequencies per marker", {
  skip_if_not_installed("pedtools")
  ped <- pedtools::linearPed(2)
  m1 <- pedtools::marker(ped, afreq = c(A = 0.3, B = 0.7), name = "LocusA")
  m2 <- pedtools::marker(ped, afreq = c(x = 0.1, y = 0.2, z = 0.7), name = "LocusB")
  ped <- pedtools::setMarkers(ped, list(m1, m2))

  pbn <- from_pedtools(ped)
  expect_equal(length(pbn$alelFreq), 2)
  expect_equal(names(pbn$alelFreq), c("LocusA", "LocusB"))
  expect_equal(as.numeric(pbn$alelFreq$LocusA), c(0.3, 0.7))
  expect_equal(names(pbn$alelFreq$LocusA), c("A", "B"))
  expect_equal(as.numeric(pbn$alelFreq$LocusB), c(0.1, 0.2, 0.7))
})

test_that("from_pedtools detects typed individuals automatically", {
  skip_if_not_installed("pedtools")
  ped <- pedtools::linearPed(2)
  m <- pedtools::marker(ped, afreq = c(A = 0.5, B = 0.5), name = "M1")
  ped <- pedtools::setMarkers(ped, m)
  ped <- pedtools::setGenotype(ped, id = "1", marker = 1, geno = "A/B")

  pbn <- from_pedtools(ped)
  expect_true("1" %in% pbn$ped$available)
})

test_that("from_pedtools respects explicit `available` argument", {
  skip_if_not_installed("pedtools")
  ped <- pedtools::linearPed(2)
  m <- pedtools::marker(ped, afreq = c(A = 0.5, B = 0.5), name = "M1")
  ped <- pedtools::setMarkers(ped, m)

  pbn <- from_pedtools(ped, available = c("3"))
  expect_equal(pbn$ped$available, "3")
})

test_that("from_pedtools QP argument marks the pedigree AFF column", {
  skip_if_not_installed("pedtools")
  ped <- pedtools::linearPed(2)
  m <- pedtools::marker(ped, afreq = c(A = 0.5, B = 0.5), name = "M1")
  ped <- pedtools::setMarkers(ped, m)

  pbn <- from_pedtools(ped, QP = "5")
  qp_row <- which(pbn$ped$ped[, "ID"] == "5")
  expect_equal(unname(pbn$ped$ped[qp_row, "AFF"]), "2")
})

test_that("from_pedtools errors on non-ped input", {
  expect_error(from_pedtools("not a ped"), "pedtools::ped object")
  expect_error(from_pedtools(list(foo = 1)), "pedtools::ped object")
})

test_that("from_pedtools markerEvidence has the right column names", {
  skip_if_not_installed("pedtools")
  ped <- pedtools::linearPed(2)
  m1 <- pedtools::marker(ped, afreq = c(A = 0.5, B = 0.5), name = "LocusA")
  m2 <- pedtools::marker(ped, afreq = c(X = 0.5, Y = 0.5), name = "LocusB")
  ped <- pedtools::setMarkers(ped, list(m1, m2))
  ped <- pedtools::setGenotype(ped, id = "1", marker = 1, geno = "A/B")
  ped <- pedtools::setGenotype(ped, id = "1", marker = 2, geno = "X/Y")

  pbn <- from_pedtools(ped)
  expect_true(ncol(pbn$markerEvidence) == 4)
  expect_equal(colnames(pbn$markerEvidence),
               c("LocusA_p", "LocusA_m", "LocusB_p", "LocusB_m"))
})

test_that("from_pedtools markerEvidence contains observed allele labels", {
  skip_if_not_installed("pedtools")
  ped <- pedtools::linearPed(2)
  m1 <- pedtools::marker(ped, afreq = c(A = 0.5, B = 0.5), name = "LocusA")
  ped <- pedtools::setMarkers(ped, m1)
  ped <- pedtools::setGenotype(ped, id = "1", marker = 1, geno = "A/B")

  pbn <- from_pedtools(ped)
  row <- which(rownames(pbn$markerEvidence) == "1")
  expect_true(all(pbn$markerEvidence[row, c("LocusA_p", "LocusA_m")] %in% c("A", "B")))
})
