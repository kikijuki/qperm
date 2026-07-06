# Regenerate data/qae_demo.rda
#
# qae_demo is the small synthetic one-sample contrast dataset used by the interactive
# demo in the original project (10 subjects, 6 voxels). It is small enough that the
# entire sign-flip null (2^10 = 1024 relabelings) can be enumerated exactly, which is
# what lets the quantum estimate be checked against exact ground truth.

qae_demo <- matrix(c(
   1.0512,  0.2987, -0.2741, -0.8906, -0.4547, -0.9916,
   1.1101,  1.3402, -0.4922, -0.6205,  0.4898,  0.3569,
   1.1554, -0.9305, -0.0293,  0.6953, -1.3442, -0.4576,
  -0.8512, -1.2895, -1.8417, -0.2351, -1.2674,  0.2713,
   1.2068, -0.1869, -2.5168, -0.5387, -0.0485,  0.1133,
  -0.4801, -0.4778, -0.9785, -0.8088,  1.0609, -0.8075,
   1.0175,  0.8844, -0.5836, -0.1117,  0.1105,  0.0638,
  -0.1751,  0.0761,  1.3588, -1.5471,  0.8594,  0.1194,
   0.4085,  2.0004,  0.7623, -1.1993,  0.0745,  0.5767,
   0.8612,  0.6829, -0.0665,  0.6672,  1.4385, -0.6757),
  nrow = 10, ncol = 6, byrow = TRUE,
  dimnames = list(paste0("subject", 1:10), paste0("voxel", 1:6)))

save(qae_demo, file = "data/qae_demo.rda", version = 2)
cat("wrote data/qae_demo.rda\n")
