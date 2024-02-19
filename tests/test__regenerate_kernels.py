import xtrack as xt

def test_init():
    xt.prebuild_kernels.regenerate_kernels(kernels=[
        "default_xcoll_only_absorbers", "default_xcoll",
        "default_xcoll_crystals"])
