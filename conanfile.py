from conan import ConanFile
from conan.tools.cmake import CMakeDeps, CMakeToolchain


class MetrikoConan(ConanFile):
    settings = "os", "arch", "compiler", "build_type"

    def requirements(self):
        self.requires("nanobind/2.12.0")
        self.requires("cli11/2.6.0")
        if self.settings.os == "Macos":
            self.requires("llvm-openmp/18.1.8")
        else:
            self.requires("openblas/0.3.30")

    def configure(self):
        if self.settings.os != "Macos":
            self.options["openblas"].shared = False
            msvc = self.settings.compiler == "msvc"
            self.options["openblas"].dynamic_arch = not msvc
            if msvc:
                self.options["openblas"].target = "HASWELL"
            self.options["openblas"].build_lapack = True

    def generate(self):
        CMakeToolchain(self).generate()
        CMakeDeps(self).generate()
