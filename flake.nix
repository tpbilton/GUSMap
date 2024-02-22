{
  description = "GUSMap";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
    GUSbase = {
      url = github:tpbilton/GUSbase;
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
    };
  };

  outputs = { self, nixpkgs, flake-utils, GUSbase }:
    flake-utils.lib.eachDefaultSystem
      (system:
        let
          pkgs = import nixpkgs {
            inherit system;
          };

          flakePkgs = {
            GUSbase = GUSbase.packages.${system}.default;
          };

          GUSMap = with pkgs;
            rPackages.buildRPackage {
              name = "GUSMap";
              src = ./.;
              propagatedBuildInputs = with rPackages;
                [ R6 smacof princurve Rdpack data_table parallel doParallel iterators foreach Rcpp plotly flakePkgs.GUSbase ];
            };

          R-with-GUSMap = with pkgs;
            rWrapper.override {
              packages = with rPackages;
                [ GUSMap ];
            };
        in
          with pkgs;
          {
            devShells.default = mkShell {
              buildInputs = [ R-with-GUSMap ];
              shellHook = ''
            mkdir -p "$(pwd)/_libs"
            export R_LIBS_USER="$(pwd)/_libs"
          '';
            };

            packages.default = GUSMap;
          });
}
