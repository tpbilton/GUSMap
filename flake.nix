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

          # We call this myAGHmatrix, because AGHmatrix is already packaged in Nix,
          # and we're only doing it like this for the purpose of having an example.
          #
          # The rev is obtained by copying the full SHA for the required commit on GitHub.
          # The sha256 is obtained from that rev like this:
          # nix-prefetch-url --unpack https://github.com/prmunoz/AGHmatrix/archive/73e7b743cd57b81957a1aec20126955717507241.tar.gz
          #
          # The propagated build inputs are copied from the Imports sections of the R package DESCRIPTION.
          # If any of these are not packaged in Nix already, they must be made available like we are doing for myAGHmatrix.  Recursively.
          myAGHmatrix = with pkgs;
            rPackages.buildRPackage {
              name = "AGHmatrix";
              src = fetchFromGitHub {
                owner = "prmunoz";
                repo = "AGHmatrix";
                rev = "73e7b743cd57b81957a1aec20126955717507241";
                sha256 = "0whfkflxjxp3qgf4332scz4b3yrvb7v3j0x9iy3vhbl4fxmdapdj";
              };
              propagatedBuildInputs = with rPackages;
                [ Matrix zoo ];
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
                [ GUSMap myAGHmatrix ];
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
