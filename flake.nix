{
  description = "Basic R, Python, and Julia Flake";
  inputs = {
    system-flake.url = "path:/etc/nixos";
    nixpkgs.follows = "system-flake/nixpkgs";
  };

  outputs =
  { nixpkgs, ... }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { 
          inherit system;
          config.allowUnfree = true;
        };
    in
    {
      devShells.${system}.default = pkgs.mkShellNoCC {
        name = "Julia Flake";

        buildInputs = with pkgs; [
          (julia-bin.withPackages [
            "DifferentialEquations"
            "ModelingToolkit"
            "Plots"
          ])
        ];

        shellHook = ''
        '';
      };
    };
}

