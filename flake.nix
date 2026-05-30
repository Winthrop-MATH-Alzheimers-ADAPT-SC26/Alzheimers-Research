{
  description = "Basic Julia Flake";
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
      devShells.${system}.default = pkgs.mkShell {
        name = "Julia Flake";

        packages = with pkgs; [
            julia-bin
        ];

        shellHook = ''
            export JULIA_PROJECT=$PWD
        '';
      };
    };
}

