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

      julia-env = pkgs.buildFHSEnv {
        name = "julia-env";
        targetPkgs = pkgs: with pkgs; [
          julia-bin
          # GR/gksqt dependencies
          libGL
          libGLU
          libxt
          libx11
          libxrender
          libxext
          glfw
          freetype
          qt5.qtbase
          qt5.qtwayland
        ];
        runScript = "julia";
      };

    in
    {
      devShells.${system}.default = pkgs.mkShell {
        name = "Julia Flake";

        packages = [
          julia-env
        ];

        shellHook = ''
          export JULIA_PROJECT=$PWD
          export QT_PLUGIN_PATH="${pkgs.qt5.qtbase}/${pkgs.qt5.qtbase.qtPluginPrefix}"
        '';
      };
    };
}
