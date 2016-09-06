# author: Bo Han (bhan@pacb.com)

declare -xr COLOR_SET=1

declare -xr FONT_COLOR_RED="\e[31;40m";
declare -xr FONT_COLOR_GREEN="\e[32;40m";
declare -xr FONT_COLOR_YELLOW="\e[33;40m";
declare -xr FONT_COLOR_BLUE="\e[34;40m";
declare -xr FONT_COLOR_MAGENTA="\e[35;40m";
declare -xr FONT_COLOR_CYAN="\e[36;40m";
declare -xr FONT_COLOR_RESET="\e[0m";

declare -xr FONT_STYLE_BOLD='\e[1m';
declare -xr FONT_STYLE_UNDERLINE='\e[4m';
declare -xr FONT_STYLE_BLINK='\e[5m';
declare -xr FONT_STYLE_RESET='\e[0m';

declare -xr COMMENT=`echo -ne $FONT_COLOR_BLUE`;
declare -xr REQUIRED=`echo -ne $FONT_COLOR_RED`;
declare -xr OPTIONAL=`echo -ne $FONT_COLOR_GREEN`;
declare -xr ADVANCED=`echo -ne $FONT_COLOR_MAGENTA`;
declare -xr UNDERDEV=`echo -ne $FONT_COLOR_CYAN`;
