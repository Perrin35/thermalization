OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91312042) q[0];
sx q[0];
rz(-1.0364391) q[0];
sx q[0];
rz(1.6057462) q[0];
rz(3.8029382) q[1];
sx q[1];
rz(4.3284419) q[1];
sx q[1];
rz(8.9882589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.614994) q[0];
sx q[0];
rz(-2.8078571) q[0];
sx q[0];
rz(-1.8769391) q[0];
rz(-pi) q[1];
rz(-0.5646602) q[2];
sx q[2];
rz(-1.9604948) q[2];
sx q[2];
rz(-1.0124568) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6923209) q[1];
sx q[1];
rz(-2.2733449) q[1];
sx q[1];
rz(-0.0535123) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7868452) q[3];
sx q[3];
rz(-2.4461544) q[3];
sx q[3];
rz(-1.9239294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46301699) q[2];
sx q[2];
rz(-2.7028694) q[2];
sx q[2];
rz(-2.2185745) q[2];
rz(3.0758514) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(-1.8008697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1317257) q[0];
sx q[0];
rz(-1.0634402) q[0];
sx q[0];
rz(-3.104082) q[0];
rz(-2.5610979) q[1];
sx q[1];
rz(-2.4721804) q[1];
sx q[1];
rz(-2.4494749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0789675) q[0];
sx q[0];
rz(-2.7699372) q[0];
sx q[0];
rz(1.070648) q[0];
rz(-pi) q[1];
rz(2.0742998) q[2];
sx q[2];
rz(-1.8177336) q[2];
sx q[2];
rz(1.4114789) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.925732) q[1];
sx q[1];
rz(-2.070148) q[1];
sx q[1];
rz(0.66470048) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3896319) q[3];
sx q[3];
rz(-1.6979361) q[3];
sx q[3];
rz(2.3746735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6633501) q[2];
sx q[2];
rz(-1.1967836) q[2];
sx q[2];
rz(-1.7459858) q[2];
rz(-0.1869525) q[3];
sx q[3];
rz(-1.170265) q[3];
sx q[3];
rz(2.601534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66889399) q[0];
sx q[0];
rz(-0.63303328) q[0];
sx q[0];
rz(2.582666) q[0];
rz(0.82866296) q[1];
sx q[1];
rz(-1.5907954) q[1];
sx q[1];
rz(2.8899946) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4959481) q[0];
sx q[0];
rz(-2.6584266) q[0];
sx q[0];
rz(-0.42074411) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1220268) q[2];
sx q[2];
rz(-1.2740268) q[2];
sx q[2];
rz(-3.0574329) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4797552) q[1];
sx q[1];
rz(-0.91714707) q[1];
sx q[1];
rz(2.8089671) q[1];
rz(2.2280327) q[3];
sx q[3];
rz(-1.0042448) q[3];
sx q[3];
rz(-2.4853137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.33637968) q[2];
sx q[2];
rz(-2.7153335) q[2];
sx q[2];
rz(-1.4317929) q[2];
rz(-2.924887) q[3];
sx q[3];
rz(-1.6101937) q[3];
sx q[3];
rz(-1.3592892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25877229) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(1.850542) q[0];
rz(2.3123815) q[1];
sx q[1];
rz(-1.1081089) q[1];
sx q[1];
rz(2.1270027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9982581) q[0];
sx q[0];
rz(-0.59276544) q[0];
sx q[0];
rz(0.44500764) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9564187) q[2];
sx q[2];
rz(-2.1004538) q[2];
sx q[2];
rz(-2.5593064) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7493128) q[1];
sx q[1];
rz(-1.5754414) q[1];
sx q[1];
rz(-2.1455812) q[1];
x q[2];
rz(-2.3865256) q[3];
sx q[3];
rz(-1.4658548) q[3];
sx q[3];
rz(-1.5615095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0097051) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(0.76060549) q[2];
rz(-0.46918121) q[3];
sx q[3];
rz(-1.4696591) q[3];
sx q[3];
rz(-0.73452264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0483265) q[0];
sx q[0];
rz(-2.6272197) q[0];
sx q[0];
rz(0.58116466) q[0];
rz(-1.5900853) q[1];
sx q[1];
rz(-0.64684144) q[1];
sx q[1];
rz(0.69492984) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0804176) q[0];
sx q[0];
rz(-1.8616849) q[0];
sx q[0];
rz(-0.24334749) q[0];
rz(-1.476711) q[2];
sx q[2];
rz(-1.861899) q[2];
sx q[2];
rz(-1.7895607) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.123053) q[1];
sx q[1];
rz(-0.54529069) q[1];
sx q[1];
rz(-3.0569949) q[1];
rz(-pi) q[2];
rz(-1.7534689) q[3];
sx q[3];
rz(-1.338306) q[3];
sx q[3];
rz(-2.7894265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.27802262) q[2];
sx q[2];
rz(-1.6081955) q[2];
sx q[2];
rz(0.30280534) q[2];
rz(1.4263724) q[3];
sx q[3];
rz(-0.36077603) q[3];
sx q[3];
rz(0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645709) q[0];
sx q[0];
rz(-0.40957054) q[0];
sx q[0];
rz(2.4466178) q[0];
rz(2.8547844) q[1];
sx q[1];
rz(-2.5359055) q[1];
sx q[1];
rz(-1.274775) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0114637) q[0];
sx q[0];
rz(-2.0507567) q[0];
sx q[0];
rz(-2.7897054) q[0];
x q[1];
rz(-1.1026682) q[2];
sx q[2];
rz(-2.063836) q[2];
sx q[2];
rz(-1.4725034) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3466323) q[1];
sx q[1];
rz(-2.2169276) q[1];
sx q[1];
rz(0.45171311) q[1];
x q[2];
rz(-0.35437472) q[3];
sx q[3];
rz(-2.5297926) q[3];
sx q[3];
rz(2.967088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0605165) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(-0.9673574) q[2];
rz(-2.4216912) q[3];
sx q[3];
rz(-2.1981809) q[3];
sx q[3];
rz(1.2758183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28441456) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(-1.4469294) q[0];
rz(-2.6743496) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(1.1955059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0366096) q[0];
sx q[0];
rz(-0.65147841) q[0];
sx q[0];
rz(0.81654136) q[0];
x q[1];
rz(2.823916) q[2];
sx q[2];
rz(-0.84140771) q[2];
sx q[2];
rz(-0.22401545) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0250281) q[1];
sx q[1];
rz(-2.0691074) q[1];
sx q[1];
rz(-0.070498585) q[1];
rz(-1.6067664) q[3];
sx q[3];
rz(-1.5835053) q[3];
sx q[3];
rz(-0.1606687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82470992) q[2];
sx q[2];
rz(-2.3066545) q[2];
sx q[2];
rz(2.6173124) q[2];
rz(0.25238642) q[3];
sx q[3];
rz(-3.0350244) q[3];
sx q[3];
rz(3.0805123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1687932) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(1.2267145) q[0];
rz(0.75617689) q[1];
sx q[1];
rz(-2.1935479) q[1];
sx q[1];
rz(0.39047584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062138012) q[0];
sx q[0];
rz(-1.3083757) q[0];
sx q[0];
rz(-0.45824893) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0608589) q[2];
sx q[2];
rz(-0.44534031) q[2];
sx q[2];
rz(-0.75917086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9041216) q[1];
sx q[1];
rz(-1.3573705) q[1];
sx q[1];
rz(-1.2280812) q[1];
rz(-pi) q[2];
x q[2];
rz(2.195695) q[3];
sx q[3];
rz(-2.3970553) q[3];
sx q[3];
rz(-2.8685399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32435027) q[2];
sx q[2];
rz(-2.4411026) q[2];
sx q[2];
rz(1.3520757) q[2];
rz(2.9663441) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(1.2043183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4125724) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(2.4793258) q[0];
rz(2.5114255) q[1];
sx q[1];
rz(-1.2799193) q[1];
sx q[1];
rz(-2.9060649) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2199687) q[0];
sx q[0];
rz(-1.3597836) q[0];
sx q[0];
rz(-1.7485662) q[0];
rz(-pi) q[1];
rz(-0.14562433) q[2];
sx q[2];
rz(-2.5674985) q[2];
sx q[2];
rz(0.51806322) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.16727) q[1];
sx q[1];
rz(-3.0758) q[1];
sx q[1];
rz(-1.9386824) q[1];
rz(-pi) q[2];
rz(3.0074545) q[3];
sx q[3];
rz(-0.87867498) q[3];
sx q[3];
rz(1.9798756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0887289) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(-2.5843184) q[2];
rz(2.3142464) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(1.6415589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92852229) q[0];
sx q[0];
rz(-2.3445573) q[0];
sx q[0];
rz(2.573977) q[0];
rz(-0.40191832) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(1.3622805) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6797736) q[0];
sx q[0];
rz(-1.1704233) q[0];
sx q[0];
rz(0.43433614) q[0];
x q[1];
rz(0.2593597) q[2];
sx q[2];
rz(-2.1742714) q[2];
sx q[2];
rz(-0.063613907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6580672) q[1];
sx q[1];
rz(-2.1999252) q[1];
sx q[1];
rz(0.097112819) q[1];
x q[2];
rz(0.79093334) q[3];
sx q[3];
rz(-0.91703992) q[3];
sx q[3];
rz(-0.66822754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7207429) q[2];
sx q[2];
rz(-1.9385312) q[2];
sx q[2];
rz(0.26025772) q[2];
rz(-2.4506954) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(-1.6183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33547587) q[0];
sx q[0];
rz(-1.5609043) q[0];
sx q[0];
rz(1.5326395) q[0];
rz(0.43329049) q[1];
sx q[1];
rz(-1.3229803) q[1];
sx q[1];
rz(1.9569474) q[1];
rz(0.70788371) q[2];
sx q[2];
rz(-2.08692) q[2];
sx q[2];
rz(-0.89649123) q[2];
rz(0.79037249) q[3];
sx q[3];
rz(-1.4740667) q[3];
sx q[3];
rz(2.4355751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
