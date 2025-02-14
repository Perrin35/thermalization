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
rz(2.2284722) q[0];
sx q[0];
rz(4.1780317) q[0];
sx q[0];
rz(10.960624) q[0];
rz(-2.4802471) q[1];
sx q[1];
rz(-1.9547434) q[1];
sx q[1];
rz(-0.43651906) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.614994) q[0];
sx q[0];
rz(-0.33373555) q[0];
sx q[0];
rz(-1.8769391) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5646602) q[2];
sx q[2];
rz(-1.1810978) q[2];
sx q[2];
rz(-2.1291358) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4492717) q[1];
sx q[1];
rz(-2.2733449) q[1];
sx q[1];
rz(-0.0535123) q[1];
rz(2.2546886) q[3];
sx q[3];
rz(-1.4330079) q[3];
sx q[3];
rz(2.6215214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6785757) q[2];
sx q[2];
rz(-2.7028694) q[2];
sx q[2];
rz(-2.2185745) q[2];
rz(-0.065741278) q[3];
sx q[3];
rz(-1.8140503) q[3];
sx q[3];
rz(1.8008697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098669212) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(-0.03751066) q[0];
rz(0.58049479) q[1];
sx q[1];
rz(-2.4721804) q[1];
sx q[1];
rz(-2.4494749) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0626252) q[0];
sx q[0];
rz(-0.37165549) q[0];
sx q[0];
rz(1.070648) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28022556) q[2];
sx q[2];
rz(-2.0576653) q[2];
sx q[2];
rz(-0.025472783) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.925732) q[1];
sx q[1];
rz(-2.070148) q[1];
sx q[1];
rz(-0.66470048) q[1];
rz(-2.9565937) q[3];
sx q[3];
rz(-0.76055148) q[3];
sx q[3];
rz(2.4724677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6633501) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(1.7459858) q[2];
rz(0.1869525) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(2.601534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66889399) q[0];
sx q[0];
rz(-2.5085594) q[0];
sx q[0];
rz(0.5589267) q[0];
rz(-2.3129297) q[1];
sx q[1];
rz(-1.5507973) q[1];
sx q[1];
rz(-2.8899946) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0280625) q[0];
sx q[0];
rz(-2.0087272) q[0];
sx q[0];
rz(-1.3597041) q[0];
rz(-pi) q[1];
rz(-2.1220268) q[2];
sx q[2];
rz(-1.8675659) q[2];
sx q[2];
rz(3.0574329) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66183749) q[1];
sx q[1];
rz(-2.2244456) q[1];
sx q[1];
rz(2.8089671) q[1];
x q[2];
rz(2.4647242) q[3];
sx q[3];
rz(-2.1123611) q[3];
sx q[3];
rz(-1.8344017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33637968) q[2];
sx q[2];
rz(-2.7153335) q[2];
sx q[2];
rz(-1.7097998) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.5313989) q[3];
sx q[3];
rz(-1.3592892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25877229) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(-1.2910507) q[0];
rz(-2.3123815) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(-1.0145899) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4764164) q[0];
sx q[0];
rz(-1.0422857) q[0];
sx q[0];
rz(1.2885874) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9564187) q[2];
sx q[2];
rz(-1.0411388) q[2];
sx q[2];
rz(-0.5822863) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9559052) q[1];
sx q[1];
rz(-2.5667911) q[1];
sx q[1];
rz(-1.5622524) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3865256) q[3];
sx q[3];
rz(-1.4658548) q[3];
sx q[3];
rz(1.5615095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1318876) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(-2.3809872) q[2];
rz(-2.6724114) q[3];
sx q[3];
rz(-1.4696591) q[3];
sx q[3];
rz(-2.40707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0483265) q[0];
sx q[0];
rz(-0.51437298) q[0];
sx q[0];
rz(2.560428) q[0];
rz(-1.5900853) q[1];
sx q[1];
rz(-0.64684144) q[1];
sx q[1];
rz(-2.4466628) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3669156) q[0];
sx q[0];
rz(-0.37702501) q[0];
sx q[0];
rz(2.24848) q[0];
rz(1.6648817) q[2];
sx q[2];
rz(-1.2796937) q[2];
sx q[2];
rz(1.7895607) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.123053) q[1];
sx q[1];
rz(-2.596302) q[1];
sx q[1];
rz(-3.0569949) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4871553) q[3];
sx q[3];
rz(-2.8469466) q[3];
sx q[3];
rz(-2.1135995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.86357) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(2.8387873) q[2];
rz(1.4263724) q[3];
sx q[3];
rz(-0.36077603) q[3];
sx q[3];
rz(0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645709) q[0];
sx q[0];
rz(-2.7320221) q[0];
sx q[0];
rz(0.69497481) q[0];
rz(0.28680828) q[1];
sx q[1];
rz(-0.60568714) q[1];
sx q[1];
rz(1.8668176) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0114637) q[0];
sx q[0];
rz(-2.0507567) q[0];
sx q[0];
rz(-0.35188727) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4430577) q[2];
sx q[2];
rz(-2.4753127) q[2];
sx q[2];
rz(-2.4874788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50837738) q[1];
sx q[1];
rz(-1.214809) q[1];
sx q[1];
rz(0.87320019) q[1];
x q[2];
rz(-0.35437472) q[3];
sx q[3];
rz(-2.5297926) q[3];
sx q[3];
rz(2.967088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0605165) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(-2.1742353) q[2];
rz(-2.4216912) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(-1.2758183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28441456) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(-1.6946633) q[0];
rz(-0.46724304) q[1];
sx q[1];
rz(-1.2930361) q[1];
sx q[1];
rz(1.1955059) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0366096) q[0];
sx q[0];
rz(-2.4901142) q[0];
sx q[0];
rz(0.81654136) q[0];
rz(2.823916) q[2];
sx q[2];
rz(-0.84140771) q[2];
sx q[2];
rz(-0.22401545) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0250281) q[1];
sx q[1];
rz(-2.0691074) q[1];
sx q[1];
rz(-3.0710941) q[1];
rz(1.5348263) q[3];
sx q[3];
rz(-1.5580873) q[3];
sx q[3];
rz(-2.980924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82470992) q[2];
sx q[2];
rz(-2.3066545) q[2];
sx q[2];
rz(-2.6173124) q[2];
rz(-0.25238642) q[3];
sx q[3];
rz(-0.10656825) q[3];
sx q[3];
rz(-0.061080385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97279945) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(-1.9148781) q[0];
rz(-0.75617689) q[1];
sx q[1];
rz(-0.94804472) q[1];
sx q[1];
rz(0.39047584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3813864) q[0];
sx q[0];
rz(-1.1293656) q[0];
sx q[0];
rz(-1.2797829) q[0];
rz(-pi) q[1];
rz(3.0608589) q[2];
sx q[2];
rz(-2.6962523) q[2];
sx q[2];
rz(-0.75917086) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8836964) q[1];
sx q[1];
rz(-1.2361649) q[1];
sx q[1];
rz(0.22617126) q[1];
rz(2.195695) q[3];
sx q[3];
rz(-0.74453738) q[3];
sx q[3];
rz(-0.27305279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8172424) q[2];
sx q[2];
rz(-0.70049006) q[2];
sx q[2];
rz(-1.3520757) q[2];
rz(-2.9663441) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(1.9372743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7290203) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(0.66226688) q[0];
rz(0.63016713) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(0.23552775) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5109638) q[0];
sx q[0];
rz(-0.27505829) q[0];
sx q[0];
rz(-0.69013005) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.477219) q[2];
sx q[2];
rz(-1.0035328) q[2];
sx q[2];
rz(2.7964489) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79865658) q[1];
sx q[1];
rz(-1.5094116) q[1];
sx q[1];
rz(3.117901) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13413818) q[3];
sx q[3];
rz(-2.2629177) q[3];
sx q[3];
rz(1.1617171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0528637) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(2.5843184) q[2];
rz(-2.3142464) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(-1.6415589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2130704) q[0];
sx q[0];
rz(-2.3445573) q[0];
sx q[0];
rz(-0.56761566) q[0];
rz(0.40191832) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(1.7793122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80776725) q[0];
sx q[0];
rz(-2.5596715) q[0];
sx q[0];
rz(0.78827095) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1902607) q[2];
sx q[2];
rz(-1.3580322) q[2];
sx q[2];
rz(-1.484953) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6475923) q[1];
sx q[1];
rz(-0.63557445) q[1];
sx q[1];
rz(1.4383491) q[1];
rz(-2.3990223) q[3];
sx q[3];
rz(-0.97108632) q[3];
sx q[3];
rz(0.35123435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7207429) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(-2.8813349) q[2];
rz(0.69089729) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(-1.6183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(0.92930195) q[2];
sx q[2];
rz(-2.1718696) q[2];
sx q[2];
rz(1.0739506) q[2];
rz(2.3512202) q[3];
sx q[3];
rz(-1.6675259) q[3];
sx q[3];
rz(-0.70601757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
