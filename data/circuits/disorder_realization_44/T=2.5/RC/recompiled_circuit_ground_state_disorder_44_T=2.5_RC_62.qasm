OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39785102) q[0];
sx q[0];
rz(-2.1704817) q[0];
sx q[0];
rz(0.71075332) q[0];
rz(0.96762586) q[1];
sx q[1];
rz(-3.1275446) q[1];
sx q[1];
rz(0.79298055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4913038) q[0];
sx q[0];
rz(-2.2772335) q[0];
sx q[0];
rz(2.6992348) q[0];
rz(-pi) q[1];
x q[1];
rz(0.014292467) q[2];
sx q[2];
rz(-1.5231175) q[2];
sx q[2];
rz(-2.45243) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88411623) q[1];
sx q[1];
rz(-2.3168644) q[1];
sx q[1];
rz(0.60666879) q[1];
x q[2];
rz(0.96785424) q[3];
sx q[3];
rz(-2.0857028) q[3];
sx q[3];
rz(0.076857278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8636318) q[2];
sx q[2];
rz(-2.7148254) q[2];
sx q[2];
rz(-2.5521736) q[2];
rz(-0.74728084) q[3];
sx q[3];
rz(-0.092844754) q[3];
sx q[3];
rz(0.59982991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079161949) q[0];
sx q[0];
rz(-2.9099162) q[0];
sx q[0];
rz(2.2404501) q[0];
rz(-0.98183739) q[1];
sx q[1];
rz(-0.58126175) q[1];
sx q[1];
rz(-2.1493886) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81047786) q[0];
sx q[0];
rz(-1.4338636) q[0];
sx q[0];
rz(2.1934319) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7478701) q[2];
sx q[2];
rz(-2.3858527) q[2];
sx q[2];
rz(2.9301639) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7957669) q[1];
sx q[1];
rz(-2.2123033) q[1];
sx q[1];
rz(-0.87916763) q[1];
x q[2];
rz(-0.91199947) q[3];
sx q[3];
rz(-2.1423295) q[3];
sx q[3];
rz(2.6603572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.046752669) q[2];
sx q[2];
rz(-2.6111626) q[2];
sx q[2];
rz(0.025010427) q[2];
rz(-2.181634) q[3];
sx q[3];
rz(-3.0955866) q[3];
sx q[3];
rz(0.068232603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18441021) q[0];
sx q[0];
rz(-0.010951696) q[0];
sx q[0];
rz(-2.4175194) q[0];
rz(-3.030576) q[1];
sx q[1];
rz(-2.5357775) q[1];
sx q[1];
rz(0.012880005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52791053) q[0];
sx q[0];
rz(-0.22697313) q[0];
sx q[0];
rz(1.2140973) q[0];
rz(-2.3518042) q[2];
sx q[2];
rz(-0.26726535) q[2];
sx q[2];
rz(-2.546026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15808039) q[1];
sx q[1];
rz(-1.3338102) q[1];
sx q[1];
rz(1.3999585) q[1];
rz(0.37498388) q[3];
sx q[3];
rz(-0.91911784) q[3];
sx q[3];
rz(0.76256547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24590242) q[2];
sx q[2];
rz(-0.7220214) q[2];
sx q[2];
rz(-0.32430172) q[2];
rz(-2.6445828) q[3];
sx q[3];
rz(-0.23247601) q[3];
sx q[3];
rz(2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67320353) q[0];
sx q[0];
rz(-0.16981801) q[0];
sx q[0];
rz(-2.66535) q[0];
rz(2.6561123) q[1];
sx q[1];
rz(-2.6464033) q[1];
sx q[1];
rz(2.9229497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63594288) q[0];
sx q[0];
rz(-2.6723862) q[0];
sx q[0];
rz(-2.0903265) q[0];
x q[1];
rz(-2.4179732) q[2];
sx q[2];
rz(-1.5800522) q[2];
sx q[2];
rz(1.1030359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87590295) q[1];
sx q[1];
rz(-0.69702083) q[1];
sx q[1];
rz(-0.69209309) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5095482) q[3];
sx q[3];
rz(-2.7586652) q[3];
sx q[3];
rz(2.6296089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.789088) q[2];
sx q[2];
rz(-2.3829491) q[2];
sx q[2];
rz(-0.57141203) q[2];
rz(2.2032951) q[3];
sx q[3];
rz(-0.58656991) q[3];
sx q[3];
rz(-2.9686484) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5649696) q[0];
sx q[0];
rz(-2.7374856) q[0];
sx q[0];
rz(0.47082666) q[0];
rz(-2.2305523) q[1];
sx q[1];
rz(-0.43993479) q[1];
sx q[1];
rz(2.7104673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.39931) q[0];
sx q[0];
rz(-0.78680187) q[0];
sx q[0];
rz(1.6086701) q[0];
x q[1];
rz(1.9091102) q[2];
sx q[2];
rz(-2.4801284) q[2];
sx q[2];
rz(-2.4544883) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2700165) q[1];
sx q[1];
rz(-1.5432913) q[1];
sx q[1];
rz(1.5020976) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60089941) q[3];
sx q[3];
rz(-2.89866) q[3];
sx q[3];
rz(3.0593135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8016781) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(-0.3413631) q[2];
rz(-0.3723799) q[3];
sx q[3];
rz(-2.5058993) q[3];
sx q[3];
rz(2.671833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895198) q[0];
sx q[0];
rz(-2.2884123) q[0];
sx q[0];
rz(3.0186655) q[0];
rz(2.7561103) q[1];
sx q[1];
rz(-0.79817927) q[1];
sx q[1];
rz(0.72365671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48945828) q[0];
sx q[0];
rz(-1.3969027) q[0];
sx q[0];
rz(3.0956718) q[0];
rz(-pi) q[1];
rz(-2.0711818) q[2];
sx q[2];
rz(-0.37795174) q[2];
sx q[2];
rz(-2.9674781) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8116709) q[1];
sx q[1];
rz(-0.99024665) q[1];
sx q[1];
rz(0.73402053) q[1];
rz(-0.21435634) q[3];
sx q[3];
rz(-1.949614) q[3];
sx q[3];
rz(-2.5942868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.75189292) q[2];
sx q[2];
rz(-0.60201001) q[2];
sx q[2];
rz(2.4683118) q[2];
rz(0.26345396) q[3];
sx q[3];
rz(-2.6665688) q[3];
sx q[3];
rz(0.30509216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4398956) q[0];
sx q[0];
rz(-2.3450527) q[0];
sx q[0];
rz(2.6252966) q[0];
rz(2.3856178) q[1];
sx q[1];
rz(-0.95840234) q[1];
sx q[1];
rz(0.36852536) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4004423) q[0];
sx q[0];
rz(-1.0545122) q[0];
sx q[0];
rz(2.3325066) q[0];
rz(-pi) q[1];
rz(3.0958648) q[2];
sx q[2];
rz(-0.70529443) q[2];
sx q[2];
rz(-0.7208342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75223756) q[1];
sx q[1];
rz(-0.19986831) q[1];
sx q[1];
rz(-1.0081069) q[1];
rz(-pi) q[2];
rz(0.51492274) q[3];
sx q[3];
rz(-0.92149261) q[3];
sx q[3];
rz(-2.2442371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4964909) q[2];
sx q[2];
rz(-0.056702159) q[2];
sx q[2];
rz(0.48218316) q[2];
rz(-2.9218946) q[3];
sx q[3];
rz(-2.4441661) q[3];
sx q[3];
rz(-2.4612332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41736233) q[0];
sx q[0];
rz(-0.14886947) q[0];
sx q[0];
rz(-2.505488) q[0];
rz(-2.1813189) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(2.2669534) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68976328) q[0];
sx q[0];
rz(-1.4216107) q[0];
sx q[0];
rz(-0.058663603) q[0];
rz(0.12107559) q[2];
sx q[2];
rz(-1.6796475) q[2];
sx q[2];
rz(-1.4360652) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0046152701) q[1];
sx q[1];
rz(-0.57356131) q[1];
sx q[1];
rz(0.80570813) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88040027) q[3];
sx q[3];
rz(-1.0171618) q[3];
sx q[3];
rz(0.88144377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8735698) q[2];
sx q[2];
rz(-0.41646725) q[2];
sx q[2];
rz(-2.2250309) q[2];
rz(-0.77740866) q[3];
sx q[3];
rz(-2.58367) q[3];
sx q[3];
rz(-0.53434813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1189608) q[0];
sx q[0];
rz(-0.73771483) q[0];
sx q[0];
rz(-2.90888) q[0];
rz(-2.4932056) q[1];
sx q[1];
rz(-2.5189923) q[1];
sx q[1];
rz(-2.5737305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083160087) q[0];
sx q[0];
rz(-0.55480236) q[0];
sx q[0];
rz(-2.2834999) q[0];
rz(-pi) q[1];
rz(-1.2987192) q[2];
sx q[2];
rz(-1.7781742) q[2];
sx q[2];
rz(0.033471154) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76946867) q[1];
sx q[1];
rz(-2.8903277) q[1];
sx q[1];
rz(2.1970046) q[1];
rz(0.81553163) q[3];
sx q[3];
rz(-2.14912) q[3];
sx q[3];
rz(2.2165143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76967543) q[2];
sx q[2];
rz(-0.19352517) q[2];
sx q[2];
rz(0.5522716) q[2];
rz(2.8596089) q[3];
sx q[3];
rz(-0.43006399) q[3];
sx q[3];
rz(-3.0388487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.929739) q[0];
sx q[0];
rz(-0.064082853) q[0];
sx q[0];
rz(0.1317568) q[0];
rz(-0.53648221) q[1];
sx q[1];
rz(-2.9832612) q[1];
sx q[1];
rz(-0.90824711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2836766) q[0];
sx q[0];
rz(-0.87661298) q[0];
sx q[0];
rz(-0.038900872) q[0];
x q[1];
rz(-1.0824049) q[2];
sx q[2];
rz(-2.6219212) q[2];
sx q[2];
rz(-2.1888417) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6301292) q[1];
sx q[1];
rz(-2.2591233) q[1];
sx q[1];
rz(-0.87911112) q[1];
rz(-pi) q[2];
rz(1.578701) q[3];
sx q[3];
rz(-2.9156682) q[3];
sx q[3];
rz(0.28858063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6324255) q[2];
sx q[2];
rz(-2.2278892) q[2];
sx q[2];
rz(0.2779648) q[2];
rz(0.5667423) q[3];
sx q[3];
rz(-2.9394579) q[3];
sx q[3];
rz(-2.22866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81000281) q[0];
sx q[0];
rz(-1.4053874) q[0];
sx q[0];
rz(-0.94595861) q[0];
rz(-0.48238659) q[1];
sx q[1];
rz(-1.9079897) q[1];
sx q[1];
rz(2.447396) q[1];
rz(1.9271156) q[2];
sx q[2];
rz(-0.6206442) q[2];
sx q[2];
rz(0.71833687) q[2];
rz(1.0645772) q[3];
sx q[3];
rz(-1.5983221) q[3];
sx q[3];
rz(1.7479001) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
