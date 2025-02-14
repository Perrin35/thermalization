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
rz(1.7914766) q[0];
sx q[0];
rz(-0.67576367) q[0];
sx q[0];
rz(-0.0078553353) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(-1.5239198) q[1];
sx q[1];
rz(2.89892) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9822916) q[0];
sx q[0];
rz(-0.86573273) q[0];
sx q[0];
rz(-0.16601913) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6818993) q[2];
sx q[2];
rz(-1.9003344) q[2];
sx q[2];
rz(2.7609563) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85181273) q[1];
sx q[1];
rz(-0.35408005) q[1];
sx q[1];
rz(0.068239958) q[1];
rz(1.3245565) q[3];
sx q[3];
rz(-0.9199577) q[3];
sx q[3];
rz(3.0512187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0365888) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(-0.71199065) q[2];
rz(2.8365734) q[3];
sx q[3];
rz(-2.9192393) q[3];
sx q[3];
rz(0.045507889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63118339) q[0];
sx q[0];
rz(-1.2774066) q[0];
sx q[0];
rz(1.7741868) q[0];
rz(3.101688) q[1];
sx q[1];
rz(-2.3343562) q[1];
sx q[1];
rz(2.5423999) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78135787) q[0];
sx q[0];
rz(-1.0971945) q[0];
sx q[0];
rz(-1.4417157) q[0];
rz(-2.5810601) q[2];
sx q[2];
rz(-2.2403702) q[2];
sx q[2];
rz(-2.5472484) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6015905) q[1];
sx q[1];
rz(-0.7379325) q[1];
sx q[1];
rz(-2.6327391) q[1];
x q[2];
rz(-2.8017524) q[3];
sx q[3];
rz(-2.426894) q[3];
sx q[3];
rz(-3.0106737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0590608) q[2];
sx q[2];
rz(-0.56190562) q[2];
sx q[2];
rz(1.833029) q[2];
rz(-0.023905309) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(1.5786952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6858653) q[0];
sx q[0];
rz(-2.4754334) q[0];
sx q[0];
rz(2.300793) q[0];
rz(-1.0288382) q[1];
sx q[1];
rz(-2.7327171) q[1];
sx q[1];
rz(-1.6811446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90636364) q[0];
sx q[0];
rz(-2.0301986) q[0];
sx q[0];
rz(-1.4109341) q[0];
rz(-pi) q[1];
rz(2.6681294) q[2];
sx q[2];
rz(-1.9046648) q[2];
sx q[2];
rz(-2.3893181) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7531794) q[1];
sx q[1];
rz(-2.3824661) q[1];
sx q[1];
rz(-2.6469346) q[1];
rz(3.06006) q[3];
sx q[3];
rz(-0.60890261) q[3];
sx q[3];
rz(1.0074774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.2568405) q[2];
sx q[2];
rz(-1.2105056) q[2];
sx q[2];
rz(-2.9849226) q[2];
rz(-1.6414075) q[3];
sx q[3];
rz(-1.2431966) q[3];
sx q[3];
rz(-2.959804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0750065) q[0];
sx q[0];
rz(-1.8445419) q[0];
sx q[0];
rz(-0.080667607) q[0];
rz(2.0196041) q[1];
sx q[1];
rz(-0.64067084) q[1];
sx q[1];
rz(-1.5544308) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98469668) q[0];
sx q[0];
rz(-1.5050355) q[0];
sx q[0];
rz(0.82729152) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8055259) q[2];
sx q[2];
rz(-3.0185351) q[2];
sx q[2];
rz(2.6390136) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39496468) q[1];
sx q[1];
rz(-2.0680799) q[1];
sx q[1];
rz(-2.4056466) q[1];
rz(-pi) q[2];
rz(0.46458475) q[3];
sx q[3];
rz(-2.4813093) q[3];
sx q[3];
rz(-2.6990776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2598205) q[2];
sx q[2];
rz(-2.0505003) q[2];
sx q[2];
rz(-2.1176977) q[2];
rz(2.3648868) q[3];
sx q[3];
rz(-1.1506162) q[3];
sx q[3];
rz(-1.0030494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7316932) q[0];
sx q[0];
rz(-2.2875146) q[0];
sx q[0];
rz(-2.5140629) q[0];
rz(-0.69951406) q[1];
sx q[1];
rz(-2.1414089) q[1];
sx q[1];
rz(2.2342009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8555657) q[0];
sx q[0];
rz(-0.81840179) q[0];
sx q[0];
rz(1.6256385) q[0];
rz(0.19030119) q[2];
sx q[2];
rz(-1.7277328) q[2];
sx q[2];
rz(0.67853329) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.34037922) q[1];
sx q[1];
rz(-1.0579018) q[1];
sx q[1];
rz(2.7199634) q[1];
x q[2];
rz(0.058919546) q[3];
sx q[3];
rz(-1.0353966) q[3];
sx q[3];
rz(-0.31188341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0044535) q[2];
sx q[2];
rz(-1.8567825) q[2];
sx q[2];
rz(-2.2966906) q[2];
rz(0.6984624) q[3];
sx q[3];
rz(-2.4013077) q[3];
sx q[3];
rz(-2.3499878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7591105) q[0];
sx q[0];
rz(-1.6805205) q[0];
sx q[0];
rz(1.2680898) q[0];
rz(-1.6146487) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(-0.39438927) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52452625) q[0];
sx q[0];
rz(-1.6163278) q[0];
sx q[0];
rz(0.92638735) q[0];
x q[1];
rz(-3.0494681) q[2];
sx q[2];
rz(-0.73706223) q[2];
sx q[2];
rz(-1.2110405) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.907669) q[1];
sx q[1];
rz(-0.65429293) q[1];
sx q[1];
rz(0.073542417) q[1];
rz(-pi) q[2];
rz(-2.7809393) q[3];
sx q[3];
rz(-1.4980199) q[3];
sx q[3];
rz(0.78950715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4322728) q[2];
sx q[2];
rz(-3.0797112) q[2];
sx q[2];
rz(2.5800932) q[2];
rz(1.1310486) q[3];
sx q[3];
rz(-0.80315042) q[3];
sx q[3];
rz(2.6530182) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5442218) q[0];
sx q[0];
rz(-2.9591296) q[0];
sx q[0];
rz(0.23319787) q[0];
rz(-0.11416642) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(2.9679969) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6703807) q[0];
sx q[0];
rz(-2.2224373) q[0];
sx q[0];
rz(2.6782413) q[0];
rz(-pi) q[1];
rz(0.87145135) q[2];
sx q[2];
rz(-0.55784278) q[2];
sx q[2];
rz(-0.43684549) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44960231) q[1];
sx q[1];
rz(-0.27562818) q[1];
sx q[1];
rz(-0.93494934) q[1];
rz(-pi) q[2];
rz(1.5830481) q[3];
sx q[3];
rz(-1.3821332) q[3];
sx q[3];
rz(1.0808627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0357828) q[2];
sx q[2];
rz(-0.95869392) q[2];
sx q[2];
rz(-2.2868273) q[2];
rz(0.0828951) q[3];
sx q[3];
rz(-1.1707183) q[3];
sx q[3];
rz(-0.054072592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6524803) q[0];
sx q[0];
rz(-1.880045) q[0];
sx q[0];
rz(-1.1789119) q[0];
rz(1.0014125) q[1];
sx q[1];
rz(-1.2636355) q[1];
sx q[1];
rz(2.9875535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6545367) q[0];
sx q[0];
rz(-2.3628919) q[0];
sx q[0];
rz(0.47327431) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70239046) q[2];
sx q[2];
rz(-0.99961126) q[2];
sx q[2];
rz(-2.1399501) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5370085) q[1];
sx q[1];
rz(-1.2717383) q[1];
sx q[1];
rz(-1.0826151) q[1];
rz(-pi) q[2];
rz(-1.834112) q[3];
sx q[3];
rz(-1.8510518) q[3];
sx q[3];
rz(1.6681886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68086326) q[2];
sx q[2];
rz(-2.0622084) q[2];
sx q[2];
rz(-2.4737849) q[2];
rz(-1.7364511) q[3];
sx q[3];
rz(-1.3677771) q[3];
sx q[3];
rz(0.63647979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77968303) q[0];
sx q[0];
rz(-1.6661665) q[0];
sx q[0];
rz(2.5478126) q[0];
rz(1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(-1.8437754) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5820044) q[0];
sx q[0];
rz(-2.1540717) q[0];
sx q[0];
rz(-1.6571664) q[0];
x q[1];
rz(-0.4164575) q[2];
sx q[2];
rz(-0.67085941) q[2];
sx q[2];
rz(2.351298) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5173856) q[1];
sx q[1];
rz(-0.38100699) q[1];
sx q[1];
rz(2.7966649) q[1];
rz(-1.3993949) q[3];
sx q[3];
rz(-2.3706782) q[3];
sx q[3];
rz(-0.41219974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7456776) q[2];
sx q[2];
rz(-0.021947689) q[2];
sx q[2];
rz(1.6321261) q[2];
rz(-0.11219003) q[3];
sx q[3];
rz(-1.0271007) q[3];
sx q[3];
rz(1.7621015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6701732) q[0];
sx q[0];
rz(-1.0059953) q[0];
sx q[0];
rz(0.26563409) q[0];
rz(-2.1521125) q[1];
sx q[1];
rz(-1.8786636) q[1];
sx q[1];
rz(0.33214733) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50221793) q[0];
sx q[0];
rz(-1.2373588) q[0];
sx q[0];
rz(0.0040472814) q[0];
x q[1];
rz(2.0004326) q[2];
sx q[2];
rz(-0.62295914) q[2];
sx q[2];
rz(-1.0819544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0451584) q[1];
sx q[1];
rz(-1.9109042) q[1];
sx q[1];
rz(-0.051326871) q[1];
rz(-pi) q[2];
rz(1.321855) q[3];
sx q[3];
rz(-2.948108) q[3];
sx q[3];
rz(1.8850808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0997448) q[2];
sx q[2];
rz(-2.0691278) q[2];
sx q[2];
rz(-2.9912046) q[2];
rz(-0.8485052) q[3];
sx q[3];
rz(-2.7917807) q[3];
sx q[3];
rz(1.295804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083241845) q[0];
sx q[0];
rz(-1.3180838) q[0];
sx q[0];
rz(-1.1710118) q[0];
rz(-1.3311483) q[1];
sx q[1];
rz(-2.34453) q[1];
sx q[1];
rz(-1.1042368) q[1];
rz(-1.4990357) q[2];
sx q[2];
rz(-2.4007779) q[2];
sx q[2];
rz(3.130198) q[2];
rz(1.5965309) q[3];
sx q[3];
rz(-2.2774057) q[3];
sx q[3];
rz(-2.0792815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
