OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.42216766) q[0];
sx q[0];
rz(7.0081975) q[0];
sx q[0];
rz(10.260697) q[0];
rz(-5.4556203) q[1];
sx q[1];
rz(6.5536708) q[1];
sx q[1];
rz(9.8319893) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70449984) q[0];
sx q[0];
rz(-1.8665534) q[0];
sx q[0];
rz(-2.9192135) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9303481) q[2];
sx q[2];
rz(-1.7952732) q[2];
sx q[2];
rz(-1.4663855) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5558124) q[1];
sx q[1];
rz(-1.4364873) q[1];
sx q[1];
rz(2.7463169) q[1];
rz(-pi) q[2];
rz(-0.9933957) q[3];
sx q[3];
rz(-2.2334986) q[3];
sx q[3];
rz(0.71383324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.779458) q[2];
sx q[2];
rz(-2.802765) q[2];
sx q[2];
rz(-0.96620488) q[2];
rz(0.90159121) q[3];
sx q[3];
rz(-2.2880771) q[3];
sx q[3];
rz(-2.0961659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40039429) q[0];
sx q[0];
rz(-0.61120954) q[0];
sx q[0];
rz(0.10261593) q[0];
rz(-1.9727033) q[1];
sx q[1];
rz(-2.1714307) q[1];
sx q[1];
rz(0.51365596) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1875604) q[0];
sx q[0];
rz(-2.6488049) q[0];
sx q[0];
rz(1.378118) q[0];
x q[1];
rz(-2.9048237) q[2];
sx q[2];
rz(-2.5920923) q[2];
sx q[2];
rz(2.7236746) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8508965) q[1];
sx q[1];
rz(-0.65819955) q[1];
sx q[1];
rz(1.0653516) q[1];
rz(2.0046141) q[3];
sx q[3];
rz(-1.8469212) q[3];
sx q[3];
rz(-1.9406589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1932842) q[2];
sx q[2];
rz(-2.9479492) q[2];
sx q[2];
rz(-0.2085169) q[2];
rz(-0.68566132) q[3];
sx q[3];
rz(-1.0582358) q[3];
sx q[3];
rz(2.9764777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0971646) q[0];
sx q[0];
rz(-1.8655638) q[0];
sx q[0];
rz(-1.2318508) q[0];
rz(2.7992898) q[1];
sx q[1];
rz(-2.0399703) q[1];
sx q[1];
rz(-1.5922348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78944713) q[0];
sx q[0];
rz(-1.88228) q[0];
sx q[0];
rz(1.3793263) q[0];
rz(-pi) q[1];
rz(-0.64022417) q[2];
sx q[2];
rz(-2.9415543) q[2];
sx q[2];
rz(2.6286516) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.07330896) q[1];
sx q[1];
rz(-1.4335634) q[1];
sx q[1];
rz(2.2925966) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54247673) q[3];
sx q[3];
rz(-2.0869792) q[3];
sx q[3];
rz(0.27942785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0937097) q[2];
sx q[2];
rz(-1.9423395) q[2];
sx q[2];
rz(2.1236146) q[2];
rz(1.7300946) q[3];
sx q[3];
rz(-2.394815) q[3];
sx q[3];
rz(1.5311034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3282851) q[0];
sx q[0];
rz(-1.3628553) q[0];
sx q[0];
rz(3.0174729) q[0];
rz(2.2165551) q[1];
sx q[1];
rz(-1.8412453) q[1];
sx q[1];
rz(-0.73840028) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8014073) q[0];
sx q[0];
rz(-0.75447318) q[0];
sx q[0];
rz(-1.4800001) q[0];
rz(-pi) q[1];
rz(1.1843119) q[2];
sx q[2];
rz(-2.1447561) q[2];
sx q[2];
rz(1.5476162) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55220375) q[1];
sx q[1];
rz(-2.5393746) q[1];
sx q[1];
rz(-2.7715384) q[1];
rz(-pi) q[2];
rz(2.5656126) q[3];
sx q[3];
rz(-1.4864731) q[3];
sx q[3];
rz(-2.4654441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.063252123) q[2];
sx q[2];
rz(-1.1880778) q[2];
sx q[2];
rz(1.1032907) q[2];
rz(0.25857806) q[3];
sx q[3];
rz(-2.0662112) q[3];
sx q[3];
rz(-0.28178373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89541188) q[0];
sx q[0];
rz(-0.87765944) q[0];
sx q[0];
rz(-1.3377162) q[0];
rz(-2.5216263) q[1];
sx q[1];
rz(-1.4620616) q[1];
sx q[1];
rz(-1.5043129) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1670187) q[0];
sx q[0];
rz(-1.3822379) q[0];
sx q[0];
rz(-1.855464) q[0];
rz(-pi) q[1];
rz(-0.70233924) q[2];
sx q[2];
rz(-1.5389256) q[2];
sx q[2];
rz(1.396871) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33584174) q[1];
sx q[1];
rz(-1.0066792) q[1];
sx q[1];
rz(-1.1104904) q[1];
rz(-pi) q[2];
rz(-0.83801954) q[3];
sx q[3];
rz(-1.6494284) q[3];
sx q[3];
rz(-1.7746314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.52516809) q[2];
sx q[2];
rz(-1.9889571) q[2];
sx q[2];
rz(-1.1241414) q[2];
rz(0.21640402) q[3];
sx q[3];
rz(-2.3085322) q[3];
sx q[3];
rz(-2.0719349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44806099) q[0];
sx q[0];
rz(-0.8140642) q[0];
sx q[0];
rz(2.6779209) q[0];
rz(-0.13557735) q[1];
sx q[1];
rz(-2.1453073) q[1];
sx q[1];
rz(-0.3840951) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4736872) q[0];
sx q[0];
rz(-2.0376922) q[0];
sx q[0];
rz(1.4422732) q[0];
rz(2.3771387) q[2];
sx q[2];
rz(-1.7181509) q[2];
sx q[2];
rz(-2.5558228) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0277703) q[1];
sx q[1];
rz(-0.40317391) q[1];
sx q[1];
rz(-1.1330963) q[1];
rz(-1.4469524) q[3];
sx q[3];
rz(-2.5899124) q[3];
sx q[3];
rz(-1.967103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7163081) q[2];
sx q[2];
rz(-0.3766489) q[2];
sx q[2];
rz(-2.6728805) q[2];
rz(-2.5675755) q[3];
sx q[3];
rz(-1.6711957) q[3];
sx q[3];
rz(-0.65042692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18844093) q[0];
sx q[0];
rz(-1.783239) q[0];
sx q[0];
rz(0.69001946) q[0];
rz(-2.4192339) q[1];
sx q[1];
rz(-1.1411618) q[1];
sx q[1];
rz(0.77675995) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1952161) q[0];
sx q[0];
rz(-1.6456983) q[0];
sx q[0];
rz(-2.9086065) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6138814) q[2];
sx q[2];
rz(-0.65831682) q[2];
sx q[2];
rz(2.0331613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1550786) q[1];
sx q[1];
rz(-2.6342794) q[1];
sx q[1];
rz(2.581937) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2508936) q[3];
sx q[3];
rz(-1.5855224) q[3];
sx q[3];
rz(-0.17531987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7338099) q[2];
sx q[2];
rz(-1.5747384) q[2];
sx q[2];
rz(-3.0877647) q[2];
rz(2.5194061) q[3];
sx q[3];
rz(-0.48431188) q[3];
sx q[3];
rz(-3.0622862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6513885) q[0];
sx q[0];
rz(-1.1505928) q[0];
sx q[0];
rz(1.7405317) q[0];
rz(0.98848629) q[1];
sx q[1];
rz(-1.8732312) q[1];
sx q[1];
rz(0.40774694) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893129) q[0];
sx q[0];
rz(-1.3527217) q[0];
sx q[0];
rz(-0.41357354) q[0];
rz(-pi) q[1];
rz(-0.12699126) q[2];
sx q[2];
rz(-2.2417129) q[2];
sx q[2];
rz(0.95580703) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0474075) q[1];
sx q[1];
rz(-2.2360206) q[1];
sx q[1];
rz(0.96324222) q[1];
rz(-2.678773) q[3];
sx q[3];
rz(-0.30202391) q[3];
sx q[3];
rz(-2.8800396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72655788) q[2];
sx q[2];
rz(-2.0099535) q[2];
sx q[2];
rz(0.69990194) q[2];
rz(-2.0578201) q[3];
sx q[3];
rz(-2.2743069) q[3];
sx q[3];
rz(2.9054902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.20339762) q[0];
sx q[0];
rz(-1.45881) q[0];
sx q[0];
rz(-0.38988018) q[0];
rz(1.9436721) q[1];
sx q[1];
rz(-0.094466297) q[1];
sx q[1];
rz(0.87497154) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.724204) q[0];
sx q[0];
rz(-1.5932788) q[0];
sx q[0];
rz(-2.1940986) q[0];
x q[1];
rz(2.0535371) q[2];
sx q[2];
rz(-0.74626479) q[2];
sx q[2];
rz(-1.7604699) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3564878) q[1];
sx q[1];
rz(-2.572959) q[1];
sx q[1];
rz(-2.4604753) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1853774) q[3];
sx q[3];
rz(-0.66171911) q[3];
sx q[3];
rz(-2.1780518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8209057) q[2];
sx q[2];
rz(-0.69183886) q[2];
sx q[2];
rz(-2.4917277) q[2];
rz(1.1167022) q[3];
sx q[3];
rz(-1.2973123) q[3];
sx q[3];
rz(0.19702774) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7277302) q[0];
sx q[0];
rz(-0.05302269) q[0];
sx q[0];
rz(-0.19456385) q[0];
rz(2.3692756) q[1];
sx q[1];
rz(-1.1338502) q[1];
sx q[1];
rz(-0.17759855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1188162) q[0];
sx q[0];
rz(-1.5701862) q[0];
sx q[0];
rz(1.5724036) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3760711) q[2];
sx q[2];
rz(-2.3280848) q[2];
sx q[2];
rz(-0.080527079) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0655443) q[1];
sx q[1];
rz(-2.0418344) q[1];
sx q[1];
rz(0.7008497) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94915139) q[3];
sx q[3];
rz(-2.2165934) q[3];
sx q[3];
rz(-0.27816712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3557768) q[2];
sx q[2];
rz(-2.5444701) q[2];
sx q[2];
rz(0.43006483) q[2];
rz(-1.1108584) q[3];
sx q[3];
rz(-1.6836124) q[3];
sx q[3];
rz(-2.7250169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149902) q[0];
sx q[0];
rz(-1.5530598) q[0];
sx q[0];
rz(-0.1834827) q[0];
rz(-1.1732187) q[1];
sx q[1];
rz(-1.219974) q[1];
sx q[1];
rz(-3.0164607) q[1];
rz(-1.0452034) q[2];
sx q[2];
rz(-1.3202616) q[2];
sx q[2];
rz(3.1332982) q[2];
rz(-0.96121721) q[3];
sx q[3];
rz(-2.3280795) q[3];
sx q[3];
rz(2.1272492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
