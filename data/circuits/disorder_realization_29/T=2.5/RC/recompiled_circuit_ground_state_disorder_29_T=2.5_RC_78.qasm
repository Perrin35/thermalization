OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.27958265) q[0];
sx q[0];
rz(-2.6065338) q[0];
sx q[0];
rz(2.9076599) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(-1.7276126) q[1];
sx q[1];
rz(1.6230621) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6739963) q[0];
sx q[0];
rz(-0.61450555) q[0];
sx q[0];
rz(-0.55802457) q[0];
rz(-pi) q[1];
rz(1.1351003) q[2];
sx q[2];
rz(-0.47347906) q[2];
sx q[2];
rz(-1.4896637) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17659345) q[1];
sx q[1];
rz(-1.0142583) q[1];
sx q[1];
rz(0.92265572) q[1];
rz(-pi) q[2];
rz(-2.2330448) q[3];
sx q[3];
rz(-1.9534142) q[3];
sx q[3];
rz(-0.97272129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5412377) q[2];
sx q[2];
rz(-1.4404094) q[2];
sx q[2];
rz(0.8055299) q[2];
rz(0.39189288) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(2.384757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58981744) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(-2.7031194) q[0];
rz(-1.27502) q[1];
sx q[1];
rz(-1.4507111) q[1];
sx q[1];
rz(-3.0335887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40168821) q[0];
sx q[0];
rz(-1.5956912) q[0];
sx q[0];
rz(1.6817001) q[0];
rz(-pi) q[1];
rz(-2.2453868) q[2];
sx q[2];
rz(-0.87088481) q[2];
sx q[2];
rz(0.43255478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0846363) q[1];
sx q[1];
rz(-0.26727391) q[1];
sx q[1];
rz(1.2483424) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9505902) q[3];
sx q[3];
rz(-2.0536978) q[3];
sx q[3];
rz(-2.1306899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3652304) q[2];
sx q[2];
rz(-2.6680816) q[2];
sx q[2];
rz(3.0253809) q[2];
rz(-2.727437) q[3];
sx q[3];
rz(-1.073758) q[3];
sx q[3];
rz(-2.3381086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47495833) q[0];
sx q[0];
rz(-1.3131498) q[0];
sx q[0];
rz(-2.3057002) q[0];
rz(-1.5180786) q[1];
sx q[1];
rz(-1.514879) q[1];
sx q[1];
rz(-1.9035043) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19003632) q[0];
sx q[0];
rz(-2.073003) q[0];
sx q[0];
rz(-0.42597187) q[0];
x q[1];
rz(0.61254259) q[2];
sx q[2];
rz(-0.41191891) q[2];
sx q[2];
rz(-1.4985794) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.064621335) q[1];
sx q[1];
rz(-0.13084743) q[1];
sx q[1];
rz(-2.0276643) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2830623) q[3];
sx q[3];
rz(-1.2997174) q[3];
sx q[3];
rz(0.5036186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8572924) q[2];
sx q[2];
rz(-0.14480536) q[2];
sx q[2];
rz(2.4227552) q[2];
rz(2.5607064) q[3];
sx q[3];
rz(-1.418117) q[3];
sx q[3];
rz(2.412793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.523664) q[0];
sx q[0];
rz(-1.5476462) q[0];
sx q[0];
rz(2.0148328) q[0];
rz(1.2976546) q[1];
sx q[1];
rz(-1.490386) q[1];
sx q[1];
rz(1.0430956) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27985036) q[0];
sx q[0];
rz(-1.6011097) q[0];
sx q[0];
rz(1.8413196) q[0];
x q[1];
rz(1.9773433) q[2];
sx q[2];
rz(-0.57465982) q[2];
sx q[2];
rz(2.3324147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9981737) q[1];
sx q[1];
rz(-2.7980248) q[1];
sx q[1];
rz(2.0433389) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6564441) q[3];
sx q[3];
rz(-2.2243735) q[3];
sx q[3];
rz(-1.9854922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8635233) q[2];
sx q[2];
rz(-1.7570644) q[2];
sx q[2];
rz(-1.0370022) q[2];
rz(2.1841124) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(-0.83468848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1400414) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(3.0539883) q[0];
rz(-2.7032848) q[1];
sx q[1];
rz(-1.0821082) q[1];
sx q[1];
rz(1.8278488) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3739819) q[0];
sx q[0];
rz(-1.0896519) q[0];
sx q[0];
rz(-2.8767881) q[0];
rz(-pi) q[1];
x q[1];
rz(1.472166) q[2];
sx q[2];
rz(-1.6079788) q[2];
sx q[2];
rz(1.1929026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48291656) q[1];
sx q[1];
rz(-0.19397846) q[1];
sx q[1];
rz(-2.0493755) q[1];
rz(0.99788061) q[3];
sx q[3];
rz(-2.8393203) q[3];
sx q[3];
rz(-0.59461601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6614512) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(0.56337774) q[2];
rz(-1.2207458) q[3];
sx q[3];
rz(-1.7505587) q[3];
sx q[3];
rz(0.63108546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1221984) q[0];
sx q[0];
rz(-0.65827426) q[0];
sx q[0];
rz(0.011938183) q[0];
rz(2.7519233) q[1];
sx q[1];
rz(-2.4395112) q[1];
sx q[1];
rz(-2.8964002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.747904) q[0];
sx q[0];
rz(-2.1576799) q[0];
sx q[0];
rz(-2.6176207) q[0];
x q[1];
rz(-0.25845627) q[2];
sx q[2];
rz(-1.3535796) q[2];
sx q[2];
rz(-0.9612135) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4820741) q[1];
sx q[1];
rz(-2.4122752) q[1];
sx q[1];
rz(-2.4575255) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40254503) q[3];
sx q[3];
rz(-1.3005101) q[3];
sx q[3];
rz(-2.9077934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37990722) q[2];
sx q[2];
rz(-1.3096389) q[2];
sx q[2];
rz(1.7725819) q[2];
rz(0.53064972) q[3];
sx q[3];
rz(-2.4574418) q[3];
sx q[3];
rz(-2.0556889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.948792) q[0];
sx q[0];
rz(-1.8185607) q[0];
sx q[0];
rz(0.2970933) q[0];
rz(1.5104712) q[1];
sx q[1];
rz(-2.511697) q[1];
sx q[1];
rz(0.43103257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1566136) q[0];
sx q[0];
rz(-2.7811858) q[0];
sx q[0];
rz(-2.5094338) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5181584) q[2];
sx q[2];
rz(-1.8013012) q[2];
sx q[2];
rz(-2.8090257) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7546852) q[1];
sx q[1];
rz(-1.0656989) q[1];
sx q[1];
rz(-0.37270697) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0991574) q[3];
sx q[3];
rz(-1.8260806) q[3];
sx q[3];
rz(2.0057037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.83583528) q[2];
sx q[2];
rz(-0.77360669) q[2];
sx q[2];
rz(0.26710278) q[2];
rz(-3.109572) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(3.035868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06374643) q[0];
sx q[0];
rz(-1.2910605) q[0];
sx q[0];
rz(-0.63999501) q[0];
rz(-2.2867639) q[1];
sx q[1];
rz(-1.6565485) q[1];
sx q[1];
rz(1.4124426) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0542568) q[0];
sx q[0];
rz(-1.0414413) q[0];
sx q[0];
rz(-2.0482778) q[0];
rz(0.051338685) q[2];
sx q[2];
rz(-1.079139) q[2];
sx q[2];
rz(1.7887539) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58929491) q[1];
sx q[1];
rz(-1.4865685) q[1];
sx q[1];
rz(2.8003947) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37483172) q[3];
sx q[3];
rz(-1.171798) q[3];
sx q[3];
rz(0.77038902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.60014805) q[2];
sx q[2];
rz(-1.4358127) q[2];
sx q[2];
rz(-0.42204648) q[2];
rz(0.47354928) q[3];
sx q[3];
rz(-2.519042) q[3];
sx q[3];
rz(-0.1951018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37050978) q[0];
sx q[0];
rz(-2.2042553) q[0];
sx q[0];
rz(-1.7899845) q[0];
rz(0.31556684) q[1];
sx q[1];
rz(-1.6866997) q[1];
sx q[1];
rz(-1.9884761) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2446072) q[0];
sx q[0];
rz(-0.073052064) q[0];
sx q[0];
rz(-0.79985072) q[0];
x q[1];
rz(-2.6062327) q[2];
sx q[2];
rz(-1.6125814) q[2];
sx q[2];
rz(2.1891862) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7866384) q[1];
sx q[1];
rz(-1.4435205) q[1];
sx q[1];
rz(-0.58357012) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82984583) q[3];
sx q[3];
rz(-2.6072558) q[3];
sx q[3];
rz(2.9754834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6844668) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(-0.37377629) q[2];
rz(1.2260381) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(1.3396243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95947295) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(-1.3207588) q[0];
rz(2.2044115) q[1];
sx q[1];
rz(-1.3359741) q[1];
sx q[1];
rz(0.46930596) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73399765) q[0];
sx q[0];
rz(-3.0014801) q[0];
sx q[0];
rz(-0.99262832) q[0];
x q[1];
rz(1.2861757) q[2];
sx q[2];
rz(-2.5409343) q[2];
sx q[2];
rz(0.17135581) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.093897029) q[1];
sx q[1];
rz(-0.17054955) q[1];
sx q[1];
rz(-2.6490414) q[1];
x q[2];
rz(1.3239884) q[3];
sx q[3];
rz(-1.2780407) q[3];
sx q[3];
rz(0.61652641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.70539537) q[2];
sx q[2];
rz(-2.2787978) q[2];
sx q[2];
rz(-1.2104642) q[2];
rz(-0.56810275) q[3];
sx q[3];
rz(-1.7567239) q[3];
sx q[3];
rz(-0.97584045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2198915) q[0];
sx q[0];
rz(-2.2911063) q[0];
sx q[0];
rz(-1.6794857) q[0];
rz(-2.2484491) q[1];
sx q[1];
rz(-1.8291263) q[1];
sx q[1];
rz(-1.6222454) q[1];
rz(1.2814796) q[2];
sx q[2];
rz(-1.791496) q[2];
sx q[2];
rz(1.6481177) q[2];
rz(-0.4332581) q[3];
sx q[3];
rz(-1.0436202) q[3];
sx q[3];
rz(2.9544261) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
