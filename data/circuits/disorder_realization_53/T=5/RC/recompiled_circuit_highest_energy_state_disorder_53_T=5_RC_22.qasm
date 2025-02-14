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
rz(-0.1102912) q[0];
sx q[0];
rz(4.9977367) q[0];
sx q[0];
rz(9.1060299) q[0];
rz(-1.9593852) q[1];
sx q[1];
rz(-1.477834) q[1];
sx q[1];
rz(1.1629265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6642374) q[0];
sx q[0];
rz(-1.2384982) q[0];
sx q[0];
rz(0.24881451) q[0];
rz(-1.5561582) q[2];
sx q[2];
rz(-1.7649289) q[2];
sx q[2];
rz(3.0513939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0434462) q[1];
sx q[1];
rz(-2.7061228) q[1];
sx q[1];
rz(0.79186073) q[1];
x q[2];
rz(-1.5507023) q[3];
sx q[3];
rz(-0.99599518) q[3];
sx q[3];
rz(-1.8679096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9691539) q[2];
sx q[2];
rz(-2.6814851) q[2];
sx q[2];
rz(3.0058506) q[2];
rz(-2.8067449) q[3];
sx q[3];
rz(-1.2232774) q[3];
sx q[3];
rz(-0.39723435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44884509) q[0];
sx q[0];
rz(-2.1362342) q[0];
sx q[0];
rz(-2.1951065) q[0];
rz(1.954156) q[1];
sx q[1];
rz(-0.58381909) q[1];
sx q[1];
rz(-0.28396398) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0683354) q[0];
sx q[0];
rz(-0.26480246) q[0];
sx q[0];
rz(-0.7629443) q[0];
rz(-pi) q[1];
rz(-1.8503485) q[2];
sx q[2];
rz(-0.7182622) q[2];
sx q[2];
rz(0.58092434) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52580035) q[1];
sx q[1];
rz(-0.4471356) q[1];
sx q[1];
rz(-2.8195803) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7534629) q[3];
sx q[3];
rz(-1.2468657) q[3];
sx q[3];
rz(-2.7517954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2076608) q[2];
sx q[2];
rz(-1.8457103) q[2];
sx q[2];
rz(0.80061039) q[2];
rz(-1.8396395) q[3];
sx q[3];
rz(-2.0079565) q[3];
sx q[3];
rz(0.058066644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12450739) q[0];
sx q[0];
rz(-2.7094816) q[0];
sx q[0];
rz(-2.3486163) q[0];
rz(-2.4555581) q[1];
sx q[1];
rz(-1.425309) q[1];
sx q[1];
rz(-2.3763903) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25178037) q[0];
sx q[0];
rz(-1.7818101) q[0];
sx q[0];
rz(3.0821256) q[0];
rz(-pi) q[1];
rz(-3.0900218) q[2];
sx q[2];
rz(-1.602939) q[2];
sx q[2];
rz(2.3598537) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1204234) q[1];
sx q[1];
rz(-2.3115054) q[1];
sx q[1];
rz(-2.7629025) q[1];
rz(-pi) q[2];
rz(-0.92242494) q[3];
sx q[3];
rz(-0.88705813) q[3];
sx q[3];
rz(-2.2278663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5746295) q[2];
sx q[2];
rz(-2.4506863) q[2];
sx q[2];
rz(-1.2464657) q[2];
rz(1.9555107) q[3];
sx q[3];
rz(-1.7639152) q[3];
sx q[3];
rz(-2.932909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74801159) q[0];
sx q[0];
rz(-1.3205386) q[0];
sx q[0];
rz(-3.0082974) q[0];
rz(2.8721299) q[1];
sx q[1];
rz(-0.91454426) q[1];
sx q[1];
rz(0.46636137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8036706) q[0];
sx q[0];
rz(-0.9804157) q[0];
sx q[0];
rz(1.4268111) q[0];
rz(-pi) q[1];
rz(-1.3106636) q[2];
sx q[2];
rz(-2.4974653) q[2];
sx q[2];
rz(-2.3144503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.129648) q[1];
sx q[1];
rz(-2.5135871) q[1];
sx q[1];
rz(2.8824214) q[1];
rz(-pi) q[2];
rz(-0.5138054) q[3];
sx q[3];
rz(-1.6282968) q[3];
sx q[3];
rz(2.1138885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.93126297) q[2];
sx q[2];
rz(-0.42079058) q[2];
sx q[2];
rz(-0.24410625) q[2];
rz(1.3611475) q[3];
sx q[3];
rz(-2.1972392) q[3];
sx q[3];
rz(1.9349499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162783) q[0];
sx q[0];
rz(-2.2940574) q[0];
sx q[0];
rz(0.086294802) q[0];
rz(1.7591954) q[1];
sx q[1];
rz(-1.0603797) q[1];
sx q[1];
rz(1.28654) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4327398) q[0];
sx q[0];
rz(-1.7403495) q[0];
sx q[0];
rz(1.9643334) q[0];
x q[1];
rz(-2.3628144) q[2];
sx q[2];
rz(-1.4113103) q[2];
sx q[2];
rz(-2.3159112) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4205192) q[1];
sx q[1];
rz(-2.2700078) q[1];
sx q[1];
rz(-0.10743227) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6161353) q[3];
sx q[3];
rz(-2.1920929) q[3];
sx q[3];
rz(2.3182825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4061819) q[2];
sx q[2];
rz(-1.4124796) q[2];
sx q[2];
rz(-2.0162876) q[2];
rz(-0.79648894) q[3];
sx q[3];
rz(-2.4060566) q[3];
sx q[3];
rz(-2.9430732) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35561246) q[0];
sx q[0];
rz(-0.2934083) q[0];
sx q[0];
rz(0.32817131) q[0];
rz(-0.38901058) q[1];
sx q[1];
rz(-1.5704472) q[1];
sx q[1];
rz(-1.4555812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7195014) q[0];
sx q[0];
rz(-1.4900643) q[0];
sx q[0];
rz(-2.9305365) q[0];
x q[1];
rz(-2.3244072) q[2];
sx q[2];
rz(-2.9574988) q[2];
sx q[2];
rz(0.56841401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38992369) q[1];
sx q[1];
rz(-1.3397386) q[1];
sx q[1];
rz(0.91737813) q[1];
rz(-pi) q[2];
rz(-1.2120655) q[3];
sx q[3];
rz(-1.609612) q[3];
sx q[3];
rz(-1.2286548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.97633156) q[2];
sx q[2];
rz(-0.92504048) q[2];
sx q[2];
rz(-0.46249214) q[2];
rz(2.1180604) q[3];
sx q[3];
rz(-2.0868802) q[3];
sx q[3];
rz(2.3033477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3384712) q[0];
sx q[0];
rz(-1.4305038) q[0];
sx q[0];
rz(2.9679003) q[0];
rz(-2.9202785) q[1];
sx q[1];
rz(-2.4559655) q[1];
sx q[1];
rz(1.7783222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6741329) q[0];
sx q[0];
rz(-1.9741575) q[0];
sx q[0];
rz(2.9441071) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4837711) q[2];
sx q[2];
rz(-1.489893) q[2];
sx q[2];
rz(-1.9593585) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0371907) q[1];
sx q[1];
rz(-0.96417448) q[1];
sx q[1];
rz(-0.9888852) q[1];
rz(-pi) q[2];
rz(-0.24226455) q[3];
sx q[3];
rz(-1.8349832) q[3];
sx q[3];
rz(1.3380877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26019105) q[2];
sx q[2];
rz(-1.1611725) q[2];
sx q[2];
rz(1.6960404) q[2];
rz(-0.77406231) q[3];
sx q[3];
rz(-0.71662199) q[3];
sx q[3];
rz(1.6707481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5792907) q[0];
sx q[0];
rz(-0.93920541) q[0];
sx q[0];
rz(-2.8939409) q[0];
rz(0.66894764) q[1];
sx q[1];
rz(-1.9525783) q[1];
sx q[1];
rz(-2.155969) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7108817) q[0];
sx q[0];
rz(-2.1257002) q[0];
sx q[0];
rz(1.4355833) q[0];
rz(-0.47565461) q[2];
sx q[2];
rz(-0.61810571) q[2];
sx q[2];
rz(-2.3007768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0528763) q[1];
sx q[1];
rz(-0.50584882) q[1];
sx q[1];
rz(0.90152503) q[1];
rz(-0.59496112) q[3];
sx q[3];
rz(-1.56968) q[3];
sx q[3];
rz(-3.0240455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.140427) q[2];
sx q[2];
rz(-0.75158921) q[2];
sx q[2];
rz(-1.2933732) q[2];
rz(-1.2398531) q[3];
sx q[3];
rz(-0.69609061) q[3];
sx q[3];
rz(2.4333439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86466113) q[0];
sx q[0];
rz(-1.7129352) q[0];
sx q[0];
rz(0.96555936) q[0];
rz(-2.7652265) q[1];
sx q[1];
rz(-1.3220359) q[1];
sx q[1];
rz(-1.2300864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7853785) q[0];
sx q[0];
rz(-1.5206771) q[0];
sx q[0];
rz(-2.8756407) q[0];
rz(1.240991) q[2];
sx q[2];
rz(-1.9422741) q[2];
sx q[2];
rz(0.35201752) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9113637) q[1];
sx q[1];
rz(-1.6411643) q[1];
sx q[1];
rz(-0.60649782) q[1];
rz(0.053080245) q[3];
sx q[3];
rz(-0.81392589) q[3];
sx q[3];
rz(-2.2102714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85912117) q[2];
sx q[2];
rz(-2.4716447) q[2];
sx q[2];
rz(-0.13151375) q[2];
rz(0.48111835) q[3];
sx q[3];
rz(-2.2115464) q[3];
sx q[3];
rz(0.49829811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1199101) q[0];
sx q[0];
rz(-0.57761884) q[0];
sx q[0];
rz(-2.6525894) q[0];
rz(1.8148212) q[1];
sx q[1];
rz(-1.2811456) q[1];
sx q[1];
rz(-0.16924032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6833333) q[0];
sx q[0];
rz(-1.1348327) q[0];
sx q[0];
rz(0.13157121) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6851114) q[2];
sx q[2];
rz(-2.940753) q[2];
sx q[2];
rz(-2.0690837) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0161017) q[1];
sx q[1];
rz(-0.91052848) q[1];
sx q[1];
rz(-1.0339526) q[1];
rz(-pi) q[2];
rz(-2.333913) q[3];
sx q[3];
rz(-2.7480304) q[3];
sx q[3];
rz(-0.53787947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5883611) q[2];
sx q[2];
rz(-1.0855805) q[2];
sx q[2];
rz(0.21151839) q[2];
rz(1.616098) q[3];
sx q[3];
rz(-1.0001405) q[3];
sx q[3];
rz(2.8741527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3612628) q[0];
sx q[0];
rz(-1.4365256) q[0];
sx q[0];
rz(-2.8509675) q[0];
rz(2.3296539) q[1];
sx q[1];
rz(-2.5135136) q[1];
sx q[1];
rz(1.5807349) q[1];
rz(2.7079034) q[2];
sx q[2];
rz(-1.5530752) q[2];
sx q[2];
rz(-1.8032522) q[2];
rz(2.3216861) q[3];
sx q[3];
rz(-1.0884566) q[3];
sx q[3];
rz(1.4940445) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
