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
rz(2.3417176) q[0];
sx q[0];
rz(-2.3246111) q[0];
sx q[0];
rz(-2.6843827) q[0];
rz(0.41181052) q[1];
sx q[1];
rz(-1.5195941) q[1];
sx q[1];
rz(-2.5817459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8053539) q[0];
sx q[0];
rz(-1.8913664) q[0];
sx q[0];
rz(0.47236116) q[0];
rz(-pi) q[1];
x q[1];
rz(2.855515) q[2];
sx q[2];
rz(-2.0564579) q[2];
sx q[2];
rz(-0.81862517) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.35750439) q[1];
sx q[1];
rz(-0.58114806) q[1];
sx q[1];
rz(0.52304348) q[1];
x q[2];
rz(1.718097) q[3];
sx q[3];
rz(-1.3527186) q[3];
sx q[3];
rz(0.17373057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7655699) q[2];
sx q[2];
rz(-2.1799808) q[2];
sx q[2];
rz(1.8454856) q[2];
rz(2.5206595) q[3];
sx q[3];
rz(-2.4758078) q[3];
sx q[3];
rz(-0.92500979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6821482) q[0];
sx q[0];
rz(-2.5542673) q[0];
sx q[0];
rz(-0.71429724) q[0];
rz(-0.722305) q[1];
sx q[1];
rz(-2.0562833) q[1];
sx q[1];
rz(-1.8169656) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6178169) q[0];
sx q[0];
rz(-1.1259698) q[0];
sx q[0];
rz(2.8033048) q[0];
rz(3.0804964) q[2];
sx q[2];
rz(-2.4078566) q[2];
sx q[2];
rz(0.82759418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3199504) q[1];
sx q[1];
rz(-2.7378203) q[1];
sx q[1];
rz(-2.7266704) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0194195) q[3];
sx q[3];
rz(-0.88910149) q[3];
sx q[3];
rz(-0.091191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.033919949) q[2];
sx q[2];
rz(-1.718797) q[2];
sx q[2];
rz(2.1495492) q[2];
rz(2.3658559) q[3];
sx q[3];
rz(-0.78191596) q[3];
sx q[3];
rz(1.0824664) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7139605) q[0];
sx q[0];
rz(-1.3466703) q[0];
sx q[0];
rz(-2.2295075) q[0];
rz(0.38463792) q[1];
sx q[1];
rz(-2.1802528) q[1];
sx q[1];
rz(0.49547637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81264436) q[0];
sx q[0];
rz(-1.616713) q[0];
sx q[0];
rz(-1.6434692) q[0];
rz(-pi) q[1];
x q[1];
rz(1.762359) q[2];
sx q[2];
rz(-0.83070101) q[2];
sx q[2];
rz(-0.57410115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8532995) q[1];
sx q[1];
rz(-1.5202731) q[1];
sx q[1];
rz(0.72447296) q[1];
rz(1.552565) q[3];
sx q[3];
rz(-0.62199253) q[3];
sx q[3];
rz(-2.7014159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2428525) q[2];
sx q[2];
rz(-2.0169368) q[2];
sx q[2];
rz(-2.7070572) q[2];
rz(0.09856002) q[3];
sx q[3];
rz(-1.7193272) q[3];
sx q[3];
rz(2.3677473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74715215) q[0];
sx q[0];
rz(-2.9065865) q[0];
sx q[0];
rz(2.8727942) q[0];
rz(1.0768184) q[1];
sx q[1];
rz(-1.7348758) q[1];
sx q[1];
rz(-1.1928308) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6679316) q[0];
sx q[0];
rz(-2.8369378) q[0];
sx q[0];
rz(-2.5308035) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7014163) q[2];
sx q[2];
rz(-1.7057749) q[2];
sx q[2];
rz(-0.93739742) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0145906) q[1];
sx q[1];
rz(-1.4140097) q[1];
sx q[1];
rz(-0.21611045) q[1];
rz(-1.773473) q[3];
sx q[3];
rz(-0.87966387) q[3];
sx q[3];
rz(-2.1512669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8863135) q[2];
sx q[2];
rz(-0.89526075) q[2];
sx q[2];
rz(2.8423584) q[2];
rz(2.8507161) q[3];
sx q[3];
rz(-1.8894922) q[3];
sx q[3];
rz(2.4860184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84151477) q[0];
sx q[0];
rz(-2.1617007) q[0];
sx q[0];
rz(-0.078068659) q[0];
rz(2.1878751) q[1];
sx q[1];
rz(-0.51906145) q[1];
sx q[1];
rz(-1.5740707) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8918607) q[0];
sx q[0];
rz(-1.7513236) q[0];
sx q[0];
rz(-0.8892699) q[0];
rz(-pi) q[1];
rz(-0.6965397) q[2];
sx q[2];
rz(-1.8395506) q[2];
sx q[2];
rz(1.112325) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.138254) q[1];
sx q[1];
rz(-1.5383771) q[1];
sx q[1];
rz(0.95862548) q[1];
x q[2];
rz(-1.5007581) q[3];
sx q[3];
rz(-2.2626855) q[3];
sx q[3];
rz(-0.50229154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2786402) q[2];
sx q[2];
rz(-2.7095257) q[2];
sx q[2];
rz(-0.26818177) q[2];
rz(-2.6523318) q[3];
sx q[3];
rz(-2.0475976) q[3];
sx q[3];
rz(-1.6930273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(2.0972524) q[0];
sx q[0];
rz(-3.1179929) q[0];
sx q[0];
rz(2.6771255) q[0];
rz(-0.52344549) q[1];
sx q[1];
rz(-0.76952666) q[1];
sx q[1];
rz(-2.4109667) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8051315) q[0];
sx q[0];
rz(-1.7117654) q[0];
sx q[0];
rz(2.3192498) q[0];
x q[1];
rz(0.36412698) q[2];
sx q[2];
rz(-2.4362323) q[2];
sx q[2];
rz(0.10157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3174322) q[1];
sx q[1];
rz(-2.3329314) q[1];
sx q[1];
rz(-0.59989329) q[1];
x q[2];
rz(0.2517638) q[3];
sx q[3];
rz(-1.6986966) q[3];
sx q[3];
rz(-2.3826016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0940493) q[2];
sx q[2];
rz(-2.0619679) q[2];
sx q[2];
rz(0.13761061) q[2];
rz(-2.008647) q[3];
sx q[3];
rz(-1.9652941) q[3];
sx q[3];
rz(-1.8784116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9889744) q[0];
sx q[0];
rz(-1.512383) q[0];
sx q[0];
rz(-0.033893943) q[0];
rz(2.7821817) q[1];
sx q[1];
rz(-1.6172599) q[1];
sx q[1];
rz(0.83438897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8501668) q[0];
sx q[0];
rz(-1.8803673) q[0];
sx q[0];
rz(-0.39637027) q[0];
x q[1];
rz(1.0831725) q[2];
sx q[2];
rz(-2.1994123) q[2];
sx q[2];
rz(0.64338976) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31454489) q[1];
sx q[1];
rz(-2.8706708) q[1];
sx q[1];
rz(1.2071868) q[1];
rz(-3.062617) q[3];
sx q[3];
rz(-0.67310909) q[3];
sx q[3];
rz(0.94386327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.413201) q[2];
sx q[2];
rz(-2.8381556) q[2];
sx q[2];
rz(-1.0580074) q[2];
rz(0.47438619) q[3];
sx q[3];
rz(-2.3460903) q[3];
sx q[3];
rz(-2.0856196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0742842) q[0];
sx q[0];
rz(-1.5165167) q[0];
sx q[0];
rz(1.3577331) q[0];
rz(-0.43074295) q[1];
sx q[1];
rz(-1.6669225) q[1];
sx q[1];
rz(2.6745083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7829091) q[0];
sx q[0];
rz(-3.0984146) q[0];
sx q[0];
rz(1.5367277) q[0];
rz(-pi) q[1];
rz(-1.7632887) q[2];
sx q[2];
rz(-0.87710947) q[2];
sx q[2];
rz(-2.326593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39462806) q[1];
sx q[1];
rz(-2.10022) q[1];
sx q[1];
rz(-1.5922597) q[1];
rz(-2.2086618) q[3];
sx q[3];
rz(-1.6137505) q[3];
sx q[3];
rz(-1.651498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0908541) q[2];
sx q[2];
rz(-0.41768062) q[2];
sx q[2];
rz(1.0806855) q[2];
rz(-2.4596227) q[3];
sx q[3];
rz(-0.8050279) q[3];
sx q[3];
rz(-2.4431156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68898106) q[0];
sx q[0];
rz(-2.6230951) q[0];
sx q[0];
rz(2.3338351) q[0];
rz(1.0001596) q[1];
sx q[1];
rz(-2.2554485) q[1];
sx q[1];
rz(-0.79661405) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33045443) q[0];
sx q[0];
rz(-1.66798) q[0];
sx q[0];
rz(3.099346) q[0];
rz(-0.32795017) q[2];
sx q[2];
rz(-0.65153507) q[2];
sx q[2];
rz(-2.8095989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.12789179) q[1];
sx q[1];
rz(-1.7096448) q[1];
sx q[1];
rz(-1.9816887) q[1];
rz(-pi) q[2];
rz(-1.8851938) q[3];
sx q[3];
rz(-1.3290429) q[3];
sx q[3];
rz(0.35255656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54958582) q[2];
sx q[2];
rz(-2.0506115) q[2];
sx q[2];
rz(2.8263212) q[2];
rz(-1.573805) q[3];
sx q[3];
rz(-1.9939634) q[3];
sx q[3];
rz(1.1228336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.10298097) q[0];
sx q[0];
rz(-2.4760315) q[0];
sx q[0];
rz(0.82157201) q[0];
rz(-2.6577677) q[1];
sx q[1];
rz(-1.7211823) q[1];
sx q[1];
rz(-0.22463591) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064817682) q[0];
sx q[0];
rz(-2.7833287) q[0];
sx q[0];
rz(1.129877) q[0];
x q[1];
rz(-2.5098652) q[2];
sx q[2];
rz(-2.1251107) q[2];
sx q[2];
rz(-0.71507031) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3514631) q[1];
sx q[1];
rz(-0.7143414) q[1];
sx q[1];
rz(-1.2495561) q[1];
rz(1.0472222) q[3];
sx q[3];
rz(-2.8955799) q[3];
sx q[3];
rz(0.29073989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3839174) q[2];
sx q[2];
rz(-0.11666798) q[2];
sx q[2];
rz(0.095001027) q[2];
rz(1.4875686) q[3];
sx q[3];
rz(-0.58949685) q[3];
sx q[3];
rz(-2.3768363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8944396) q[0];
sx q[0];
rz(-2.7338487) q[0];
sx q[0];
rz(-0.26640531) q[0];
rz(-2.5447625) q[1];
sx q[1];
rz(-1.6886371) q[1];
sx q[1];
rz(1.4773038) q[1];
rz(2.0043787) q[2];
sx q[2];
rz(-2.4162393) q[2];
sx q[2];
rz(-2.2888714) q[2];
rz(-0.0048051759) q[3];
sx q[3];
rz(-2.7660696) q[3];
sx q[3];
rz(0.39733359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
