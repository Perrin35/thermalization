OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(-0.87061849) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(2.9918848) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4275442) q[0];
sx q[0];
rz(-0.85435003) q[0];
sx q[0];
rz(0.9057522) q[0];
x q[1];
rz(-2.3157273) q[2];
sx q[2];
rz(-0.68544938) q[2];
sx q[2];
rz(-2.4862188) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.183061) q[1];
sx q[1];
rz(-1.5325938) q[1];
sx q[1];
rz(-1.7512291) q[1];
rz(-pi) q[2];
rz(-1.0204131) q[3];
sx q[3];
rz(-1.7341511) q[3];
sx q[3];
rz(-2.6640716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.25847882) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(2.4374938) q[2];
rz(0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(0.63252226) q[0];
rz(-2.6951492) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(2.4893563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61440496) q[0];
sx q[0];
rz(-1.5831278) q[0];
sx q[0];
rz(3.1129818) q[0];
rz(-pi) q[1];
rz(-2.2488238) q[2];
sx q[2];
rz(-0.36913482) q[2];
sx q[2];
rz(-0.54112753) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3995041) q[1];
sx q[1];
rz(-1.3474476) q[1];
sx q[1];
rz(-0.29466596) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2324123) q[3];
sx q[3];
rz(-1.0717857) q[3];
sx q[3];
rz(-0.0024851174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5941045) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(-1.9799505) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(-1.6903711) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0911672) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(-2.2913349) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(1.7920378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36950668) q[0];
sx q[0];
rz(-1.8096576) q[0];
sx q[0];
rz(2.4580965) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5047968) q[2];
sx q[2];
rz(-1.608035) q[2];
sx q[2];
rz(2.6413692) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.459356) q[1];
sx q[1];
rz(-0.56486928) q[1];
sx q[1];
rz(-2.6911246) q[1];
rz(-pi) q[2];
rz(3.1317741) q[3];
sx q[3];
rz(-1.1447284) q[3];
sx q[3];
rz(-0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(-2.8404625) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(-1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(1.7310463) q[0];
rz(-0.63181216) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(-3.1052123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8481537) q[0];
sx q[0];
rz(-1.7833976) q[0];
sx q[0];
rz(0.92377499) q[0];
x q[1];
rz(-1.8965917) q[2];
sx q[2];
rz(-2.1278283) q[2];
sx q[2];
rz(-1.4404802) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.034678) q[1];
sx q[1];
rz(-0.87768302) q[1];
sx q[1];
rz(0.17130674) q[1];
rz(-pi) q[2];
rz(-2.9070204) q[3];
sx q[3];
rz(-1.4751225) q[3];
sx q[3];
rz(-2.5535339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2146384) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(1.1882163) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(-0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.7017512) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(-0.23434815) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4011824) q[0];
sx q[0];
rz(-0.64490841) q[0];
sx q[0];
rz(2.5572204) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.064803) q[2];
sx q[2];
rz(-1.4865808) q[2];
sx q[2];
rz(-2.6780724) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.105141) q[1];
sx q[1];
rz(-2.3804133) q[1];
sx q[1];
rz(-0.21703227) q[1];
rz(-0.078951051) q[3];
sx q[3];
rz(-1.2478932) q[3];
sx q[3];
rz(-2.3372646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4328737) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(-0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7261312) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-2.3244526) q[0];
rz(0.56610402) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(-1.9979427) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4139347) q[0];
sx q[0];
rz(-1.9415138) q[0];
sx q[0];
rz(1.6137705) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32240378) q[2];
sx q[2];
rz(-2.44256) q[2];
sx q[2];
rz(-2.2088745) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3559349) q[1];
sx q[1];
rz(-0.4520843) q[1];
sx q[1];
rz(-1.7788586) q[1];
x q[2];
rz(-2.6753622) q[3];
sx q[3];
rz(-0.56785184) q[3];
sx q[3];
rz(-3.1371491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75366655) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(2.6521519) q[2];
rz(-2.9135381) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5722826) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(1.9082665) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39499261) q[0];
sx q[0];
rz(-1.1867503) q[0];
sx q[0];
rz(0.011944255) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85031021) q[2];
sx q[2];
rz(-1.0734953) q[2];
sx q[2];
rz(-0.26956272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8334956) q[1];
sx q[1];
rz(-1.6819681) q[1];
sx q[1];
rz(1.584946) q[1];
rz(-pi) q[2];
x q[2];
rz(2.535378) q[3];
sx q[3];
rz(-2.5698235) q[3];
sx q[3];
rz(-1.1796463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.037501637) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(0.070090381) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(-2.0751374) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(1.4455147) q[0];
rz(0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(1.8833556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90342605) q[0];
sx q[0];
rz(-1.4362207) q[0];
sx q[0];
rz(-1.1358791) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6446722) q[2];
sx q[2];
rz(-2.1053227) q[2];
sx q[2];
rz(-2.1246186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82842365) q[1];
sx q[1];
rz(-1.0272044) q[1];
sx q[1];
rz(1.4491175) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6120841) q[3];
sx q[3];
rz(-2.2141075) q[3];
sx q[3];
rz(-0.95638004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(0.56274596) q[2];
rz(0.19966666) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0257618) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(2.8905706) q[0];
rz(0.42731467) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(0.02773157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4902892) q[0];
sx q[0];
rz(-2.2130744) q[0];
sx q[0];
rz(2.519033) q[0];
rz(-pi) q[1];
rz(-1.8856144) q[2];
sx q[2];
rz(-2.417649) q[2];
sx q[2];
rz(0.9418504) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8393644) q[1];
sx q[1];
rz(-1.6703969) q[1];
sx q[1];
rz(-0.60217963) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5042217) q[3];
sx q[3];
rz(-1.5449617) q[3];
sx q[3];
rz(-1.0915826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1256844) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(-1.1431747) q[2];
rz(-2.9987191) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.7509572) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(2.6877158) q[0];
rz(2.4699396) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(-0.25751105) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89850241) q[0];
sx q[0];
rz(-1.9542964) q[0];
sx q[0];
rz(-0.48454185) q[0];
rz(-0.54458877) q[2];
sx q[2];
rz(-2.0047744) q[2];
sx q[2];
rz(-1.7314272) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35044893) q[1];
sx q[1];
rz(-2.2669683) q[1];
sx q[1];
rz(-1.552156) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1274687) q[3];
sx q[3];
rz(-0.52237288) q[3];
sx q[3];
rz(2.4362107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8132849) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(0.60662398) q[2];
rz(-0.47484067) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941147) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(-0.22656245) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(0.67129927) q[2];
sx q[2];
rz(-0.57325267) q[2];
sx q[2];
rz(-0.43043955) q[2];
rz(0.26279454) q[3];
sx q[3];
rz(-2.0398718) q[3];
sx q[3];
rz(-0.3077988) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];