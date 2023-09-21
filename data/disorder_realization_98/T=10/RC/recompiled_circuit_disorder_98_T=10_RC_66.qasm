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
rz(2.2709742) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71404845) q[0];
sx q[0];
rz(-2.2872426) q[0];
sx q[0];
rz(-2.2358405) q[0];
rz(-pi) q[1];
rz(-0.50617354) q[2];
sx q[2];
rz(-2.0548327) q[2];
sx q[2];
rz(1.5278221) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7357199) q[1];
sx q[1];
rz(-2.957203) q[1];
sx q[1];
rz(-1.3609481) q[1];
rz(-pi) q[2];
x q[2];
rz(1.876086) q[3];
sx q[3];
rz(-2.5698834) q[3];
sx q[3];
rz(1.7892464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.25847882) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(0.70409888) q[2];
rz(0.95300931) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(-2.5090704) q[0];
rz(2.6951492) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(0.65223637) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1855542) q[0];
sx q[0];
rz(-1.599405) q[0];
sx q[0];
rz(1.5584598) q[0];
rz(-pi) q[1];
rz(1.863443) q[2];
sx q[2];
rz(-1.3424982) q[2];
sx q[2];
rz(-1.6739068) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3398509) q[1];
sx q[1];
rz(-0.36777126) q[1];
sx q[1];
rz(2.4778609) q[1];
rz(-pi) q[2];
rz(0.84516256) q[3];
sx q[3];
rz(-0.80544986) q[3];
sx q[3];
rz(-2.1243387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5474881) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(1.15796) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(1.6903711) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0911672) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(2.2913349) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.3495548) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.772086) q[0];
sx q[0];
rz(-1.8096576) q[0];
sx q[0];
rz(2.4580965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.037319855) q[2];
sx q[2];
rz(-1.5048426) q[2];
sx q[2];
rz(2.0685591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2022484) q[1];
sx q[1];
rz(-1.0679686) q[1];
sx q[1];
rz(1.3015675) q[1];
rz(0.0098185929) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(-0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(0.30113014) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49144739) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(0.63181216) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(0.036380336) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2934389) q[0];
sx q[0];
rz(-1.358195) q[0];
sx q[0];
rz(0.92377499) q[0];
rz(-pi) q[1];
rz(1.245001) q[2];
sx q[2];
rz(-2.1278283) q[2];
sx q[2];
rz(1.7011124) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1069146) q[1];
sx q[1];
rz(-0.87768302) q[1];
sx q[1];
rz(-0.17130674) q[1];
rz(-1.6691469) q[3];
sx q[3];
rz(-1.8042759) q[3];
sx q[3];
rz(-2.1360306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(-2.4711117) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(-2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(-2.3855551) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(2.9072445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74041022) q[0];
sx q[0];
rz(-2.4966842) q[0];
sx q[0];
rz(0.58437225) q[0];
rz(-pi) q[1];
rz(1.7469823) q[2];
sx q[2];
rz(-0.50054769) q[2];
sx q[2];
rz(1.2622152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.036451642) q[1];
sx q[1];
rz(-0.76117939) q[1];
sx q[1];
rz(0.21703227) q[1];
rz(1.8022728) q[3];
sx q[3];
rz(-0.33208648) q[3];
sx q[3];
rz(-2.092923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4328737) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(2.7045414) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(-0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7261312) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(2.3244526) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.9979427) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691539) q[0];
sx q[0];
rz(-1.6108496) q[0];
sx q[0];
rz(-2.7705631) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67310682) q[2];
sx q[2];
rz(-1.7761201) q[2];
sx q[2];
rz(-0.38773195) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3559349) q[1];
sx q[1];
rz(-2.6895084) q[1];
sx q[1];
rz(1.7788586) q[1];
rz(-pi) q[2];
rz(2.6236344) q[3];
sx q[3];
rz(-1.3266139) q[3];
sx q[3];
rz(1.9675919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-2.6521519) q[2];
rz(2.9135381) q[3];
sx q[3];
rz(-1.8830048) q[3];
sx q[3];
rz(-0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.56931) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(1.9082665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7147303) q[0];
sx q[0];
rz(-0.38422248) q[0];
sx q[0];
rz(-1.6003438) q[0];
rz(-pi) q[1];
rz(-2.2912824) q[2];
sx q[2];
rz(-2.0680973) q[2];
sx q[2];
rz(2.8720299) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7066321) q[1];
sx q[1];
rz(-3.0295277) q[1];
sx q[1];
rz(-3.0155165) q[1];
rz(-0.4864278) q[3];
sx q[3];
rz(-1.8842116) q[3];
sx q[3];
rz(-0.13669554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.037501637) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(1.0160149) q[2];
rz(0.070090381) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(1.0664553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(1.4455147) q[0];
rz(2.9267172) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(-1.8833556) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94851516) q[0];
sx q[0];
rz(-2.6876039) q[0];
sx q[0];
rz(1.8817188) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6446722) q[2];
sx q[2];
rz(-1.03627) q[2];
sx q[2];
rz(-2.1246186) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3360577) q[1];
sx q[1];
rz(-1.6748669) q[1];
sx q[1];
rz(-2.59471) q[1];
rz(2.2860252) q[3];
sx q[3];
rz(-1.9868402) q[3];
sx q[3];
rz(2.1895529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.68226472) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(2.941926) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(-2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0257618) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-0.2510221) q[0];
rz(-0.42731467) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(-0.02773157) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4902892) q[0];
sx q[0];
rz(-2.2130744) q[0];
sx q[0];
rz(0.62255967) q[0];
rz(-pi) q[1];
rz(-1.8856144) q[2];
sx q[2];
rz(-0.72394365) q[2];
sx q[2];
rz(2.1997423) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8393644) q[1];
sx q[1];
rz(-1.6703969) q[1];
sx q[1];
rz(2.539413) q[1];
rz(-pi) q[2];
rz(-0.025891993) q[3];
sx q[3];
rz(-1.6373487) q[3];
sx q[3];
rz(2.6641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(-1.1431747) q[2];
rz(2.9987191) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(-2.8840816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47782183) q[0];
sx q[0];
rz(-2.0174694) q[0];
sx q[0];
rz(1.1429943) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5970039) q[2];
sx q[2];
rz(-1.1368183) q[2];
sx q[2];
rz(1.4101654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.232302) q[1];
sx q[1];
rz(-1.5850987) q[1];
sx q[1];
rz(-0.69625744) q[1];
rz(-pi) q[2];
rz(0.24210838) q[3];
sx q[3];
rz(-1.1032411) q[3];
sx q[3];
rz(2.9374591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3283078) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(2.5349687) q[2];
rz(2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7941147) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
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
rz(2.0444617) q[3];
sx q[3];
rz(-0.53285014) q[3];
sx q[3];
rz(0.22900029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];