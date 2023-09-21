OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(-2.9601331) q[0];
rz(2.060086) q[1];
sx q[1];
rz(-0.67343155) q[1];
sx q[1];
rz(-1.0531309) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15381972) q[0];
sx q[0];
rz(-1.1115371) q[0];
sx q[0];
rz(-2.905373) q[0];
rz(0.40197576) q[2];
sx q[2];
rz(-0.66265124) q[2];
sx q[2];
rz(-1.0647917) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3095113) q[1];
sx q[1];
rz(-2.32825) q[1];
sx q[1];
rz(-1.3002212) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1374723) q[3];
sx q[3];
rz(-0.37131272) q[3];
sx q[3];
rz(2.5759047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16333214) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(1.367761) q[2];
rz(1.0129499) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(0.071454123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.0435836) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(2.5464771) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(-1.5140623) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13578116) q[0];
sx q[0];
rz(-1.0732871) q[0];
sx q[0];
rz(0.46732975) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83346955) q[2];
sx q[2];
rz(-1.5499299) q[2];
sx q[2];
rz(-1.2534864) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6914312) q[1];
sx q[1];
rz(-1.7427923) q[1];
sx q[1];
rz(1.7683692) q[1];
rz(-1.4684832) q[3];
sx q[3];
rz(-1.6204837) q[3];
sx q[3];
rz(-1.697584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65511584) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(0.79616037) q[2];
rz(0.97186175) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(-2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765091) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(3.0774975) q[0];
rz(2.8308716) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(1.6832738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1761988) q[0];
sx q[0];
rz(-1.6178693) q[0];
sx q[0];
rz(1.7182299) q[0];
rz(-pi) q[1];
rz(2.4887423) q[2];
sx q[2];
rz(-1.0703501) q[2];
sx q[2];
rz(-0.99393883) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7160733) q[1];
sx q[1];
rz(-0.41285535) q[1];
sx q[1];
rz(3.0070261) q[1];
x q[2];
rz(-1.5879257) q[3];
sx q[3];
rz(-1.2560085) q[3];
sx q[3];
rz(-1.1312315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98214275) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(1.7689765) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(-0.94846559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(1.1244208) q[0];
rz(-1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(1.6569998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35614466) q[0];
sx q[0];
rz(-1.6011642) q[0];
sx q[0];
rz(1.6214451) q[0];
x q[1];
rz(-2.0882656) q[2];
sx q[2];
rz(-1.7382858) q[2];
sx q[2];
rz(0.53211624) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.209219) q[1];
sx q[1];
rz(-1.1353496) q[1];
sx q[1];
rz(-2.2030764) q[1];
rz(-pi) q[2];
rz(-1.2403537) q[3];
sx q[3];
rz(-2.2664865) q[3];
sx q[3];
rz(-2.2479923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0174039) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(-2.1253288) q[2];
rz(-1.4034363) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(1.0884292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50773412) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(2.741709) q[0];
rz(-1.1625066) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-0.17366017) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7090209) q[0];
sx q[0];
rz(-3.0637494) q[0];
sx q[0];
rz(2.7709333) q[0];
x q[1];
rz(2.747614) q[2];
sx q[2];
rz(-1.7658965) q[2];
sx q[2];
rz(-0.56001284) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.793321) q[1];
sx q[1];
rz(-1.6760567) q[1];
sx q[1];
rz(-3.1131016) q[1];
rz(-pi) q[2];
rz(-0.33186121) q[3];
sx q[3];
rz(-1.3734986) q[3];
sx q[3];
rz(0.33533898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48012039) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(-2.373467) q[2];
rz(0.85401946) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1059234) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.4703898) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(-1.1434198) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9914046) q[0];
sx q[0];
rz(-0.17230573) q[0];
sx q[0];
rz(-0.50921391) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5383699) q[2];
sx q[2];
rz(-1.2957186) q[2];
sx q[2];
rz(0.075866931) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1089576) q[1];
sx q[1];
rz(-2.7320478) q[1];
sx q[1];
rz(1.1213379) q[1];
rz(-pi) q[2];
rz(2.8361736) q[3];
sx q[3];
rz(-1.050204) q[3];
sx q[3];
rz(-0.49314317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4986971) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(-1.6112304) q[2];
rz(-1.4536084) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(-2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(3.0221508) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(-1.4402333) q[0];
rz(0.72921905) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(-1.1332606) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2204809) q[0];
sx q[0];
rz(-2.8163914) q[0];
sx q[0];
rz(-0.13334206) q[0];
rz(-pi) q[1];
rz(1.5845756) q[2];
sx q[2];
rz(-0.85859495) q[2];
sx q[2];
rz(-2.6028002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1843695) q[1];
sx q[1];
rz(-1.8436699) q[1];
sx q[1];
rz(2.2437614) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4134737) q[3];
sx q[3];
rz(-1.8295349) q[3];
sx q[3];
rz(1.3820005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40522727) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-2.6531632) q[2];
rz(1.3119665) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(-2.5002938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064780386) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(0.07117614) q[0];
rz(-0.03216234) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(1.9326928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29584822) q[0];
sx q[0];
rz(-0.83202067) q[0];
sx q[0];
rz(2.8315298) q[0];
x q[1];
rz(-2.2279943) q[2];
sx q[2];
rz(-0.29619869) q[2];
sx q[2];
rz(-1.9343455) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.68122411) q[1];
sx q[1];
rz(-1.7011233) q[1];
sx q[1];
rz(-3.1254915) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62187059) q[3];
sx q[3];
rz(-1.8727881) q[3];
sx q[3];
rz(2.5627476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.20748392) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(-1.570638) q[2];
rz(0.87336826) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(-2.7895555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97380012) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(-0.31162509) q[0];
rz(0.82178003) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(1.6315546) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8845997) q[0];
sx q[0];
rz(-1.25367) q[0];
sx q[0];
rz(-1.3604128) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0636343) q[2];
sx q[2];
rz(-0.38804752) q[2];
sx q[2];
rz(1.0840814) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9540625) q[1];
sx q[1];
rz(-0.64779753) q[1];
sx q[1];
rz(-2.8026583) q[1];
rz(-2.4771677) q[3];
sx q[3];
rz(-2.8492152) q[3];
sx q[3];
rz(-1.2053306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8682378) q[2];
sx q[2];
rz(-2.6028825) q[2];
sx q[2];
rz(-2.3256425) q[2];
rz(-0.50968918) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(-1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508535) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(-0.2302641) q[0];
rz(0.62581217) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(-0.65840107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31357665) q[0];
sx q[0];
rz(-2.1470634) q[0];
sx q[0];
rz(-2.1616031) q[0];
rz(-pi) q[1];
rz(1.9071104) q[2];
sx q[2];
rz(-2.2173777) q[2];
sx q[2];
rz(-2.1474349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.534348) q[1];
sx q[1];
rz(-1.4048647) q[1];
sx q[1];
rz(-2.0545309) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9202616) q[3];
sx q[3];
rz(-2.0798827) q[3];
sx q[3];
rz(-1.5290608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.6752424) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(-2.8038483) q[2];
rz(-1.0234458) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52453775) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(1.2383923) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(-0.94454371) q[2];
sx q[2];
rz(-1.7164451) q[2];
sx q[2];
rz(-3.0838983) q[2];
rz(3.0425439) q[3];
sx q[3];
rz(-2.5669813) q[3];
sx q[3];
rz(3.0909227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
