OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52842927) q[0];
sx q[0];
rz(2.0818721) q[0];
sx q[0];
rz(11.835397) q[0];
rz(-1.5001186) q[1];
sx q[1];
rz(-2.1067696) q[1];
sx q[1];
rz(0.94354454) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4306513) q[0];
sx q[0];
rz(-1.6057771) q[0];
sx q[0];
rz(-0.067023858) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.655683) q[2];
sx q[2];
rz(-0.92157084) q[2];
sx q[2];
rz(2.6390586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.01602068) q[1];
sx q[1];
rz(-1.9106094) q[1];
sx q[1];
rz(-1.9858951) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30480095) q[3];
sx q[3];
rz(-1.3918575) q[3];
sx q[3];
rz(1.6182871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.25508183) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(1.250766) q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0086867) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(2.5426478) q[0];
rz(-1.8006181) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(2.1751931) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7705298) q[0];
sx q[0];
rz(-2.3822228) q[0];
sx q[0];
rz(2.4735527) q[0];
x q[1];
rz(3.0722926) q[2];
sx q[2];
rz(-0.6233223) q[2];
sx q[2];
rz(-0.22340439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4233154) q[1];
sx q[1];
rz(-1.4916972) q[1];
sx q[1];
rz(3.0273816) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1411243) q[3];
sx q[3];
rz(-0.86887348) q[3];
sx q[3];
rz(0.3598635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0559343) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(0.30109626) q[2];
rz(1.1931233) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(-1.1863856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10467228) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(-1.954129) q[0];
rz(-1.2359515) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(-1.8240066) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0061958) q[0];
sx q[0];
rz(-2.4164696) q[0];
sx q[0];
rz(-0.32307415) q[0];
x q[1];
rz(1.7682398) q[2];
sx q[2];
rz(-1.3609481) q[2];
sx q[2];
rz(0.91453493) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4730075) q[1];
sx q[1];
rz(-0.43273941) q[1];
sx q[1];
rz(-0.24344484) q[1];
rz(-pi) q[2];
rz(2.0882323) q[3];
sx q[3];
rz(-1.468717) q[3];
sx q[3];
rz(1.3641588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(-2.9349566) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(-2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3669423) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(0.88622093) q[0];
rz(2.1318502) q[1];
sx q[1];
rz(-0.90615288) q[1];
sx q[1];
rz(1.9151691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3564295) q[0];
sx q[0];
rz(-2.4117081) q[0];
sx q[0];
rz(-1.4714144) q[0];
rz(-pi) q[1];
rz(-0.88422758) q[2];
sx q[2];
rz(-1.9366169) q[2];
sx q[2];
rz(-0.58931749) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1434506) q[1];
sx q[1];
rz(-1.7805903) q[1];
sx q[1];
rz(1.7721121) q[1];
rz(-pi) q[2];
rz(1.1981443) q[3];
sx q[3];
rz(-2.028392) q[3];
sx q[3];
rz(-2.3082993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6440789) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(-2.148596) q[2];
rz(1.3126866) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-0.48373568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023225697) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(-1.9556048) q[0];
rz(-1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-2.5591154) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7459481) q[0];
sx q[0];
rz(-2.1714604) q[0];
sx q[0];
rz(3.0546741) q[0];
rz(0.8930348) q[2];
sx q[2];
rz(-1.6871916) q[2];
sx q[2];
rz(-1.900577) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9936258) q[1];
sx q[1];
rz(-2.2851351) q[1];
sx q[1];
rz(1.5453969) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4062823) q[3];
sx q[3];
rz(-1.5092106) q[3];
sx q[3];
rz(-0.53698925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4867268) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(-3.0598818) q[2];
rz(2.667526) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(1.4985532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.896647) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(-0.0078049302) q[0];
rz(1.7410949) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(-2.0369464) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305599) q[0];
sx q[0];
rz(-1.7337807) q[0];
sx q[0];
rz(2.4736604) q[0];
rz(-pi) q[1];
rz(2.8203037) q[2];
sx q[2];
rz(-2.4761204) q[2];
sx q[2];
rz(1.9208391) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.066597477) q[1];
sx q[1];
rz(-2.1146333) q[1];
sx q[1];
rz(-0.72504136) q[1];
x q[2];
rz(2.036318) q[3];
sx q[3];
rz(-1.421531) q[3];
sx q[3];
rz(-2.535459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2509987) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(-0.55523038) q[2];
rz(-2.9688719) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5531439) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(-3.0798262) q[0];
rz(0.24208367) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(2.0297091) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88469807) q[0];
sx q[0];
rz(-1.9848794) q[0];
sx q[0];
rz(-2.3143682) q[0];
rz(-3.1277666) q[2];
sx q[2];
rz(-1.6484043) q[2];
sx q[2];
rz(1.3071878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55918499) q[1];
sx q[1];
rz(-1.8473986) q[1];
sx q[1];
rz(0.72699593) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63515969) q[3];
sx q[3];
rz(-0.40024647) q[3];
sx q[3];
rz(1.1349585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8093439) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(-0.81364441) q[2];
rz(1.7371197) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5161045) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(-1.425449) q[0];
rz(-1.5215993) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(2.5040748) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0901776) q[0];
sx q[0];
rz(-2.1438103) q[0];
sx q[0];
rz(-0.90419241) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39890639) q[2];
sx q[2];
rz(-2.3762694) q[2];
sx q[2];
rz(1.4590291) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4379165) q[1];
sx q[1];
rz(-2.1494467) q[1];
sx q[1];
rz(2.0520567) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5708837) q[3];
sx q[3];
rz(-1.4548886) q[3];
sx q[3];
rz(-1.7403719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0104388) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(1.1361702) q[2];
rz(-1.6561967) q[3];
sx q[3];
rz(-2.6172726) q[3];
sx q[3];
rz(1.9320528) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6256325) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(2.4654454) q[0];
rz(-2.3162084) q[1];
sx q[1];
rz(-0.4239347) q[1];
sx q[1];
rz(1.0151781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0494173) q[0];
sx q[0];
rz(-1.9855238) q[0];
sx q[0];
rz(-1.7134922) q[0];
rz(-0.19998156) q[2];
sx q[2];
rz(-2.4790384) q[2];
sx q[2];
rz(2.8560864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.059541313) q[1];
sx q[1];
rz(-0.50652981) q[1];
sx q[1];
rz(-2.3561213) q[1];
rz(-pi) q[2];
rz(2.3750651) q[3];
sx q[3];
rz(-1.3228647) q[3];
sx q[3];
rz(0.31162308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72145808) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-0.96735111) q[2];
rz(-1.5445276) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(-0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5230781) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(3.0850947) q[0];
rz(2.0691195) q[1];
sx q[1];
rz(-1.2697376) q[1];
sx q[1];
rz(-1.4046232) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4182189) q[0];
sx q[0];
rz(-1.2105816) q[0];
sx q[0];
rz(3.0035613) q[0];
rz(-0.3041515) q[2];
sx q[2];
rz(-1.6870105) q[2];
sx q[2];
rz(-0.2955557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2056634) q[1];
sx q[1];
rz(-2.231039) q[1];
sx q[1];
rz(2.5388989) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9910562) q[3];
sx q[3];
rz(-1.5794465) q[3];
sx q[3];
rz(1.0516402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.848032) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(-0.30612293) q[2];
rz(0.39811578) q[3];
sx q[3];
rz(-1.476215) q[3];
sx q[3];
rz(-1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33481471) q[0];
sx q[0];
rz(-2.0756742) q[0];
sx q[0];
rz(-2.6609127) q[0];
rz(1.1595935) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(2.2904916) q[2];
sx q[2];
rz(-2.7534178) q[2];
sx q[2];
rz(-0.22321246) q[2];
rz(-1.0032734) q[3];
sx q[3];
rz(-0.56837396) q[3];
sx q[3];
rz(3.0336998) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];