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
rz(-2.1980481) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65981728) q[0];
sx q[0];
rz(-0.075591139) q[0];
sx q[0];
rz(0.48150058) q[0];
rz(-1.4859096) q[2];
sx q[2];
rz(-2.2200218) q[2];
sx q[2];
rz(2.6390586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.01602068) q[1];
sx q[1];
rz(-1.2309832) q[1];
sx q[1];
rz(-1.9858951) q[1];
x q[2];
rz(-2.5991873) q[3];
sx q[3];
rz(-0.35202682) q[3];
sx q[3];
rz(0.56234081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(1.8908267) q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(2.1616518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0086867) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(0.59894484) q[0];
rz(-1.3409746) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(0.96639955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7705298) q[0];
sx q[0];
rz(-2.3822228) q[0];
sx q[0];
rz(-0.66803996) q[0];
rz(-pi) q[1];
rz(-0.069300058) q[2];
sx q[2];
rz(-2.5182704) q[2];
sx q[2];
rz(0.22340439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84345531) q[1];
sx q[1];
rz(-1.4569439) q[1];
sx q[1];
rz(-1.650412) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3923799) q[3];
sx q[3];
rz(-1.2470761) q[3];
sx q[3];
rz(-0.92326984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0559343) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(-0.30109626) q[2];
rz(1.1931233) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10467228) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(-1.1874636) q[0];
rz(-1.2359515) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(1.3175861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0061958) q[0];
sx q[0];
rz(-0.72512308) q[0];
sx q[0];
rz(-0.32307415) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3733528) q[2];
sx q[2];
rz(-1.7806446) q[2];
sx q[2];
rz(0.91453493) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.8981478) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1173238) q[3];
sx q[3];
rz(-1.0563207) q[3];
sx q[3];
rz(-0.26457149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1304156) q[2];
sx q[2];
rz(-1.3935564) q[2];
sx q[2];
rz(0.2066361) q[2];
rz(0.7080428) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(-2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77465039) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(-2.2553717) q[0];
rz(2.1318502) q[1];
sx q[1];
rz(-0.90615288) q[1];
sx q[1];
rz(1.9151691) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7851631) q[0];
sx q[0];
rz(-2.4117081) q[0];
sx q[0];
rz(-1.4714144) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88422758) q[2];
sx q[2];
rz(-1.2049757) q[2];
sx q[2];
rz(-0.58931749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.671771) q[1];
sx q[1];
rz(-1.3739532) q[1];
sx q[1];
rz(-0.21398869) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63663441) q[3];
sx q[3];
rz(-0.58167471) q[3];
sx q[3];
rz(3.0326774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6440789) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(-2.148596) q[2];
rz(-1.3126866) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(0.48373568) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023225697) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(1.1859878) q[0];
rz(-1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-2.5591154) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7459481) q[0];
sx q[0];
rz(-2.1714604) q[0];
sx q[0];
rz(-0.086918513) q[0];
rz(-pi) q[1];
rz(-2.2485579) q[2];
sx q[2];
rz(-1.6871916) q[2];
sx q[2];
rz(1.2410156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7354048) q[1];
sx q[1];
rz(-1.5899854) q[1];
sx q[1];
rz(-2.4270942) q[1];
x q[2];
rz(-1.5037687) q[3];
sx q[3];
rz(-1.9762632) q[3];
sx q[3];
rz(1.0073347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65486583) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(-3.0598818) q[2];
rz(-0.47406667) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(-1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24494568) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(-0.0078049302) q[0];
rz(-1.7410949) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(-1.1046462) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35690755) q[0];
sx q[0];
rz(-0.68455682) q[0];
sx q[0];
rz(2.8820769) q[0];
x q[1];
rz(-1.8137663) q[2];
sx q[2];
rz(-2.1967078) q[2];
sx q[2];
rz(1.6210131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.9756411) q[1];
sx q[1];
rz(-0.8756606) q[1];
sx q[1];
rz(2.402311) q[1];
rz(2.9748627) q[3];
sx q[3];
rz(-2.03074) q[3];
sx q[3];
rz(-2.102364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8905939) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(0.55523038) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(-1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5884488) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(-0.061766457) q[0];
rz(-2.899509) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(1.1118836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33169532) q[0];
sx q[0];
rz(-2.2391717) q[0];
sx q[0];
rz(-2.6033127) q[0];
rz(-pi) q[1];
rz(1.7467473) q[2];
sx q[2];
rz(-0.07882747) q[2];
sx q[2];
rz(1.4836756) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55918499) q[1];
sx q[1];
rz(-1.8473986) q[1];
sx q[1];
rz(2.4145967) q[1];
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
rz(-0.76433864) q[2];
sx q[2];
rz(0.81364441) q[2];
rz(-1.404473) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(-0.61541921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62548816) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(1.7161436) q[0];
rz(-1.5215993) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(-0.63751784) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0901776) q[0];
sx q[0];
rz(-0.99778236) q[0];
sx q[0];
rz(-2.2374002) q[0];
rz(-pi) q[1];
rz(-1.9279187) q[2];
sx q[2];
rz(-2.2631858) q[2];
sx q[2];
rz(1.9879607) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9962822) q[1];
sx q[1];
rz(-1.1728371) q[1];
sx q[1];
rz(2.5064962) q[1];
rz(-pi) q[2];
rz(0.57070891) q[3];
sx q[3];
rz(-1.4548886) q[3];
sx q[3];
rz(1.4012208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0104388) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(-1.1361702) q[2];
rz(-1.6561967) q[3];
sx q[3];
rz(-2.6172726) q[3];
sx q[3];
rz(-1.2095399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(1.6256325) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(2.4654454) q[0];
rz(0.82538429) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(-1.0151781) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7069106) q[0];
sx q[0];
rz(-0.43723956) q[0];
sx q[0];
rz(0.31243639) q[0];
rz(0.65281547) q[2];
sx q[2];
rz(-1.6932994) q[2];
sx q[2];
rz(-1.1268238) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.059541313) q[1];
sx q[1];
rz(-2.6350628) q[1];
sx q[1];
rz(2.3561213) q[1];
rz(-pi) q[2];
rz(2.7916662) q[3];
sx q[3];
rz(-2.3438128) q[3];
sx q[3];
rz(1.5087023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72145808) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-0.96735111) q[2];
rz(-1.597065) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(-2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5230781) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(-0.05649795) q[0];
rz(2.0691195) q[1];
sx q[1];
rz(-1.2697376) q[1];
sx q[1];
rz(-1.4046232) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79849762) q[0];
sx q[0];
rz(-1.6999177) q[0];
sx q[0];
rz(1.207418) q[0];
x q[1];
rz(-1.4490453) q[2];
sx q[2];
rz(-1.2687614) q[2];
sx q[2];
rz(1.9027325) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7787331) q[1];
sx q[1];
rz(-0.86231386) q[1];
sx q[1];
rz(0.94019903) q[1];
rz(-pi) q[2];
rz(-2.9910562) q[3];
sx q[3];
rz(-1.5794465) q[3];
sx q[3];
rz(1.0516402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29356062) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(2.8354697) q[2];
rz(0.39811578) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067779) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(-1.9819992) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(-2.878306) q[2];
sx q[2];
rz(-1.2821715) q[2];
sx q[2];
rz(0.53496219) q[2];
rz(2.8108834) q[3];
sx q[3];
rz(-1.0996795) q[3];
sx q[3];
rz(-0.75547937) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];