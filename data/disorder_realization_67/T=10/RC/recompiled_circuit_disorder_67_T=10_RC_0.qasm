OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6131634) q[0];
sx q[0];
rz(-2.0818721) q[0];
sx q[0];
rz(-0.73097316) q[0];
rz(1.641474) q[1];
sx q[1];
rz(-1.0348231) q[1];
sx q[1];
rz(2.1980481) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9991) q[0];
sx q[0];
rz(-1.5038135) q[0];
sx q[0];
rz(-1.5357369) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.655683) q[2];
sx q[2];
rz(-2.2200218) q[2];
sx q[2];
rz(-2.6390586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.125572) q[1];
sx q[1];
rz(-1.2309832) q[1];
sx q[1];
rz(-1.1556975) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5991873) q[3];
sx q[3];
rz(-0.35202682) q[3];
sx q[3];
rz(0.56234081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.25508183) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(-1.8908267) q[2];
rz(-1.7154153) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0086867) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(2.5426478) q[0];
rz(1.8006181) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(0.96639955) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37106284) q[0];
sx q[0];
rz(-0.75936985) q[0];
sx q[0];
rz(-0.66803996) q[0];
rz(-pi) q[1];
rz(1.5210549) q[2];
sx q[2];
rz(-0.94919862) q[2];
sx q[2];
rz(0.13812401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6858789) q[1];
sx q[1];
rz(-3.0027632) q[1];
sx q[1];
rz(2.5338737) q[1];
rz(-2.6838449) q[3];
sx q[3];
rz(-0.80349892) q[3];
sx q[3];
rz(0.97704923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0559343) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(2.8404964) q[2];
rz(-1.9484693) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(-1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-3.0369204) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(1.954129) q[0];
rz(-1.9056412) q[1];
sx q[1];
rz(-1.0373479) q[1];
sx q[1];
rz(-1.8240066) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7147489) q[0];
sx q[0];
rz(-2.2509529) q[0];
sx q[0];
rz(-1.8450518) q[0];
x q[1];
rz(1.3733528) q[2];
sx q[2];
rz(-1.7806446) q[2];
sx q[2];
rz(-2.2270577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31955645) q[1];
sx q[1];
rz(-1.469538) q[1];
sx q[1];
rz(-0.42145573) q[1];
rz(1.0533603) q[3];
sx q[3];
rz(-1.6728757) q[3];
sx q[3];
rz(1.3641588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1304156) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(0.2066361) q[2];
rz(0.7080428) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(-0.89282435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77465039) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(-2.2553717) q[0];
rz(-1.0097424) q[1];
sx q[1];
rz(-0.90615288) q[1];
sx q[1];
rz(1.9151691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14019379) q[0];
sx q[0];
rz(-1.6370019) q[0];
sx q[0];
rz(-0.84336908) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1143772) q[2];
sx q[2];
rz(-0.76374861) q[2];
sx q[2];
rz(-0.56994146) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36818477) q[1];
sx q[1];
rz(-2.8518624) q[1];
sx q[1];
rz(2.3875931) q[1];
rz(2.5049582) q[3];
sx q[3];
rz(-2.5599179) q[3];
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
rz(-1.3138694) q[2];
sx q[2];
rz(2.148596) q[2];
rz(1.8289061) q[3];
sx q[3];
rz(-2.1195181) q[3];
sx q[3];
rz(-0.48373568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(1.9556048) q[0];
rz(1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(2.5591154) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3956446) q[0];
sx q[0];
rz(-2.1714604) q[0];
sx q[0];
rz(3.0546741) q[0];
x q[1];
rz(-2.992606) q[2];
sx q[2];
rz(-2.2431231) q[2];
sx q[2];
rz(2.7186) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7354048) q[1];
sx q[1];
rz(-1.5516073) q[1];
sx q[1];
rz(-0.71449844) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5037687) q[3];
sx q[3];
rz(-1.1653295) q[3];
sx q[3];
rz(2.1342579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65486583) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(-0.081710903) q[2];
rz(-2.667526) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(1.6430395) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24494568) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(-0.0078049302) q[0];
rz(1.4004978) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(1.1046462) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0110328) q[0];
sx q[0];
rz(-1.7337807) q[0];
sx q[0];
rz(2.4736604) q[0];
x q[1];
rz(1.8137663) q[2];
sx q[2];
rz(-2.1967078) q[2];
sx q[2];
rz(1.5205795) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9340583) q[1];
sx q[1];
rz(-0.96734069) q[1];
sx q[1];
rz(-0.8912837) q[1];
x q[2];
rz(-1.2475345) q[3];
sx q[3];
rz(-2.6544016) q[3];
sx q[3];
rz(0.67684735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8905939) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(2.5863623) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(-1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5884488) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(-0.061766457) q[0];
rz(2.899509) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(2.0297091) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88469807) q[0];
sx q[0];
rz(-1.9848794) q[0];
sx q[0];
rz(2.3143682) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7467473) q[2];
sx q[2];
rz(-0.07882747) q[2];
sx q[2];
rz(-1.4836756) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4278533) q[1];
sx q[1];
rz(-0.76875988) q[1];
sx q[1];
rz(2.7379235) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63515969) q[3];
sx q[3];
rz(-0.40024647) q[3];
sx q[3];
rz(-1.1349585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.8093439) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(2.3279482) q[2];
rz(-1.7371197) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(-2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5161045) q[0];
sx q[0];
rz(-1.3502716) q[0];
sx q[0];
rz(-1.7161436) q[0];
rz(-1.5215993) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(-2.5040748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051415074) q[0];
sx q[0];
rz(-2.1438103) q[0];
sx q[0];
rz(-0.90419241) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9279187) q[2];
sx q[2];
rz(-0.87840688) q[2];
sx q[2];
rz(1.9879607) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0829518) q[1];
sx q[1];
rz(-0.73459) q[1];
sx q[1];
rz(-0.61648468) q[1];
rz(-0.57070891) q[3];
sx q[3];
rz(-1.4548886) q[3];
sx q[3];
rz(1.7403719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1311538) q[2];
sx q[2];
rz(-0.70019478) q[2];
sx q[2];
rz(-1.1361702) q[2];
rz(1.6561967) q[3];
sx q[3];
rz(-2.6172726) q[3];
sx q[3];
rz(-1.9320528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6256325) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(2.4654454) q[0];
rz(2.3162084) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(1.0151781) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092175352) q[0];
sx q[0];
rz(-1.1560688) q[0];
sx q[0];
rz(-1.4281005) q[0];
x q[1];
rz(-1.7245618) q[2];
sx q[2];
rz(-2.2178938) q[2];
sx q[2];
rz(0.53714067) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.059541313) q[1];
sx q[1];
rz(-2.6350628) q[1];
sx q[1];
rz(-0.7854714) q[1];
x q[2];
rz(-0.34992643) q[3];
sx q[3];
rz(-2.3438128) q[3];
sx q[3];
rz(1.5087023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4201346) q[2];
sx q[2];
rz(-2.9253503) q[2];
sx q[2];
rz(-0.96735111) q[2];
rz(1.5445276) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(-2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5230781) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(-3.0850947) q[0];
rz(2.0691195) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.4046232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.343095) q[0];
sx q[0];
rz(-1.6999177) q[0];
sx q[0];
rz(-1.9341747) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4490453) q[2];
sx q[2];
rz(-1.8728313) q[2];
sx q[2];
rz(1.2388602) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9060668) q[1];
sx q[1];
rz(-2.0350254) q[1];
sx q[1];
rz(2.3266351) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5620473) q[3];
sx q[3];
rz(-1.4202655) q[3];
sx q[3];
rz(2.6237486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.848032) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(-0.30612293) q[2];
rz(0.39811578) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(-2.117363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33481471) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(1.1595935) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(-2.2904916) q[2];
sx q[2];
rz(-0.38817482) q[2];
sx q[2];
rz(2.9183802) q[2];
rz(1.0767827) q[3];
sx q[3];
rz(-1.8643338) q[3];
sx q[3];
rz(-2.1716933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
