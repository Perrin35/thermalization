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
rz(5.2483622) q[1];
sx q[1];
rz(11.622826) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7109414) q[0];
sx q[0];
rz(-1.5358155) q[0];
sx q[0];
rz(3.0745688) q[0];
rz(-pi) q[1];
rz(3.0303454) q[2];
sx q[2];
rz(-0.6539549) q[2];
sx q[2];
rz(-2.7788869) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2343443) q[1];
sx q[1];
rz(-2.611479) q[1];
sx q[1];
rz(-0.85104533) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7581851) q[3];
sx q[3];
rz(-1.2710147) q[3];
sx q[3];
rz(0.0084458394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(-1.8908267) q[2];
rz(1.7154153) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(-0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.132906) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(2.5426478) q[0];
rz(-1.3409746) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(-0.96639955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37106284) q[0];
sx q[0];
rz(-2.3822228) q[0];
sx q[0];
rz(-0.66803996) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5210549) q[2];
sx q[2];
rz(-0.94919862) q[2];
sx q[2];
rz(-3.0034686) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84345531) q[1];
sx q[1];
rz(-1.6846488) q[1];
sx q[1];
rz(1.650412) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6838449) q[3];
sx q[3];
rz(-0.80349892) q[3];
sx q[3];
rz(-0.97704923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0559343) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(2.8404964) q[2];
rz(-1.1931233) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(-1.955207) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10467228) q[0];
sx q[0];
rz(-1.5043229) q[0];
sx q[0];
rz(1.1874636) q[0];
rz(-1.9056412) q[1];
sx q[1];
rz(-1.0373479) q[1];
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
rz(2.8185185) q[0];
x q[1];
rz(0.21388162) q[2];
sx q[2];
rz(-1.7638532) q[2];
sx q[2];
rz(-0.69790998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2059523) q[1];
sx q[1];
rz(-1.1516358) q[1];
sx q[1];
rz(-1.4599035) q[1];
rz(-pi) q[2];
rz(-2.0882323) q[3];
sx q[3];
rz(-1.6728757) q[3];
sx q[3];
rz(1.3641588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(-2.9349566) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(-0.89282435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3669423) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(-2.2553717) q[0];
rz(-1.0097424) q[1];
sx q[1];
rz(-0.90615288) q[1];
sx q[1];
rz(-1.2264235) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0013989) q[0];
sx q[0];
rz(-1.6370019) q[0];
sx q[0];
rz(-0.84336908) q[0];
rz(-pi) q[1];
rz(-2.6817276) q[2];
sx q[2];
rz(-2.2042639) q[2];
sx q[2];
rz(1.8749274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46982161) q[1];
sx q[1];
rz(-1.7676395) q[1];
sx q[1];
rz(0.21398869) q[1];
x q[2];
rz(0.63663441) q[3];
sx q[3];
rz(-2.5599179) q[3];
sx q[3];
rz(0.10891529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4975138) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(2.148596) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-2.1195181) q[3];
sx q[3];
rz(-2.657857) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023225697) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(1.9556048) q[0];
rz(-1.4233308) q[1];
sx q[1];
rz(-1.490373) q[1];
sx q[1];
rz(-0.58247724) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12594189) q[0];
sx q[0];
rz(-1.6424718) q[0];
sx q[0];
rz(2.1732251) q[0];
rz(1.7551454) q[2];
sx q[2];
rz(-0.68612387) q[2];
sx q[2];
rz(-0.18649907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1867265) q[1];
sx q[1];
rz(-2.426882) q[1];
sx q[1];
rz(3.1123118) q[1];
x q[2];
rz(0.4062823) q[3];
sx q[3];
rz(-1.632382) q[3];
sx q[3];
rz(0.53698925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65486583) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(-3.0598818) q[2];
rz(-2.667526) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(-1.4985532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24494568) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(0.0078049302) q[0];
rz(1.7410949) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(2.0369464) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68708006) q[0];
sx q[0];
rz(-2.2283163) q[0];
sx q[0];
rz(-1.3643273) q[0];
rz(-pi) q[1];
rz(2.5014624) q[2];
sx q[2];
rz(-1.7670317) q[2];
sx q[2];
rz(0.093984691) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.066597477) q[1];
sx q[1];
rz(-2.1146333) q[1];
sx q[1];
rz(0.72504136) q[1];
x q[2];
rz(0.16672991) q[3];
sx q[3];
rz(-2.03074) q[3];
sx q[3];
rz(-1.0392287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2509987) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(-2.5863623) q[2];
rz(-2.9688719) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(-1.5482607) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5531439) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(-3.0798262) q[0];
rz(-0.24208367) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(-1.1118836) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0985078) q[0];
sx q[0];
rz(-0.83139172) q[0];
sx q[0];
rz(-0.99494536) q[0];
rz(1.493181) q[2];
sx q[2];
rz(-1.5845808) q[2];
sx q[2];
rz(2.8790561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71373938) q[1];
sx q[1];
rz(-2.3728328) q[1];
sx q[1];
rz(-2.7379235) q[1];
x q[2];
rz(2.8133409) q[3];
sx q[3];
rz(-1.337507) q[3];
sx q[3];
rz(1.0321898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3322488) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(0.81364441) q[2];
rz(-1.7371197) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(-2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62548816) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(1.425449) q[0];
rz(-1.5215993) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(2.5040748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051415074) q[0];
sx q[0];
rz(-0.99778236) q[0];
sx q[0];
rz(-0.90419241) q[0];
rz(1.9279187) q[2];
sx q[2];
rz(-2.2631858) q[2];
sx q[2];
rz(-1.9879607) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.70367614) q[1];
sx q[1];
rz(-0.99214593) q[1];
sx q[1];
rz(-2.0520567) q[1];
rz(-pi) q[2];
rz(1.4333126) q[3];
sx q[3];
rz(-2.1372037) q[3];
sx q[3];
rz(0.24368225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1311538) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(-1.1361702) q[2];
rz(-1.4853959) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(1.9320528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5159601) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(-2.4654454) q[0];
rz(-2.3162084) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(2.1264145) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4207941) q[0];
sx q[0];
rz(-1.4402698) q[0];
sx q[0];
rz(-0.41850787) q[0];
rz(2.9416111) q[2];
sx q[2];
rz(-0.66255424) q[2];
sx q[2];
rz(0.28550622) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91190126) q[1];
sx q[1];
rz(-1.9209314) q[1];
sx q[1];
rz(1.9446816) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3750651) q[3];
sx q[3];
rz(-1.818728) q[3];
sx q[3];
rz(2.8299696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72145808) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-2.1742415) q[2];
rz(1.5445276) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(-0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6185146) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(3.0850947) q[0];
rz(1.0724732) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(-1.4046232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0987941) q[0];
sx q[0];
rz(-2.7569175) q[0];
sx q[0];
rz(1.2205475) q[0];
x q[1];
rz(-2.769906) q[2];
sx q[2];
rz(-2.8166397) q[2];
sx q[2];
rz(1.5124958) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2056634) q[1];
sx q[1];
rz(-2.231039) q[1];
sx q[1];
rz(2.5388989) q[1];
rz(-pi) q[2];
rz(-1.5795454) q[3];
sx q[3];
rz(-1.4202655) q[3];
sx q[3];
rz(0.51784408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29356062) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(-2.8354697) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.476215) q[3];
sx q[3];
rz(1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33481471) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(-1.1595935) q[1];
sx q[1];
rz(-1.0212785) q[1];
sx q[1];
rz(-2.8550128) q[1];
rz(2.878306) q[2];
sx q[2];
rz(-1.8594212) q[2];
sx q[2];
rz(-2.6066305) q[2];
rz(-2.0648099) q[3];
sx q[3];
rz(-1.8643338) q[3];
sx q[3];
rz(-2.1716933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
