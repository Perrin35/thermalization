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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7109414) q[0];
sx q[0];
rz(-1.5358155) q[0];
sx q[0];
rz(-0.067023858) q[0];
rz(0.65096345) q[2];
sx q[2];
rz(-1.5032094) q[2];
sx q[2];
rz(-2.0219321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4409677) q[1];
sx q[1];
rz(-1.1807627) q[1];
sx q[1];
rz(2.7729211) q[1];
x q[2];
rz(2.5991873) q[3];
sx q[3];
rz(-0.35202682) q[3];
sx q[3];
rz(2.5792518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8865108) q[2];
sx q[2];
rz(-1.3811029) q[2];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0086867) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(-2.5426478) q[0];
rz(1.3409746) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(0.96639955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67993977) q[0];
sx q[0];
rz(-1.1302117) q[0];
sx q[0];
rz(0.64042129) q[0];
rz(1.6205377) q[2];
sx q[2];
rz(-0.94919862) q[2];
sx q[2];
rz(3.0034686) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2981373) q[1];
sx q[1];
rz(-1.6846488) q[1];
sx q[1];
rz(1.4911806) q[1];
x q[2];
rz(-1.1411243) q[3];
sx q[3];
rz(-2.2727192) q[3];
sx q[3];
rz(-2.7817291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.085658375) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(0.30109626) q[2];
rz(1.1931233) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(-1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(1.2359515) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(1.8240066) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353969) q[0];
sx q[0];
rz(-0.72512308) q[0];
sx q[0];
rz(-2.8185185) q[0];
rz(-pi) q[1];
rz(2.927711) q[2];
sx q[2];
rz(-1.7638532) q[2];
sx q[2];
rz(0.69790998) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8220362) q[1];
sx q[1];
rz(-1.469538) q[1];
sx q[1];
rz(-2.7201369) q[1];
x q[2];
rz(-1.3665974) q[3];
sx q[3];
rz(-2.6150828) q[3];
sx q[3];
rz(-3.112117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0111771) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(2.9349566) q[2];
rz(2.4335499) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(2.2487683) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3669423) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(-0.88622093) q[0];
rz(2.1318502) q[1];
sx q[1];
rz(-0.90615288) q[1];
sx q[1];
rz(-1.2264235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14019379) q[0];
sx q[0];
rz(-1.6370019) q[0];
sx q[0];
rz(0.84336908) q[0];
rz(-pi) q[1];
rz(-2.6817276) q[2];
sx q[2];
rz(-0.93732873) q[2];
sx q[2];
rz(-1.8749274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.46982161) q[1];
sx q[1];
rz(-1.3739532) q[1];
sx q[1];
rz(0.21398869) q[1];
rz(-1.9434483) q[3];
sx q[3];
rz(-2.028392) q[3];
sx q[3];
rz(0.83329337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4975138) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(-0.99299661) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-0.48373568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(-1.1859878) q[0];
rz(-1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-2.5591154) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3956446) q[0];
sx q[0];
rz(-2.1714604) q[0];
sx q[0];
rz(0.086918513) q[0];
rz(2.2485579) q[2];
sx q[2];
rz(-1.6871916) q[2];
sx q[2];
rz(-1.2410156) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1479668) q[1];
sx q[1];
rz(-2.2851351) q[1];
sx q[1];
rz(-1.5453969) q[1];
rz(0.1547847) q[3];
sx q[3];
rz(-0.41066658) q[3];
sx q[3];
rz(1.1759024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.65486583) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(0.081710903) q[2];
rz(2.667526) q[3];
sx q[3];
rz(-1.877955) q[3];
sx q[3];
rz(-1.4985532) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24494568) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(0.0078049302) q[0];
rz(1.7410949) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(-1.1046462) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0110328) q[0];
sx q[0];
rz(-1.4078119) q[0];
sx q[0];
rz(-2.4736604) q[0];
rz(-pi) q[1];
rz(1.8137663) q[2];
sx q[2];
rz(-2.1967078) q[2];
sx q[2];
rz(1.5205795) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1659516) q[1];
sx q[1];
rz(-2.2659321) q[1];
sx q[1];
rz(0.73928164) q[1];
x q[2];
rz(2.036318) q[3];
sx q[3];
rz(-1.7200617) q[3];
sx q[3];
rz(2.535459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2509987) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(-0.55523038) q[2];
rz(2.9688719) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(-1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.7663019) q[1];
sx q[1];
rz(1.1118836) q[1];
rz(pi/2) q[2];
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
x q[1];
rz(3.1277666) q[2];
sx q[2];
rz(-1.4931884) q[2];
sx q[2];
rz(-1.8344049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71373938) q[1];
sx q[1];
rz(-0.76875988) q[1];
sx q[1];
rz(2.7379235) q[1];
rz(-0.32825177) q[3];
sx q[3];
rz(-1.337507) q[3];
sx q[3];
rz(1.0321898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3322488) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(2.3279482) q[2];
rz(-1.7371197) q[3];
sx q[3];
rz(-2.9128894) q[3];
sx q[3];
rz(-0.61541921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5161045) q[0];
sx q[0];
rz(-1.3502716) q[0];
sx q[0];
rz(-1.425449) q[0];
rz(1.6199934) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(-2.5040748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0901776) q[0];
sx q[0];
rz(-2.1438103) q[0];
sx q[0];
rz(0.90419241) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9279187) q[2];
sx q[2];
rz(-0.87840688) q[2];
sx q[2];
rz(-1.1536319) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.058640826) q[1];
sx q[1];
rz(-2.4070027) q[1];
sx q[1];
rz(-2.525108) q[1];
x q[2];
rz(-0.2122722) q[3];
sx q[3];
rz(-0.58107725) q[3];
sx q[3];
rz(3.1329336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1311538) q[2];
sx q[2];
rz(-0.70019478) q[2];
sx q[2];
rz(1.1361702) q[2];
rz(-1.4853959) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(1.9320528) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5159601) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(-0.67614722) q[0];
rz(2.3162084) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(-2.1264145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092175352) q[0];
sx q[0];
rz(-1.9855238) q[0];
sx q[0];
rz(-1.7134922) q[0];
rz(-2.4887772) q[2];
sx q[2];
rz(-1.4482933) q[2];
sx q[2];
rz(1.1268238) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.91190126) q[1];
sx q[1];
rz(-1.2206612) q[1];
sx q[1];
rz(-1.9446816) q[1];
rz(0.76652758) q[3];
sx q[3];
rz(-1.818728) q[3];
sx q[3];
rz(0.31162308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72145808) q[2];
sx q[2];
rz(-2.9253503) q[2];
sx q[2];
rz(-2.1742415) q[2];
rz(1.597065) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5230781) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(3.0850947) q[0];
rz(-1.0724732) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.4046232) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0427986) q[0];
sx q[0];
rz(-2.7569175) q[0];
sx q[0];
rz(1.9210451) q[0];
x q[1];
rz(-2.8374412) q[2];
sx q[2];
rz(-1.4545822) q[2];
sx q[2];
rz(-0.2955557) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9359293) q[1];
sx q[1];
rz(-2.231039) q[1];
sx q[1];
rz(2.5388989) q[1];
x q[2];
rz(-1.5795454) q[3];
sx q[3];
rz(-1.7213271) q[3];
sx q[3];
rz(2.6237486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29356062) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
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
rz(1.8691312) q[2];
sx q[2];
rz(-1.8229501) q[2];
sx q[2];
rz(2.0291871) q[2];
rz(2.0648099) q[3];
sx q[3];
rz(-1.2772588) q[3];
sx q[3];
rz(0.96989934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];