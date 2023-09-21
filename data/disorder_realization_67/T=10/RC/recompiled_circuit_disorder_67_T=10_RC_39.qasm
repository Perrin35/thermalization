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
rz(-1.0597205) q[0];
sx q[0];
rz(0.73097316) q[0];
rz(1.641474) q[1];
sx q[1];
rz(-1.0348231) q[1];
sx q[1];
rz(2.1980481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4306513) q[0];
sx q[0];
rz(-1.5358155) q[0];
sx q[0];
rz(0.067023858) q[0];
rz(-pi) q[1];
rz(2.4906292) q[2];
sx q[2];
rz(-1.5032094) q[2];
sx q[2];
rz(-1.1196605) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.125572) q[1];
sx q[1];
rz(-1.9106094) q[1];
sx q[1];
rz(-1.1556975) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3834076) q[3];
sx q[3];
rz(-1.2710147) q[3];
sx q[3];
rz(0.0084458394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.25508183) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(1.250766) q[2];
rz(1.7154153) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0086867) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(-0.59894484) q[0];
rz(1.3409746) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(-0.96639955) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67993977) q[0];
sx q[0];
rz(-1.1302117) q[0];
sx q[0];
rz(0.64042129) q[0];
x q[1];
rz(2.5194089) q[2];
sx q[2];
rz(-1.6112279) q[2];
sx q[2];
rz(1.7379023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71827723) q[1];
sx q[1];
rz(-1.4916972) q[1];
sx q[1];
rz(0.11421108) q[1];
x q[2];
rz(-1.1411243) q[3];
sx q[3];
rz(-2.2727192) q[3];
sx q[3];
rz(-2.7817291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.085658375) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(2.8404964) q[2];
rz(1.1931233) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.10467228) q[0];
sx q[0];
rz(-1.5043229) q[0];
sx q[0];
rz(1.1874636) q[0];
rz(-1.9056412) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(-1.3175861) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31908195) q[0];
sx q[0];
rz(-1.7829478) q[0];
sx q[0];
rz(0.69885079) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3972829) q[2];
sx q[2];
rz(-0.28713206) q[2];
sx q[2];
rz(0.14936514) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8220362) q[1];
sx q[1];
rz(-1.6720547) q[1];
sx q[1];
rz(-0.42145573) q[1];
rz(-3.0242689) q[3];
sx q[3];
rz(-2.085272) q[3];
sx q[3];
rz(2.8770212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.3935564) q[2];
sx q[2];
rz(-0.2066361) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(-0.89282435) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-2.2354398) q[1];
sx q[1];
rz(1.2264235) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14019379) q[0];
sx q[0];
rz(-1.6370019) q[0];
sx q[0];
rz(-2.2982236) q[0];
x q[1];
rz(-0.88422758) q[2];
sx q[2];
rz(-1.2049757) q[2];
sx q[2];
rz(0.58931749) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36818477) q[1];
sx q[1];
rz(-2.8518624) q[1];
sx q[1];
rz(0.75399953) q[1];
rz(-pi) q[2];
rz(-2.6552116) q[3];
sx q[3];
rz(-1.9035305) q[3];
sx q[3];
rz(2.2330724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4975138) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(2.148596) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-2.1195181) q[3];
sx q[3];
rz(0.48373568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(1.1859878) q[0];
rz(1.7182619) q[1];
sx q[1];
rz(-1.490373) q[1];
sx q[1];
rz(2.5591154) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7459481) q[0];
sx q[0];
rz(-2.1714604) q[0];
sx q[0];
rz(0.086918513) q[0];
x q[1];
rz(-1.3864473) q[2];
sx q[2];
rz(-2.4554688) q[2];
sx q[2];
rz(0.18649907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.40618784) q[1];
sx q[1];
rz(-1.5899854) q[1];
sx q[1];
rz(2.4270942) q[1];
rz(-pi) q[2];
rz(0.4062823) q[3];
sx q[3];
rz(-1.632382) q[3];
sx q[3];
rz(-2.6046034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(0.081710903) q[2];
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
rz(-pi) q[3];
sx q[3];
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
rz(-2.896647) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(3.1337877) q[0];
rz(1.7410949) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(2.0369464) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7846851) q[0];
sx q[0];
rz(-2.4570358) q[0];
sx q[0];
rz(0.25951578) q[0];
rz(-pi) q[1];
rz(-1.3278264) q[2];
sx q[2];
rz(-0.94488482) q[2];
sx q[2];
rz(1.6210131) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9340583) q[1];
sx q[1];
rz(-2.174252) q[1];
sx q[1];
rz(-0.8912837) q[1];
rz(1.1052746) q[3];
sx q[3];
rz(-1.421531) q[3];
sx q[3];
rz(2.535459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8905939) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(-2.5863623) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5531439) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(0.061766457) q[0];
rz(0.24208367) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(-1.1118836) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0430849) q[0];
sx q[0];
rz(-2.3102009) q[0];
sx q[0];
rz(-2.1466473) q[0];
rz(1.3948453) q[2];
sx q[2];
rz(-3.0627652) q[2];
sx q[2];
rz(-1.6579171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4278533) q[1];
sx q[1];
rz(-2.3728328) q[1];
sx q[1];
rz(0.40366918) q[1];
x q[2];
rz(2.506433) q[3];
sx q[3];
rz(-0.40024647) q[3];
sx q[3];
rz(-1.1349585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8093439) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(-2.3279482) q[2];
rz(-1.7371197) q[3];
sx q[3];
rz(-2.9128894) q[3];
sx q[3];
rz(2.5261734) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5161045) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(1.7161436) q[0];
rz(-1.6199934) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(-0.63751784) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051415074) q[0];
sx q[0];
rz(-2.1438103) q[0];
sx q[0];
rz(2.2374002) q[0];
rz(-pi) q[1];
rz(-0.39890639) q[2];
sx q[2];
rz(-2.3762694) q[2];
sx q[2];
rz(-1.4590291) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0829518) q[1];
sx q[1];
rz(-2.4070027) q[1];
sx q[1];
rz(2.525108) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2122722) q[3];
sx q[3];
rz(-2.5605154) q[3];
sx q[3];
rz(-0.0086590954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0104388) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(-2.0054224) q[2];
rz(1.6561967) q[3];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6256325) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(-2.4654454) q[0];
rz(2.3162084) q[1];
sx q[1];
rz(-0.4239347) q[1];
sx q[1];
rz(2.1264145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7069106) q[0];
sx q[0];
rz(-0.43723956) q[0];
sx q[0];
rz(2.8291563) q[0];
rz(-pi) q[1];
rz(-2.4887772) q[2];
sx q[2];
rz(-1.4482933) q[2];
sx q[2];
rz(1.1268238) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2296914) q[1];
sx q[1];
rz(-1.2206612) q[1];
sx q[1];
rz(1.9446816) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76652758) q[3];
sx q[3];
rz(-1.3228647) q[3];
sx q[3];
rz(-2.8299696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72145808) q[2];
sx q[2];
rz(-2.9253503) q[2];
sx q[2];
rz(-0.96735111) q[2];
rz(1.597065) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(-0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6185146) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(0.05649795) q[0];
rz(1.0724732) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.7369695) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72337379) q[0];
sx q[0];
rz(-1.9310111) q[0];
sx q[0];
rz(3.0035613) q[0];
rz(-pi) q[1];
rz(-2.8374412) q[2];
sx q[2];
rz(-1.4545822) q[2];
sx q[2];
rz(2.846037) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7787331) q[1];
sx q[1];
rz(-2.2792788) q[1];
sx q[1];
rz(2.2013936) q[1];
x q[2];
rz(-0.057617188) q[3];
sx q[3];
rz(-2.9908097) q[3];
sx q[3];
rz(0.57612102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.848032) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(0.30612293) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(-1.0242296) q[3];
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
sx q[0];
rz(-pi) q[1];
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
rz(1.9819992) q[1];
sx q[1];
rz(-1.0212785) q[1];
sx q[1];
rz(-2.8550128) q[1];
rz(2.2904916) q[2];
sx q[2];
rz(-2.7534178) q[2];
sx q[2];
rz(-0.22321246) q[2];
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
