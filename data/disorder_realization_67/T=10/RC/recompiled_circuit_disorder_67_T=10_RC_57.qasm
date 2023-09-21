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
rz(-1.5001186) q[1];
sx q[1];
rz(-2.1067696) q[1];
sx q[1];
rz(-2.1980481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9991) q[0];
sx q[0];
rz(-1.6377791) q[0];
sx q[0];
rz(-1.5357369) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4906292) q[2];
sx q[2];
rz(-1.5032094) q[2];
sx q[2];
rz(1.1196605) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.01602068) q[1];
sx q[1];
rz(-1.9106094) q[1];
sx q[1];
rz(1.1556975) q[1];
rz(0.30480095) q[3];
sx q[3];
rz(-1.7497352) q[3];
sx q[3];
rz(1.5233056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(1.8908267) q[2];
rz(-1.4261774) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(2.1616518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.0086867) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(-2.5426478) q[0];
rz(1.8006181) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(-0.96639955) q[1];
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
x q[1];
rz(0.62218372) q[2];
sx q[2];
rz(-1.6112279) q[2];
sx q[2];
rz(-1.7379023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4557138) q[1];
sx q[1];
rz(-3.0027632) q[1];
sx q[1];
rz(0.60771897) q[1];
x q[2];
rz(2.0004683) q[3];
sx q[3];
rz(-2.2727192) q[3];
sx q[3];
rz(-2.7817291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.085658375) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(-2.8404964) q[2];
rz(1.1931233) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(1.1863856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10467228) q[0];
sx q[0];
rz(-1.5043229) q[0];
sx q[0];
rz(-1.954129) q[0];
rz(1.2359515) q[1];
sx q[1];
rz(-1.0373479) q[1];
sx q[1];
rz(-1.8240066) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31908195) q[0];
sx q[0];
rz(-1.7829478) q[0];
sx q[0];
rz(-0.69885079) q[0];
rz(1.3733528) q[2];
sx q[2];
rz(-1.3609481) q[2];
sx q[2];
rz(-0.91453493) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8220362) q[1];
sx q[1];
rz(-1.6720547) q[1];
sx q[1];
rz(0.42145573) q[1];
x q[2];
rz(-2.0882323) q[3];
sx q[3];
rz(-1.6728757) q[3];
sx q[3];
rz(-1.7774338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0111771) q[2];
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
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77465039) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(-2.2553717) q[0];
rz(1.0097424) q[1];
sx q[1];
rz(-0.90615288) q[1];
sx q[1];
rz(-1.9151691) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521586) q[0];
sx q[0];
rz(-0.84531784) q[0];
sx q[0];
rz(-0.088539601) q[0];
rz(-0.45986508) q[2];
sx q[2];
rz(-0.93732873) q[2];
sx q[2];
rz(1.8749274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.46982161) q[1];
sx q[1];
rz(-1.7676395) q[1];
sx q[1];
rz(-0.21398869) q[1];
x q[2];
rz(-2.6552116) q[3];
sx q[3];
rz(-1.2380621) q[3];
sx q[3];
rz(-2.2330724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6440789) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(0.99299661) q[2];
rz(1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(1.1859878) q[0];
rz(1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-0.58247724) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5486149) q[0];
sx q[0];
rz(-2.5354404) q[0];
sx q[0];
rz(1.4447601) q[0];
x q[1];
rz(1.3864473) q[2];
sx q[2];
rz(-0.68612387) q[2];
sx q[2];
rz(-2.9550936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7354048) q[1];
sx q[1];
rz(-1.5516073) q[1];
sx q[1];
rz(0.71449844) q[1];
rz(-1.6378239) q[3];
sx q[3];
rz(-1.9762632) q[3];
sx q[3];
rz(2.1342579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4867268) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(3.0598818) q[2];
rz(-2.667526) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(1.6430395) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24494568) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(0.0078049302) q[0];
rz(-1.7410949) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(1.1046462) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4545126) q[0];
sx q[0];
rz(-0.91327635) q[0];
sx q[0];
rz(1.7772654) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3278264) q[2];
sx q[2];
rz(-2.1967078) q[2];
sx q[2];
rz(1.6210131) q[2];
x q[3];
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
rz(-2.9748627) q[3];
sx q[3];
rz(-1.1108526) q[3];
sx q[3];
rz(-2.102364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2509987) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(0.55523038) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(1.593332) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5884488) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(0.061766457) q[0];
rz(-2.899509) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(-1.1118836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0985078) q[0];
sx q[0];
rz(-2.3102009) q[0];
sx q[0];
rz(-2.1466473) q[0];
rz(-pi) q[1];
rz(1.7467473) q[2];
sx q[2];
rz(-3.0627652) q[2];
sx q[2];
rz(-1.4836756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.249914) q[1];
sx q[1];
rz(-0.87712446) q[1];
sx q[1];
rz(-1.9338884) q[1];
x q[2];
rz(1.3248596) q[3];
sx q[3];
rz(-1.8898367) q[3];
sx q[3];
rz(-0.46003534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3322488) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(2.3279482) q[2];
rz(1.404473) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(-2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.5161045) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(1.425449) q[0];
rz(1.6199934) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(-0.63751784) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051415074) q[0];
sx q[0];
rz(-0.99778236) q[0];
sx q[0];
rz(-0.90419241) q[0];
x q[1];
rz(-2.417008) q[2];
sx q[2];
rz(-1.2983592) q[2];
sx q[2];
rz(2.9582634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70367614) q[1];
sx q[1];
rz(-0.99214593) q[1];
sx q[1];
rz(-1.089536) q[1];
rz(-2.5708837) q[3];
sx q[3];
rz(-1.686704) q[3];
sx q[3];
rz(1.7403719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1311538) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(2.0054224) q[2];
rz(-1.4853959) q[3];
sx q[3];
rz(-2.6172726) q[3];
sx q[3];
rz(1.2095399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.6256325) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(-2.4654454) q[0];
rz(-2.3162084) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(-1.0151781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43468201) q[0];
sx q[0];
rz(-2.7043531) q[0];
sx q[0];
rz(-2.8291563) q[0];
rz(-2.9416111) q[2];
sx q[2];
rz(-2.4790384) q[2];
sx q[2];
rz(-2.8560864) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3489192) q[1];
sx q[1];
rz(-1.2206077) q[1];
sx q[1];
rz(-2.7677571) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3750651) q[3];
sx q[3];
rz(-1.3228647) q[3];
sx q[3];
rz(2.8299696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72145808) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-2.1742415) q[2];
rz(-1.597065) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(-0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5230781) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(0.05649795) q[0];
rz(-1.0724732) q[1];
sx q[1];
rz(-1.2697376) q[1];
sx q[1];
rz(-1.4046232) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79849762) q[0];
sx q[0];
rz(-1.6999177) q[0];
sx q[0];
rz(-1.207418) q[0];
x q[1];
rz(-0.3041515) q[2];
sx q[2];
rz(-1.4545822) q[2];
sx q[2];
rz(-2.846037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7787331) q[1];
sx q[1];
rz(-2.2792788) q[1];
sx q[1];
rz(-0.94019903) q[1];
rz(-pi) q[2];
rz(1.5620473) q[3];
sx q[3];
rz(-1.4202655) q[3];
sx q[3];
rz(0.51784408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29356062) q[2];
sx q[2];
rz(-2.3876987) q[2];
sx q[2];
rz(2.8354697) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(-1.0242296) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(1.2724614) q[2];
sx q[2];
rz(-1.3186426) q[2];
sx q[2];
rz(-1.1124055) q[2];
rz(0.33070926) q[3];
sx q[3];
rz(-2.0419131) q[3];
sx q[3];
rz(2.3861133) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
