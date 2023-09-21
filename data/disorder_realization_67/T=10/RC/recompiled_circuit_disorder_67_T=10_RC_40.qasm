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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4817754) q[0];
sx q[0];
rz(-3.0660015) q[0];
sx q[0];
rz(0.48150058) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4906292) q[2];
sx q[2];
rz(-1.5032094) q[2];
sx q[2];
rz(-1.1196605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2343443) q[1];
sx q[1];
rz(-2.611479) q[1];
sx q[1];
rz(-2.2905473) q[1];
rz(0.30480095) q[3];
sx q[3];
rz(-1.3918575) q[3];
sx q[3];
rz(1.6182871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(-1.8908267) q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(2.1616518) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.132906) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(-2.5426478) q[0];
rz(1.3409746) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(0.96639955) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7705298) q[0];
sx q[0];
rz(-0.75936985) q[0];
sx q[0];
rz(-2.4735527) q[0];
rz(-pi) q[1];
rz(-1.5210549) q[2];
sx q[2];
rz(-0.94919862) q[2];
sx q[2];
rz(-0.13812401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.84345531) q[1];
sx q[1];
rz(-1.6846488) q[1];
sx q[1];
rz(-1.650412) q[1];
x q[2];
rz(1.1411243) q[3];
sx q[3];
rz(-2.2727192) q[3];
sx q[3];
rz(-0.3598635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.085658375) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(-2.8404964) q[2];
rz(1.9484693) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0369204) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(1.1874636) q[0];
rz(1.9056412) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(-1.8240066) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4268437) q[0];
sx q[0];
rz(-2.2509529) q[0];
sx q[0];
rz(1.8450518) q[0];
rz(0.74430978) q[2];
sx q[2];
rz(-2.8544606) q[2];
sx q[2];
rz(-2.9922275) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2059523) q[1];
sx q[1];
rz(-1.9899568) q[1];
sx q[1];
rz(1.6816891) q[1];
rz(-0.1173238) q[3];
sx q[3];
rz(-1.0563207) q[3];
sx q[3];
rz(2.8770212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.3935564) q[2];
sx q[2];
rz(2.9349566) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(-2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77465039) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(2.2553717) q[0];
rz(1.0097424) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(-1.2264235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3564295) q[0];
sx q[0];
rz(-2.4117081) q[0];
sx q[0];
rz(1.4714144) q[0];
x q[1];
rz(2.1143772) q[2];
sx q[2];
rz(-2.377844) q[2];
sx q[2];
rz(-2.5716512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9981421) q[1];
sx q[1];
rz(-1.7805903) q[1];
sx q[1];
rz(-1.7721121) q[1];
rz(1.1981443) q[3];
sx q[3];
rz(-2.028392) q[3];
sx q[3];
rz(-2.3082993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4975138) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(-2.148596) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(-1.9556048) q[0];
rz(-1.7182619) q[1];
sx q[1];
rz(-1.490373) q[1];
sx q[1];
rz(0.58247724) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5929778) q[0];
sx q[0];
rz(-2.5354404) q[0];
sx q[0];
rz(1.6968326) q[0];
rz(-pi) q[1];
rz(1.7551454) q[2];
sx q[2];
rz(-0.68612387) q[2];
sx q[2];
rz(2.9550936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.40618784) q[1];
sx q[1];
rz(-1.5899854) q[1];
sx q[1];
rz(-2.4270942) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6378239) q[3];
sx q[3];
rz(-1.1653295) q[3];
sx q[3];
rz(2.1342579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(0.081710903) q[2];
rz(0.47406667) q[3];
sx q[3];
rz(-1.877955) q[3];
sx q[3];
rz(1.4985532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.24494568) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(3.1337877) q[0];
rz(-1.4004978) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(-1.1046462) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68708006) q[0];
sx q[0];
rz(-2.2283163) q[0];
sx q[0];
rz(1.3643273) q[0];
rz(2.5014624) q[2];
sx q[2];
rz(-1.7670317) q[2];
sx q[2];
rz(-3.047608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.066597477) q[1];
sx q[1];
rz(-1.0269594) q[1];
sx q[1];
rz(-0.72504136) q[1];
x q[2];
rz(0.16672991) q[3];
sx q[3];
rz(-1.1108526) q[3];
sx q[3];
rz(1.0392287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2509987) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(-2.5863623) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(-1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5884488) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(-0.061766457) q[0];
rz(-0.24208367) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(1.1118836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0430849) q[0];
sx q[0];
rz(-0.83139172) q[0];
sx q[0];
rz(-0.99494536) q[0];
rz(-pi) q[1];
rz(-1.7467473) q[2];
sx q[2];
rz(-0.07882747) q[2];
sx q[2];
rz(1.6579171) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8916787) q[1];
sx q[1];
rz(-2.2644682) q[1];
sx q[1];
rz(-1.9338884) q[1];
rz(-1.816733) q[3];
sx q[3];
rz(-1.2517559) q[3];
sx q[3];
rz(0.46003534) q[3];
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
rz(-1.404473) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62548816) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(-1.7161436) q[0];
rz(-1.5215993) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(-0.63751784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051415074) q[0];
sx q[0];
rz(-2.1438103) q[0];
sx q[0];
rz(0.90419241) q[0];
x q[1];
rz(-2.417008) q[2];
sx q[2];
rz(-1.8432334) q[2];
sx q[2];
rz(0.1833293) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4379165) q[1];
sx q[1];
rz(-2.1494467) q[1];
sx q[1];
rz(2.0520567) q[1];
x q[2];
rz(-0.2122722) q[3];
sx q[3];
rz(-2.5605154) q[3];
sx q[3];
rz(-3.1329336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1311538) q[2];
sx q[2];
rz(-0.70019478) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.5159601) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(-2.4654454) q[0];
rz(-2.3162084) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(2.1264145) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43468201) q[0];
sx q[0];
rz(-2.7043531) q[0];
sx q[0];
rz(0.31243639) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4887772) q[2];
sx q[2];
rz(-1.6932994) q[2];
sx q[2];
rz(2.0147689) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3489192) q[1];
sx q[1];
rz(-1.9209849) q[1];
sx q[1];
rz(-2.7677571) q[1];
rz(-0.76652758) q[3];
sx q[3];
rz(-1.818728) q[3];
sx q[3];
rz(2.8299696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-0.96735111) q[2];
rz(-1.597065) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5230781) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(-3.0850947) q[0];
rz(-1.0724732) q[1];
sx q[1];
rz(-1.2697376) q[1];
sx q[1];
rz(1.7369695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.343095) q[0];
sx q[0];
rz(-1.4416749) q[0];
sx q[0];
rz(1.207418) q[0];
rz(-pi) q[1];
rz(2.769906) q[2];
sx q[2];
rz(-2.8166397) q[2];
sx q[2];
rz(-1.5124958) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2056634) q[1];
sx q[1];
rz(-2.231039) q[1];
sx q[1];
rz(2.5388989) q[1];
rz(1.5795454) q[3];
sx q[3];
rz(-1.4202655) q[3];
sx q[3];
rz(-0.51784408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.848032) q[2];
sx q[2];
rz(-2.3876987) q[2];
sx q[2];
rz(-0.30612293) q[2];
rz(-2.7434769) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(-2.117363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8067779) q[0];
sx q[0];
rz(-2.0756742) q[0];
sx q[0];
rz(-2.6609127) q[0];
rz(1.1595935) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(0.26328662) q[2];
sx q[2];
rz(-1.2821715) q[2];
sx q[2];
rz(0.53496219) q[2];
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
