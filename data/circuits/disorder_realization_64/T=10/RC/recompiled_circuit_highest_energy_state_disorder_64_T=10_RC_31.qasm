OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0567226) q[0];
sx q[0];
rz(-1.280008) q[0];
sx q[0];
rz(-0.7260538) q[0];
rz(0.59745204) q[1];
sx q[1];
rz(-2.6260881) q[1];
sx q[1];
rz(-2.1093624) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7643578) q[0];
sx q[0];
rz(-2.4774158) q[0];
sx q[0];
rz(2.9773877) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1025508) q[2];
sx q[2];
rz(-1.0532951) q[2];
sx q[2];
rz(-1.3618748) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1223678) q[1];
sx q[1];
rz(-1.6926973) q[1];
sx q[1];
rz(-2.8278964) q[1];
rz(-pi) q[2];
rz(-0.80268245) q[3];
sx q[3];
rz(-1.2000053) q[3];
sx q[3];
rz(-0.30457531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7008179) q[2];
sx q[2];
rz(-2.5371607) q[2];
sx q[2];
rz(-3.0536998) q[2];
rz(-1.4676189) q[3];
sx q[3];
rz(-1.847495) q[3];
sx q[3];
rz(-0.99883336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0198332) q[0];
sx q[0];
rz(-1.3061981) q[0];
sx q[0];
rz(1.649296) q[0];
rz(0.78798931) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(-0.96510395) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.777634) q[0];
sx q[0];
rz(-2.2978362) q[0];
sx q[0];
rz(-1.6939837) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1568753) q[2];
sx q[2];
rz(-2.3986536) q[2];
sx q[2];
rz(2.3791831) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.31307855) q[1];
sx q[1];
rz(-2.0277436) q[1];
sx q[1];
rz(-2.9345153) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4697462) q[3];
sx q[3];
rz(-1.489991) q[3];
sx q[3];
rz(1.4572136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1284156) q[2];
sx q[2];
rz(-2.8030277) q[2];
sx q[2];
rz(1.817912) q[2];
rz(-2.2549818) q[3];
sx q[3];
rz(-1.1611791) q[3];
sx q[3];
rz(0.19865856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897268) q[0];
sx q[0];
rz(-0.97049814) q[0];
sx q[0];
rz(-0.6955198) q[0];
rz(-0.8356525) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(-2.4998891) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27370975) q[0];
sx q[0];
rz(-0.84016227) q[0];
sx q[0];
rz(0.19489906) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19410816) q[2];
sx q[2];
rz(-2.420937) q[2];
sx q[2];
rz(-2.1168762) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5776331) q[1];
sx q[1];
rz(-2.0490987) q[1];
sx q[1];
rz(0.35525124) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4536127) q[3];
sx q[3];
rz(-2.0350581) q[3];
sx q[3];
rz(-1.9339428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.16331638) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(-1.8070492) q[2];
rz(-1.7415907) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(1.1388206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1432994) q[0];
sx q[0];
rz(-0.87109733) q[0];
sx q[0];
rz(0.77348462) q[0];
rz(3.0553014) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(2.6533244) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21742709) q[0];
sx q[0];
rz(-2.0323951) q[0];
sx q[0];
rz(-1.6469) q[0];
x q[1];
rz(-2.4441798) q[2];
sx q[2];
rz(-1.1367172) q[2];
sx q[2];
rz(3.0212439) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3662373) q[1];
sx q[1];
rz(-1.7960494) q[1];
sx q[1];
rz(0.83323343) q[1];
x q[2];
rz(-2.7375324) q[3];
sx q[3];
rz(-1.1935788) q[3];
sx q[3];
rz(0.67507832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.46828541) q[2];
sx q[2];
rz(-2.0279334) q[2];
sx q[2];
rz(2.9735273) q[2];
rz(-0.42260653) q[3];
sx q[3];
rz(-0.97413078) q[3];
sx q[3];
rz(0.15338038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76822686) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(-1.9510829) q[0];
rz(0.52781421) q[1];
sx q[1];
rz(-2.0595136) q[1];
sx q[1];
rz(-0.0091008069) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6884228) q[0];
sx q[0];
rz(-0.97081447) q[0];
sx q[0];
rz(1.021941) q[0];
rz(-pi) q[1];
rz(1.6318237) q[2];
sx q[2];
rz(-1.5102855) q[2];
sx q[2];
rz(0.64693816) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.121351) q[1];
sx q[1];
rz(-0.96881908) q[1];
sx q[1];
rz(-0.51687981) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9420932) q[3];
sx q[3];
rz(-1.1807943) q[3];
sx q[3];
rz(2.3671234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8268383) q[2];
sx q[2];
rz(-2.645731) q[2];
sx q[2];
rz(2.4763988) q[2];
rz(-2.1151755) q[3];
sx q[3];
rz(-2.0668991) q[3];
sx q[3];
rz(-0.78072602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6219532) q[0];
sx q[0];
rz(-2.6060947) q[0];
sx q[0];
rz(-0.39805472) q[0];
rz(-1.4705426) q[1];
sx q[1];
rz(-0.87659756) q[1];
sx q[1];
rz(-2.3614531) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0288643) q[0];
sx q[0];
rz(-2.5327589) q[0];
sx q[0];
rz(0.91797773) q[0];
rz(-0.69177912) q[2];
sx q[2];
rz(-0.31111141) q[2];
sx q[2];
rz(-1.9115314) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4698339) q[1];
sx q[1];
rz(-2.2398754) q[1];
sx q[1];
rz(1.9356739) q[1];
rz(-2.2557115) q[3];
sx q[3];
rz(-1.8198593) q[3];
sx q[3];
rz(-1.1108171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2842497) q[2];
sx q[2];
rz(-2.0390859) q[2];
sx q[2];
rz(-0.16172376) q[2];
rz(3.0297847) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(-2.4017754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.1424471) q[0];
sx q[0];
rz(-1.7113926) q[0];
sx q[0];
rz(0.69851843) q[0];
rz(-1.0922208) q[1];
sx q[1];
rz(-2.3869546) q[1];
sx q[1];
rz(-0.05365595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2630477) q[0];
sx q[0];
rz(-1.3113271) q[0];
sx q[0];
rz(2.3058913) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8648134) q[2];
sx q[2];
rz(-0.55241441) q[2];
sx q[2];
rz(-0.68589003) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3813579) q[1];
sx q[1];
rz(-1.5173727) q[1];
sx q[1];
rz(-2.670856) q[1];
rz(-1.7683623) q[3];
sx q[3];
rz(-0.65614349) q[3];
sx q[3];
rz(0.11696091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0003537) q[2];
sx q[2];
rz(-0.69811368) q[2];
sx q[2];
rz(0.21256438) q[2];
rz(2.8138748) q[3];
sx q[3];
rz(-2.1543584) q[3];
sx q[3];
rz(-1.3535961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.245529) q[0];
sx q[0];
rz(-2.0241757) q[0];
sx q[0];
rz(-2.8118706) q[0];
rz(2.3700736) q[1];
sx q[1];
rz(-2.3604269) q[1];
sx q[1];
rz(-2.3158997) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86733183) q[0];
sx q[0];
rz(-1.6578339) q[0];
sx q[0];
rz(2.4768193) q[0];
rz(-pi) q[1];
rz(0.97445455) q[2];
sx q[2];
rz(-1.6357517) q[2];
sx q[2];
rz(-0.37298733) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.030071478) q[1];
sx q[1];
rz(-0.81395212) q[1];
sx q[1];
rz(-1.4210644) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8844374) q[3];
sx q[3];
rz(-1.7622593) q[3];
sx q[3];
rz(1.934727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.19994911) q[2];
sx q[2];
rz(-1.4120833) q[2];
sx q[2];
rz(-1.4818209) q[2];
rz(-3.0485349) q[3];
sx q[3];
rz(-1.8424282) q[3];
sx q[3];
rz(0.53946462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900443) q[0];
sx q[0];
rz(-2.7261782) q[0];
sx q[0];
rz(-0.39733091) q[0];
rz(-0.98995248) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(-1.0362524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62972126) q[0];
sx q[0];
rz(-2.7366834) q[0];
sx q[0];
rz(0.029442336) q[0];
x q[1];
rz(-1.431688) q[2];
sx q[2];
rz(-1.221162) q[2];
sx q[2];
rz(1.0316499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3237988) q[1];
sx q[1];
rz(-1.7897871) q[1];
sx q[1];
rz(1.6621672) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0898706) q[3];
sx q[3];
rz(-1.241893) q[3];
sx q[3];
rz(0.84612209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.022126023) q[2];
sx q[2];
rz(-0.72212044) q[2];
sx q[2];
rz(1.3163346) q[2];
rz(-2.8890166) q[3];
sx q[3];
rz(-1.794869) q[3];
sx q[3];
rz(2.778497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0040358) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(-2.9055415) q[0];
rz(-3.0421742) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(1.258446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1374219) q[0];
sx q[0];
rz(-1.6485456) q[0];
sx q[0];
rz(-0.56044062) q[0];
rz(0.61014097) q[2];
sx q[2];
rz(-2.379866) q[2];
sx q[2];
rz(1.6118647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8736781) q[1];
sx q[1];
rz(-2.2482052) q[1];
sx q[1];
rz(-1.4235086) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80983617) q[3];
sx q[3];
rz(-2.5812529) q[3];
sx q[3];
rz(0.59691012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9749757) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(-2.100259) q[2];
rz(2.5552022) q[3];
sx q[3];
rz(-2.6643463) q[3];
sx q[3];
rz(-2.3311116) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7753684) q[0];
sx q[0];
rz(-1.0564221) q[0];
sx q[0];
rz(-1.139241) q[0];
rz(-2.1495023) q[1];
sx q[1];
rz(-1.6012259) q[1];
sx q[1];
rz(1.59902) q[1];
rz(-0.31789657) q[2];
sx q[2];
rz(-1.7902725) q[2];
sx q[2];
rz(-0.27223311) q[2];
rz(0.59320368) q[3];
sx q[3];
rz(-1.2013669) q[3];
sx q[3];
rz(-0.47040924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
