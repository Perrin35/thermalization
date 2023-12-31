OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7473937) q[0];
sx q[0];
rz(-2.6497901) q[0];
sx q[0];
rz(2.9536182) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(-1.6245276) q[1];
sx q[1];
rz(-2.7741073) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6715235) q[0];
sx q[0];
rz(-2.1052261) q[0];
sx q[0];
rz(3.021391) q[0];
rz(2.0619225) q[2];
sx q[2];
rz(-1.6899781) q[2];
sx q[2];
rz(1.0246547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5012813) q[1];
sx q[1];
rz(-0.2959364) q[1];
sx q[1];
rz(2.179115) q[1];
rz(-pi) q[2];
rz(2.6731657) q[3];
sx q[3];
rz(-2.763063) q[3];
sx q[3];
rz(-2.1551876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.964103) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(0.5509848) q[2];
rz(-1.3059113) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630163) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(0.4719032) q[0];
rz(0.42981237) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(-0.93634161) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3204526) q[0];
sx q[0];
rz(-1.8006514) q[0];
sx q[0];
rz(-1.6605404) q[0];
x q[1];
rz(1.9905375) q[2];
sx q[2];
rz(-0.83101666) q[2];
sx q[2];
rz(1.3295528) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.97570005) q[1];
sx q[1];
rz(-2.4411538) q[1];
sx q[1];
rz(2.9794934) q[1];
rz(-pi) q[2];
rz(-3.1182321) q[3];
sx q[3];
rz(-2.0733895) q[3];
sx q[3];
rz(0.7487637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3669746) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(-0.42638391) q[2];
rz(1.9042227) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24580978) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(-0.93908969) q[0];
rz(0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(-2.5476707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9868601) q[0];
sx q[0];
rz(-1.2504471) q[0];
sx q[0];
rz(2.8264168) q[0];
rz(-pi) q[1];
rz(-2.2592696) q[2];
sx q[2];
rz(-1.6154628) q[2];
sx q[2];
rz(2.4633212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8201901) q[1];
sx q[1];
rz(-2.6650975) q[1];
sx q[1];
rz(1.0914151) q[1];
rz(2.5251758) q[3];
sx q[3];
rz(-2.5865002) q[3];
sx q[3];
rz(2.4544231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64017355) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(-1.4397941) q[2];
rz(-0.38763186) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3751635) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(-2.638812) q[0];
rz(-0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(0.75685135) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6605646) q[0];
sx q[0];
rz(-1.2111944) q[0];
sx q[0];
rz(-1.9268376) q[0];
rz(-pi) q[1];
rz(2.0365305) q[2];
sx q[2];
rz(-1.5591991) q[2];
sx q[2];
rz(-2.2955745) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6824324) q[1];
sx q[1];
rz(-1.2130514) q[1];
sx q[1];
rz(2.1898502) q[1];
rz(-pi) q[2];
rz(0.73369153) q[3];
sx q[3];
rz(-1.9568223) q[3];
sx q[3];
rz(1.8611849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7148774) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(-1.654401) q[2];
rz(0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(2.5845161) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29397598) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-0.75772444) q[0];
rz(-1.853653) q[1];
sx q[1];
rz(-0.92823354) q[1];
sx q[1];
rz(-1.0505189) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3462853) q[0];
sx q[0];
rz(-1.2588358) q[0];
sx q[0];
rz(-2.8739268) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1804785) q[2];
sx q[2];
rz(-0.69176199) q[2];
sx q[2];
rz(-1.320653) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1615636) q[1];
sx q[1];
rz(-2.1254351) q[1];
sx q[1];
rz(-2.5142923) q[1];
x q[2];
rz(2.1220783) q[3];
sx q[3];
rz(-0.84421221) q[3];
sx q[3];
rz(-0.36908484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22333764) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(-2.5081432) q[2];
rz(1.9472306) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(-0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69960064) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(0.25318405) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(-1.6794499) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77308649) q[0];
sx q[0];
rz(-0.22531548) q[0];
sx q[0];
rz(0.36264514) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8842823) q[2];
sx q[2];
rz(-2.4682211) q[2];
sx q[2];
rz(1.1416669) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0837005) q[1];
sx q[1];
rz(-1.18827) q[1];
sx q[1];
rz(1.7544569) q[1];
rz(-0.31520321) q[3];
sx q[3];
rz(-1.8990371) q[3];
sx q[3];
rz(-2.4344276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.2016466) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(-0.80491006) q[2];
rz(1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(-2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.2729623) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(0.92765635) q[0];
rz(-2.1169128) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(-1.0120846) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1597848) q[0];
sx q[0];
rz(-2.4111528) q[0];
sx q[0];
rz(0.83321379) q[0];
x q[1];
rz(-0.32832844) q[2];
sx q[2];
rz(-2.7331181) q[2];
sx q[2];
rz(0.35818737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.95841366) q[1];
sx q[1];
rz(-1.603754) q[1];
sx q[1];
rz(2.5977913) q[1];
x q[2];
rz(1.6530232) q[3];
sx q[3];
rz(-1.0136908) q[3];
sx q[3];
rz(-1.132387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53081375) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(0.42993316) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(-0.53340069) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5291418) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(2.9274143) q[0];
rz(2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(-0.28373757) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8003214) q[0];
sx q[0];
rz(-0.30297908) q[0];
sx q[0];
rz(-0.11462258) q[0];
rz(-pi) q[1];
rz(-1.9010504) q[2];
sx q[2];
rz(-1.3888161) q[2];
sx q[2];
rz(0.44588003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.46854308) q[1];
sx q[1];
rz(-2.0646411) q[1];
sx q[1];
rz(1.3480575) q[1];
x q[2];
rz(-0.92054263) q[3];
sx q[3];
rz(-2.1170108) q[3];
sx q[3];
rz(1.2342681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6909137) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(1.696375) q[2];
rz(1.5971659) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(2.8022695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504869) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(-3.0723363) q[0];
rz(-1.4878558) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(1.5725296) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5208961) q[0];
sx q[0];
rz(-0.81917742) q[0];
sx q[0];
rz(-0.92185123) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6530767) q[2];
sx q[2];
rz(-1.8956208) q[2];
sx q[2];
rz(0.61330739) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20977565) q[1];
sx q[1];
rz(-1.6546384) q[1];
sx q[1];
rz(-0.2340338) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2194355) q[3];
sx q[3];
rz(-1.4033917) q[3];
sx q[3];
rz(-2.5006014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9562324) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(-0.13988477) q[2];
rz(-0.36758962) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-2.4998253) q[0];
rz(-1.9104674) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(-0.26783255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5738327) q[0];
sx q[0];
rz(-1.1288252) q[0];
sx q[0];
rz(-2.3964336) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1688813) q[2];
sx q[2];
rz(-2.9209666) q[2];
sx q[2];
rz(1.0704744) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6707582) q[1];
sx q[1];
rz(-1.8948312) q[1];
sx q[1];
rz(-1.2653989) q[1];
x q[2];
rz(1.0288826) q[3];
sx q[3];
rz(-1.8174603) q[3];
sx q[3];
rz(-2.820462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(-2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(0.5464856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(2.3616882) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(0.98942479) q[2];
sx q[2];
rz(-0.94716723) q[2];
sx q[2];
rz(-0.92826044) q[2];
rz(1.3569309) q[3];
sx q[3];
rz(-0.90448096) q[3];
sx q[3];
rz(2.7703551) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
