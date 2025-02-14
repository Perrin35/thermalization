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
rz(-0.40989447) q[0];
sx q[0];
rz(4.8978187) q[0];
sx q[0];
rz(10.495821) q[0];
rz(1.8154124) q[1];
sx q[1];
rz(-0.92547995) q[1];
sx q[1];
rz(-0.34223908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15187787) q[0];
sx q[0];
rz(-0.79024716) q[0];
sx q[0];
rz(1.0343767) q[0];
x q[1];
rz(-2.5985122) q[2];
sx q[2];
rz(-0.8870856) q[2];
sx q[2];
rz(1.8250842) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2633873) q[1];
sx q[1];
rz(-2.6942746) q[1];
sx q[1];
rz(1.8187128) q[1];
x q[2];
rz(-1.5299523) q[3];
sx q[3];
rz(-0.39642912) q[3];
sx q[3];
rz(1.2547697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5388415) q[2];
sx q[2];
rz(-2.4138236) q[2];
sx q[2];
rz(1.7068498) q[2];
rz(2.1753066) q[3];
sx q[3];
rz(-1.1937701) q[3];
sx q[3];
rz(-2.5626591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(3.0173138) q[0];
sx q[0];
rz(-0.74324981) q[0];
sx q[0];
rz(-3.0857575) q[0];
rz(2.3827379) q[1];
sx q[1];
rz(-1.1412303) q[1];
sx q[1];
rz(0.32078823) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7082328) q[0];
sx q[0];
rz(-2.367041) q[0];
sx q[0];
rz(-0.21277986) q[0];
x q[1];
rz(-0.35050138) q[2];
sx q[2];
rz(-0.81074088) q[2];
sx q[2];
rz(-1.8271918) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5949147) q[1];
sx q[1];
rz(-1.3370829) q[1];
sx q[1];
rz(-3.1114034) q[1];
rz(-pi) q[2];
rz(-1.1029712) q[3];
sx q[3];
rz(-0.89455869) q[3];
sx q[3];
rz(-0.53550628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32030973) q[2];
sx q[2];
rz(-1.2592969) q[2];
sx q[2];
rz(2.8584976) q[2];
rz(-0.29575944) q[3];
sx q[3];
rz(-1.7167973) q[3];
sx q[3];
rz(-1.504771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54842424) q[0];
sx q[0];
rz(-1.9110837) q[0];
sx q[0];
rz(3.044627) q[0];
rz(-1.7123669) q[1];
sx q[1];
rz(-1.2173419) q[1];
sx q[1];
rz(-0.355535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9622346) q[0];
sx q[0];
rz(-1.0197687) q[0];
sx q[0];
rz(-1.9105115) q[0];
x q[1];
rz(-2.9121551) q[2];
sx q[2];
rz(-0.4248473) q[2];
sx q[2];
rz(-2.4160624) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3073193) q[1];
sx q[1];
rz(-1.6641698) q[1];
sx q[1];
rz(1.1591256) q[1];
rz(1.612298) q[3];
sx q[3];
rz(-1.8222162) q[3];
sx q[3];
rz(1.6715028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43704438) q[2];
sx q[2];
rz(-1.1532249) q[2];
sx q[2];
rz(2.8538749) q[2];
rz(-2.8904166) q[3];
sx q[3];
rz(-2.0468678) q[3];
sx q[3];
rz(-1.2742961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1890202) q[0];
sx q[0];
rz(-0.59307161) q[0];
sx q[0];
rz(0.75743842) q[0];
rz(-2.3144552) q[1];
sx q[1];
rz(-2.4816781) q[1];
sx q[1];
rz(-1.3077024) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2986733) q[0];
sx q[0];
rz(-1.5184776) q[0];
sx q[0];
rz(-1.6949881) q[0];
rz(-pi) q[1];
x q[1];
rz(2.04129) q[2];
sx q[2];
rz(-1.845282) q[2];
sx q[2];
rz(0.16484552) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4003936) q[1];
sx q[1];
rz(-2.315042) q[1];
sx q[1];
rz(-1.2024906) q[1];
rz(3.0809866) q[3];
sx q[3];
rz(-1.8148141) q[3];
sx q[3];
rz(-0.054494245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7149675) q[2];
sx q[2];
rz(-1.3447309) q[2];
sx q[2];
rz(1.7460543) q[2];
rz(1.7914145) q[3];
sx q[3];
rz(-1.6387286) q[3];
sx q[3];
rz(3.1135528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79513079) q[0];
sx q[0];
rz(-2.6030354) q[0];
sx q[0];
rz(-1.441347) q[0];
rz(0.92789188) q[1];
sx q[1];
rz(-2.3065232) q[1];
sx q[1];
rz(-2.3477614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5466325) q[0];
sx q[0];
rz(-1.4671456) q[0];
sx q[0];
rz(-1.9846546) q[0];
rz(-pi) q[1];
rz(1.4537816) q[2];
sx q[2];
rz(-1.7378561) q[2];
sx q[2];
rz(2.2599205) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1542357) q[1];
sx q[1];
rz(-1.8786228) q[1];
sx q[1];
rz(2.1218142) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1855741) q[3];
sx q[3];
rz(-1.0971991) q[3];
sx q[3];
rz(-1.4095167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64357197) q[2];
sx q[2];
rz(-1.4683495) q[2];
sx q[2];
rz(0.35430655) q[2];
rz(-2.0453359) q[3];
sx q[3];
rz(-2.8772964) q[3];
sx q[3];
rz(-1.8335584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3827913) q[0];
sx q[0];
rz(-0.76171869) q[0];
sx q[0];
rz(0.85269165) q[0];
rz(2.8755152) q[1];
sx q[1];
rz(-1.0238901) q[1];
sx q[1];
rz(-1.783225) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13183188) q[0];
sx q[0];
rz(-2.1340098) q[0];
sx q[0];
rz(-2.8866803) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89766111) q[2];
sx q[2];
rz(-1.3579457) q[2];
sx q[2];
rz(-2.0124679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1193124) q[1];
sx q[1];
rz(-1.90442) q[1];
sx q[1];
rz(-2.5319478) q[1];
x q[2];
rz(0.90161277) q[3];
sx q[3];
rz(-1.3517778) q[3];
sx q[3];
rz(-2.7215093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69409662) q[2];
sx q[2];
rz(-1.9393549) q[2];
sx q[2];
rz(-2.6608432) q[2];
rz(0.9350183) q[3];
sx q[3];
rz(-0.98116773) q[3];
sx q[3];
rz(-2.7267314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9347436) q[0];
sx q[0];
rz(-0.54731363) q[0];
sx q[0];
rz(-2.4547087) q[0];
rz(-2.1737449) q[1];
sx q[1];
rz(-1.5279488) q[1];
sx q[1];
rz(-2.9639249) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2813276) q[0];
sx q[0];
rz(-1.2277368) q[0];
sx q[0];
rz(-2.706091) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2306941) q[2];
sx q[2];
rz(-1.0086806) q[2];
sx q[2];
rz(2.0695994) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8502151) q[1];
sx q[1];
rz(-1.134877) q[1];
sx q[1];
rz(0.30471431) q[1];
rz(-0.068818109) q[3];
sx q[3];
rz(-1.4806169) q[3];
sx q[3];
rz(-0.18776151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0873506) q[2];
sx q[2];
rz(-2.6998616) q[2];
sx q[2];
rz(2.1596215) q[2];
rz(2.8864077) q[3];
sx q[3];
rz(-1.1533777) q[3];
sx q[3];
rz(-1.0758146) q[3];
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
rz(2.0725237) q[0];
sx q[0];
rz(-2.531932) q[0];
sx q[0];
rz(-0.12741086) q[0];
rz(-1.472507) q[1];
sx q[1];
rz(-1.9576879) q[1];
sx q[1];
rz(1.6944108) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1217598) q[0];
sx q[0];
rz(-2.7708672) q[0];
sx q[0];
rz(2.6080568) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5980392) q[2];
sx q[2];
rz(-1.8951891) q[2];
sx q[2];
rz(-2.3860562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9738335) q[1];
sx q[1];
rz(-2.9855507) q[1];
sx q[1];
rz(2.6521801) q[1];
x q[2];
rz(-0.98663892) q[3];
sx q[3];
rz(-1.3615063) q[3];
sx q[3];
rz(1.926146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.89333263) q[2];
sx q[2];
rz(-0.97951639) q[2];
sx q[2];
rz(2.7323006) q[2];
rz(-0.69581699) q[3];
sx q[3];
rz(-0.69605637) q[3];
sx q[3];
rz(-0.61923641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5291587) q[0];
sx q[0];
rz(-2.9127064) q[0];
sx q[0];
rz(-2.9869475) q[0];
rz(1.0908499) q[1];
sx q[1];
rz(-2.2583074) q[1];
sx q[1];
rz(-1.8152016) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94448254) q[0];
sx q[0];
rz(-1.8354699) q[0];
sx q[0];
rz(0.84828429) q[0];
rz(-2.0699851) q[2];
sx q[2];
rz(-1.6130924) q[2];
sx q[2];
rz(-2.7214839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5828344) q[1];
sx q[1];
rz(-1.5213005) q[1];
sx q[1];
rz(-0.82079408) q[1];
x q[2];
rz(2.6821413) q[3];
sx q[3];
rz(-1.748845) q[3];
sx q[3];
rz(2.0083754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.078330127) q[2];
sx q[2];
rz(-0.41389725) q[2];
sx q[2];
rz(0.70270339) q[2];
rz(3.0472158) q[3];
sx q[3];
rz(-2.1636212) q[3];
sx q[3];
rz(-0.92488658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5265882) q[0];
sx q[0];
rz(-2.1115392) q[0];
sx q[0];
rz(-2.7870542) q[0];
rz(1.4226557) q[1];
sx q[1];
rz(-1.645003) q[1];
sx q[1];
rz(0.52182251) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0063534) q[0];
sx q[0];
rz(-1.420605) q[0];
sx q[0];
rz(-1.9271668) q[0];
rz(1.8759236) q[2];
sx q[2];
rz(-1.702945) q[2];
sx q[2];
rz(-0.73763359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.21510151) q[1];
sx q[1];
rz(-1.0287849) q[1];
sx q[1];
rz(2.3597329) q[1];
rz(0.75504889) q[3];
sx q[3];
rz(-0.88249373) q[3];
sx q[3];
rz(-1.2699708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.78174138) q[2];
sx q[2];
rz(-1.968911) q[2];
sx q[2];
rz(2.4594128) q[2];
rz(-1.3683246) q[3];
sx q[3];
rz(-1.6006288) q[3];
sx q[3];
rz(2.1413596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.901004) q[0];
sx q[0];
rz(-2.1841342) q[0];
sx q[0];
rz(2.8484455) q[0];
rz(2.3121569) q[1];
sx q[1];
rz(-1.71143) q[1];
sx q[1];
rz(-1.5477187) q[1];
rz(-2.2915056) q[2];
sx q[2];
rz(-2.8003393) q[2];
sx q[2];
rz(0.48366) q[2];
rz(-2.8021723) q[3];
sx q[3];
rz(-0.11133736) q[3];
sx q[3];
rz(2.8471024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
