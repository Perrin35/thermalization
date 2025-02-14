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
rz(2.8994695) q[0];
sx q[0];
rz(-2.4162633) q[0];
sx q[0];
rz(2.2348833) q[0];
rz(2.6746305) q[1];
sx q[1];
rz(-0.90277201) q[1];
sx q[1];
rz(0.24388193) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7142993) q[0];
sx q[0];
rz(-0.9511982) q[0];
sx q[0];
rz(1.0788973) q[0];
rz(-pi) q[1];
rz(0.93893013) q[2];
sx q[2];
rz(-2.1515232) q[2];
sx q[2];
rz(-1.9730568) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6441124) q[1];
sx q[1];
rz(-1.6429158) q[1];
sx q[1];
rz(0.18326541) q[1];
rz(-pi) q[2];
rz(-1.7964057) q[3];
sx q[3];
rz(-0.76588956) q[3];
sx q[3];
rz(1.8647068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68524086) q[2];
sx q[2];
rz(-0.46619236) q[2];
sx q[2];
rz(2.1166128) q[2];
rz(3.0023365) q[3];
sx q[3];
rz(-1.9459629) q[3];
sx q[3];
rz(0.16417424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3716632) q[0];
sx q[0];
rz(-2.0491845) q[0];
sx q[0];
rz(-1.2051693) q[0];
rz(1.4681762) q[1];
sx q[1];
rz(-1.2638777) q[1];
sx q[1];
rz(-1.42043) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44180621) q[0];
sx q[0];
rz(-1.4107293) q[0];
sx q[0];
rz(-0.82224857) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7718138) q[2];
sx q[2];
rz(-2.8716785) q[2];
sx q[2];
rz(2.6156099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6008792) q[1];
sx q[1];
rz(-2.7329325) q[1];
sx q[1];
rz(-2.1585682) q[1];
rz(-pi) q[2];
rz(2.6254856) q[3];
sx q[3];
rz(-0.81023216) q[3];
sx q[3];
rz(2.6339171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9517407) q[2];
sx q[2];
rz(-0.18988374) q[2];
sx q[2];
rz(-0.80113775) q[2];
rz(-0.43870157) q[3];
sx q[3];
rz(-1.8935685) q[3];
sx q[3];
rz(1.1130822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37757117) q[0];
sx q[0];
rz(-2.8442597) q[0];
sx q[0];
rz(-2.0024894) q[0];
rz(-2.8746919) q[1];
sx q[1];
rz(-1.8903939) q[1];
sx q[1];
rz(-0.36531726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26317715) q[0];
sx q[0];
rz(-0.52367822) q[0];
sx q[0];
rz(-2.9429463) q[0];
rz(-0.56198012) q[2];
sx q[2];
rz(-2.1393074) q[2];
sx q[2];
rz(-1.7719444) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.461044) q[1];
sx q[1];
rz(-0.37255424) q[1];
sx q[1];
rz(-1.7303321) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5015058) q[3];
sx q[3];
rz(-0.96154204) q[3];
sx q[3];
rz(-0.14441381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.9109362) q[2];
sx q[2];
rz(-0.74030423) q[2];
sx q[2];
rz(3.018107) q[2];
rz(2.2423045) q[3];
sx q[3];
rz(-1.6920009) q[3];
sx q[3];
rz(1.2353157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7742291) q[0];
sx q[0];
rz(-1.6772567) q[0];
sx q[0];
rz(2.2520219) q[0];
rz(0.74596897) q[1];
sx q[1];
rz(-1.4624701) q[1];
sx q[1];
rz(-1.7568582) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2413413) q[0];
sx q[0];
rz(-1.5129733) q[0];
sx q[0];
rz(-3.1067294) q[0];
rz(-0.64176527) q[2];
sx q[2];
rz(-1.3744397) q[2];
sx q[2];
rz(-1.8066483) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3900323) q[1];
sx q[1];
rz(-2.1728325) q[1];
sx q[1];
rz(-1.3436418) q[1];
x q[2];
rz(-2.607658) q[3];
sx q[3];
rz(-1.5435394) q[3];
sx q[3];
rz(1.4213002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7437462) q[2];
sx q[2];
rz(-2.3036239) q[2];
sx q[2];
rz(2.1261334) q[2];
rz(3.1189392) q[3];
sx q[3];
rz(-1.7706784) q[3];
sx q[3];
rz(-3.0034351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37462336) q[0];
sx q[0];
rz(-1.8648819) q[0];
sx q[0];
rz(-0.20481566) q[0];
rz(0.21518937) q[1];
sx q[1];
rz(-0.22061017) q[1];
sx q[1];
rz(2.2664216) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143863) q[0];
sx q[0];
rz(-1.7191186) q[0];
sx q[0];
rz(1.2640087) q[0];
rz(-1.4007447) q[2];
sx q[2];
rz(-0.74234662) q[2];
sx q[2];
rz(-0.60630732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.2339911) q[1];
sx q[1];
rz(-0.37461606) q[1];
sx q[1];
rz(-2.5020155) q[1];
rz(-pi) q[2];
rz(-1.7654083) q[3];
sx q[3];
rz(-0.64009692) q[3];
sx q[3];
rz(2.3605337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.508226) q[2];
sx q[2];
rz(-0.88819155) q[2];
sx q[2];
rz(1.7249031) q[2];
rz(-2.21777) q[3];
sx q[3];
rz(-0.97771907) q[3];
sx q[3];
rz(2.0098604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0210339) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(-2.293204) q[0];
rz(1.5757164) q[1];
sx q[1];
rz(-1.3278241) q[1];
sx q[1];
rz(-2.1956992) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1740055) q[0];
sx q[0];
rz(-1.286309) q[0];
sx q[0];
rz(-0.6993642) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20393238) q[2];
sx q[2];
rz(-2.3205767) q[2];
sx q[2];
rz(-1.4216139) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.083748181) q[1];
sx q[1];
rz(-1.9942772) q[1];
sx q[1];
rz(0.64238703) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35870798) q[3];
sx q[3];
rz(-0.51878319) q[3];
sx q[3];
rz(2.0174873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1762323) q[2];
sx q[2];
rz(-1.9051899) q[2];
sx q[2];
rz(1.7410834) q[2];
rz(1.6598353) q[3];
sx q[3];
rz(-1.3503431) q[3];
sx q[3];
rz(0.19937936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7401212) q[0];
sx q[0];
rz(-0.49982247) q[0];
sx q[0];
rz(-1.9299782) q[0];
rz(0.47677332) q[1];
sx q[1];
rz(-1.2766653) q[1];
sx q[1];
rz(1.635294) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20902987) q[0];
sx q[0];
rz(-1.8537921) q[0];
sx q[0];
rz(2.1542633) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25038211) q[2];
sx q[2];
rz(-1.8844452) q[2];
sx q[2];
rz(0.55082441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8426714) q[1];
sx q[1];
rz(-1.6940882) q[1];
sx q[1];
rz(2.4409632) q[1];
rz(-pi) q[2];
rz(-0.52495082) q[3];
sx q[3];
rz(-0.54407507) q[3];
sx q[3];
rz(1.8067738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3393042) q[2];
sx q[2];
rz(-1.5172639) q[2];
sx q[2];
rz(2.9586672) q[2];
rz(0.66169345) q[3];
sx q[3];
rz(-1.2122093) q[3];
sx q[3];
rz(-0.70484149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45109192) q[0];
sx q[0];
rz(-1.4649614) q[0];
sx q[0];
rz(-0.14725421) q[0];
rz(-2.9946949) q[1];
sx q[1];
rz(-2.2203827) q[1];
sx q[1];
rz(1.1796835) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0834107) q[0];
sx q[0];
rz(-1.9866052) q[0];
sx q[0];
rz(1.2021628) q[0];
rz(-0.83663655) q[2];
sx q[2];
rz(-0.61021611) q[2];
sx q[2];
rz(1.1532551) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2218297) q[1];
sx q[1];
rz(-2.5357863) q[1];
sx q[1];
rz(-2.9472652) q[1];
rz(-pi) q[2];
rz(0.75408077) q[3];
sx q[3];
rz(-1.6508266) q[3];
sx q[3];
rz(-2.4802993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2008721) q[2];
sx q[2];
rz(-0.93738896) q[2];
sx q[2];
rz(-0.76784039) q[2];
rz(-2.614894) q[3];
sx q[3];
rz(-1.0517164) q[3];
sx q[3];
rz(2.5832978) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.679477) q[0];
sx q[0];
rz(-2.9556892) q[0];
sx q[0];
rz(1.0900849) q[0];
rz(-1.7915626) q[1];
sx q[1];
rz(-1.4005902) q[1];
sx q[1];
rz(-2.7166691) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23608828) q[0];
sx q[0];
rz(-2.6668735) q[0];
sx q[0];
rz(1.4464054) q[0];
x q[1];
rz(-3.0504136) q[2];
sx q[2];
rz(-0.65185368) q[2];
sx q[2];
rz(-0.21287795) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7817345) q[1];
sx q[1];
rz(-1.1175851) q[1];
sx q[1];
rz(-1.9537326) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0474237) q[3];
sx q[3];
rz(-0.93195385) q[3];
sx q[3];
rz(-2.1372319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2316042) q[2];
sx q[2];
rz(-0.51077545) q[2];
sx q[2];
rz(2.6166022) q[2];
rz(-1.3548343) q[3];
sx q[3];
rz(-2.2460008) q[3];
sx q[3];
rz(1.3625328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.7056284) q[0];
sx q[0];
rz(-2.4810677) q[0];
sx q[0];
rz(0.10667644) q[0];
rz(-1.0542997) q[1];
sx q[1];
rz(-1.432447) q[1];
sx q[1];
rz(1.7342825) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6808436) q[0];
sx q[0];
rz(-1.1011194) q[0];
sx q[0];
rz(-0.90523714) q[0];
rz(-pi) q[1];
rz(-0.65347341) q[2];
sx q[2];
rz(-1.1879041) q[2];
sx q[2];
rz(0.083195638) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.095276621) q[1];
sx q[1];
rz(-1.7937978) q[1];
sx q[1];
rz(2.9848802) q[1];
rz(-pi) q[2];
rz(-2.484177) q[3];
sx q[3];
rz(-0.86674009) q[3];
sx q[3];
rz(-2.2179066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6537031) q[2];
sx q[2];
rz(-0.50450486) q[2];
sx q[2];
rz(2.8686236) q[2];
rz(-0.87131396) q[3];
sx q[3];
rz(-1.6815691) q[3];
sx q[3];
rz(3.0070846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3684261) q[0];
sx q[0];
rz(-1.9482524) q[0];
sx q[0];
rz(0.85309749) q[0];
rz(2.7586965) q[1];
sx q[1];
rz(-1.8538414) q[1];
sx q[1];
rz(-0.014898653) q[1];
rz(-0.64528633) q[2];
sx q[2];
rz(-0.66678478) q[2];
sx q[2];
rz(-0.18258584) q[2];
rz(-1.4703582) q[3];
sx q[3];
rz(-2.3681869) q[3];
sx q[3];
rz(-2.8988732) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
