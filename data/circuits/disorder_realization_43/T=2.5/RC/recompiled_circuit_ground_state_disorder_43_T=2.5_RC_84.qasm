OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3620152) q[0];
sx q[0];
rz(-2.9147122) q[0];
sx q[0];
rz(1.3751295) q[0];
rz(-0.094376266) q[1];
sx q[1];
rz(-2.2278995) q[1];
sx q[1];
rz(2.8606666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.323712) q[0];
sx q[0];
rz(-1.4120988) q[0];
sx q[0];
rz(0.90526064) q[0];
x q[1];
rz(1.9380577) q[2];
sx q[2];
rz(-0.31031552) q[2];
sx q[2];
rz(-2.4825437) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62796616) q[1];
sx q[1];
rz(-1.1556323) q[1];
sx q[1];
rz(-2.8529851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6219764) q[3];
sx q[3];
rz(-2.5027788) q[3];
sx q[3];
rz(2.6878302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1137587) q[2];
sx q[2];
rz(-1.7366624) q[2];
sx q[2];
rz(0.11715451) q[2];
rz(-0.33200085) q[3];
sx q[3];
rz(-2.3947075) q[3];
sx q[3];
rz(-0.78506708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
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
rz(2.4431045) q[0];
sx q[0];
rz(-2.3790058) q[0];
sx q[0];
rz(-2.7681328) q[0];
rz(-0.39730486) q[1];
sx q[1];
rz(-1.9970857) q[1];
sx q[1];
rz(-0.73748803) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1339379) q[0];
sx q[0];
rz(-0.93385252) q[0];
sx q[0];
rz(0.68151125) q[0];
rz(-pi) q[1];
rz(2.5869484) q[2];
sx q[2];
rz(-2.0337542) q[2];
sx q[2];
rz(2.2210768) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0651689) q[1];
sx q[1];
rz(-0.57145125) q[1];
sx q[1];
rz(2.854611) q[1];
rz(-pi) q[2];
rz(1.2190388) q[3];
sx q[3];
rz(-0.84461194) q[3];
sx q[3];
rz(1.5016457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.071659) q[2];
sx q[2];
rz(-1.5295014) q[2];
sx q[2];
rz(2.7640479) q[2];
rz(0.30970445) q[3];
sx q[3];
rz(-1.0891424) q[3];
sx q[3];
rz(-1.6449876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51652235) q[0];
sx q[0];
rz(-0.83928883) q[0];
sx q[0];
rz(1.9260433) q[0];
rz(-2.1274321) q[1];
sx q[1];
rz(-1.4233669) q[1];
sx q[1];
rz(-0.97022143) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.559186) q[0];
sx q[0];
rz(-1.0716039) q[0];
sx q[0];
rz(0.046120709) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90486352) q[2];
sx q[2];
rz(-2.3653125) q[2];
sx q[2];
rz(1.5011476) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5287703) q[1];
sx q[1];
rz(-2.7547358) q[1];
sx q[1];
rz(1.0861138) q[1];
x q[2];
rz(0.35200624) q[3];
sx q[3];
rz(-1.6791293) q[3];
sx q[3];
rz(-1.7217404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0479451) q[2];
sx q[2];
rz(-0.84257546) q[2];
sx q[2];
rz(-2.688431) q[2];
rz(1.2293182) q[3];
sx q[3];
rz(-1.2776351) q[3];
sx q[3];
rz(-0.17190988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9446843) q[0];
sx q[0];
rz(-2.4561645) q[0];
sx q[0];
rz(-0.36566439) q[0];
rz(0.91224313) q[1];
sx q[1];
rz(-1.279) q[1];
sx q[1];
rz(-2.7361187) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376108) q[0];
sx q[0];
rz(-1.77586) q[0];
sx q[0];
rz(-3.0201941) q[0];
rz(-0.75587749) q[2];
sx q[2];
rz(-1.4149932) q[2];
sx q[2];
rz(0.98658094) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4279086) q[1];
sx q[1];
rz(-2.1024087) q[1];
sx q[1];
rz(-1.3757399) q[1];
rz(-2.7248091) q[3];
sx q[3];
rz(-2.458771) q[3];
sx q[3];
rz(-2.5522751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.071986467) q[2];
sx q[2];
rz(-1.7013841) q[2];
sx q[2];
rz(-1.5988505) q[2];
rz(0.37374464) q[3];
sx q[3];
rz(-1.6662686) q[3];
sx q[3];
rz(-1.6811949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54295802) q[0];
sx q[0];
rz(-2.3668508) q[0];
sx q[0];
rz(-0.69611088) q[0];
rz(-1.2593345) q[1];
sx q[1];
rz(-2.0399317) q[1];
sx q[1];
rz(-0.36516821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29699638) q[0];
sx q[0];
rz(-2.0005324) q[0];
sx q[0];
rz(-0.55083042) q[0];
rz(-pi) q[1];
rz(-2.2155207) q[2];
sx q[2];
rz(-2.640291) q[2];
sx q[2];
rz(2.9013366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2272353) q[1];
sx q[1];
rz(-1.0252756) q[1];
sx q[1];
rz(-1.5724214) q[1];
rz(2.5138084) q[3];
sx q[3];
rz(-1.2300228) q[3];
sx q[3];
rz(0.22751156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.460707) q[2];
sx q[2];
rz(-1.978771) q[2];
sx q[2];
rz(-0.17670259) q[2];
rz(-1.3867311) q[3];
sx q[3];
rz(-1.0597798) q[3];
sx q[3];
rz(0.37469125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6143167) q[0];
sx q[0];
rz(-2.2820331) q[0];
sx q[0];
rz(1.5981307) q[0];
rz(1.3044926) q[1];
sx q[1];
rz(-2.0871128) q[1];
sx q[1];
rz(-2.3354882) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4325754) q[0];
sx q[0];
rz(-2.3516549) q[0];
sx q[0];
rz(-2.8054601) q[0];
x q[1];
rz(-0.88509934) q[2];
sx q[2];
rz(-0.72413596) q[2];
sx q[2];
rz(-0.76937461) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42000586) q[1];
sx q[1];
rz(-1.3665821) q[1];
sx q[1];
rz(-0.050724357) q[1];
x q[2];
rz(-0.73202912) q[3];
sx q[3];
rz(-0.71683466) q[3];
sx q[3];
rz(-0.88535374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1690037) q[2];
sx q[2];
rz(-2.6469595) q[2];
sx q[2];
rz(-0.063610323) q[2];
rz(-0.58498597) q[3];
sx q[3];
rz(-1.4002742) q[3];
sx q[3];
rz(1.6271648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99120283) q[0];
sx q[0];
rz(-1.3911893) q[0];
sx q[0];
rz(2.4124131) q[0];
rz(-1.9623307) q[1];
sx q[1];
rz(-0.74960342) q[1];
sx q[1];
rz(2.8470305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7594265) q[0];
sx q[0];
rz(-1.3148892) q[0];
sx q[0];
rz(-2.6271174) q[0];
rz(2.0038347) q[2];
sx q[2];
rz(-1.6707509) q[2];
sx q[2];
rz(-1.2765027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2475491) q[1];
sx q[1];
rz(-1.4534543) q[1];
sx q[1];
rz(-2.3168922) q[1];
rz(-pi) q[2];
rz(0.90055777) q[3];
sx q[3];
rz(-1.2646741) q[3];
sx q[3];
rz(-2.8365119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6509167) q[2];
sx q[2];
rz(-1.4399521) q[2];
sx q[2];
rz(0.73224625) q[2];
rz(-0.8693153) q[3];
sx q[3];
rz(-1.9711875) q[3];
sx q[3];
rz(-1.7349582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9100087) q[0];
sx q[0];
rz(-1.5861479) q[0];
sx q[0];
rz(-0.79767942) q[0];
rz(0.32294598) q[1];
sx q[1];
rz(-0.99635092) q[1];
sx q[1];
rz(-1.4083883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8033215) q[0];
sx q[0];
rz(-0.59280286) q[0];
sx q[0];
rz(0.37180255) q[0];
rz(-1.1118719) q[2];
sx q[2];
rz(-1.7844229) q[2];
sx q[2];
rz(-1.2101419) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1221315) q[1];
sx q[1];
rz(-0.82585649) q[1];
sx q[1];
rz(-1.248658) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6334535) q[3];
sx q[3];
rz(-0.36341011) q[3];
sx q[3];
rz(-1.5606708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8352167) q[2];
sx q[2];
rz(-2.0312467) q[2];
sx q[2];
rz(-0.19732538) q[2];
rz(2.3399682) q[3];
sx q[3];
rz(-2.1620731) q[3];
sx q[3];
rz(2.7200123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.1590969) q[0];
sx q[0];
rz(-1.4747341) q[0];
sx q[0];
rz(0.91484797) q[0];
rz(0.21028701) q[1];
sx q[1];
rz(-0.57840127) q[1];
sx q[1];
rz(-0.68967825) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.698477) q[0];
sx q[0];
rz(-1.7962828) q[0];
sx q[0];
rz(1.6661635) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4532926) q[2];
sx q[2];
rz(-1.5853865) q[2];
sx q[2];
rz(2.7772825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9079118) q[1];
sx q[1];
rz(-2.4123976) q[1];
sx q[1];
rz(-2.3883874) q[1];
rz(-2.7977169) q[3];
sx q[3];
rz(-1.0707885) q[3];
sx q[3];
rz(0.22658928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0476394) q[2];
sx q[2];
rz(-0.84838212) q[2];
sx q[2];
rz(-2.8116255) q[2];
rz(1.0432358) q[3];
sx q[3];
rz(-1.7236575) q[3];
sx q[3];
rz(-0.51658336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.030815) q[0];
sx q[0];
rz(-1.3309706) q[0];
sx q[0];
rz(-2.0844841) q[0];
rz(2.8874176) q[1];
sx q[1];
rz(-1.1421685) q[1];
sx q[1];
rz(0.70456299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1833558) q[0];
sx q[0];
rz(-0.57794774) q[0];
sx q[0];
rz(-2.1653752) q[0];
rz(-pi) q[1];
rz(-1.3866587) q[2];
sx q[2];
rz(-1.6312851) q[2];
sx q[2];
rz(-0.16512251) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3290289) q[1];
sx q[1];
rz(-0.48406752) q[1];
sx q[1];
rz(1.9639534) q[1];
x q[2];
rz(-0.71263616) q[3];
sx q[3];
rz(-0.27240005) q[3];
sx q[3];
rz(0.89100641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5083984) q[2];
sx q[2];
rz(-1.9693547) q[2];
sx q[2];
rz(2.634826) q[2];
rz(2.9554328) q[3];
sx q[3];
rz(-2.9045744) q[3];
sx q[3];
rz(-0.53562927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(1.7645466) q[0];
sx q[0];
rz(-1.6896387) q[0];
sx q[0];
rz(1.1402546) q[0];
rz(2.2346732) q[1];
sx q[1];
rz(-0.74105558) q[1];
sx q[1];
rz(1.8190609) q[1];
rz(2.604031) q[2];
sx q[2];
rz(-0.96502177) q[2];
sx q[2];
rz(1.6966664) q[2];
rz(-0.68578699) q[3];
sx q[3];
rz(-2.7332173) q[3];
sx q[3];
rz(0.12242534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
