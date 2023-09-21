OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(-1.1146201) q[0];
sx q[0];
rz(-0.00014076509) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(4.1058022) q[1];
sx q[1];
rz(10.618187) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0573187) q[0];
sx q[0];
rz(-2.8017375) q[0];
sx q[0];
rz(-1.0124595) q[0];
rz(-pi) q[1];
rz(-0.46618669) q[2];
sx q[2];
rz(-2.5417915) q[2];
sx q[2];
rz(-0.28238645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5692917) q[1];
sx q[1];
rz(-0.83218677) q[1];
sx q[1];
rz(2.4898847) q[1];
x q[2];
rz(1.7883349) q[3];
sx q[3];
rz(-1.4635524) q[3];
sx q[3];
rz(-1.6579962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.45941916) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(-1.9127282) q[2];
rz(-1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(-1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380149) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(-3.1139328) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(-2.0181296) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7682122) q[0];
sx q[0];
rz(-2.0313615) q[0];
sx q[0];
rz(-3.0773613) q[0];
rz(2.3617619) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(-1.1080527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9532721) q[1];
sx q[1];
rz(-2.3768432) q[1];
sx q[1];
rz(-2.3482167) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1334322) q[3];
sx q[3];
rz(-2.1481272) q[3];
sx q[3];
rz(-1.2853704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79364395) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(-0.91903764) q[2];
rz(-2.4675026) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27750257) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(-1.8664237) q[0];
rz(-0.69349849) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(-1.1330053) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1904859) q[0];
sx q[0];
rz(-1.5674942) q[0];
sx q[0];
rz(1.7130873) q[0];
rz(-pi) q[1];
rz(0.79046952) q[2];
sx q[2];
rz(-1.1970453) q[2];
sx q[2];
rz(0.30465301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36724597) q[1];
sx q[1];
rz(-2.5173752) q[1];
sx q[1];
rz(-1.063785) q[1];
x q[2];
rz(0.95136178) q[3];
sx q[3];
rz(-1.0361443) q[3];
sx q[3];
rz(1.9922647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.8481002) q[2];
rz(-3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816198) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(1.1608634) q[0];
rz(-2.2456031) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(-0.13555759) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7134705) q[0];
sx q[0];
rz(-0.8849511) q[0];
sx q[0];
rz(-1.9786406) q[0];
rz(-pi) q[1];
rz(-1.1866456) q[2];
sx q[2];
rz(-1.5015366) q[2];
sx q[2];
rz(-2.4406976) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1783274) q[1];
sx q[1];
rz(-1.7566924) q[1];
sx q[1];
rz(-2.9796897) q[1];
rz(2.9837708) q[3];
sx q[3];
rz(-1.8421679) q[3];
sx q[3];
rz(-1.2543169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(-0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0376461) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(2.1283545) q[0];
rz(-3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(2.0577046) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5686544) q[0];
sx q[0];
rz(-1.3875811) q[0];
sx q[0];
rz(-1.3230447) q[0];
rz(-0.82917825) q[2];
sx q[2];
rz(-0.36311705) q[2];
sx q[2];
rz(1.6104289) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14099546) q[1];
sx q[1];
rz(-2.6037569) q[1];
sx q[1];
rz(1.3586033) q[1];
rz(-pi) q[2];
rz(1.6211987) q[3];
sx q[3];
rz(-2.084123) q[3];
sx q[3];
rz(0.36171519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9087387) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(0.24442913) q[2];
rz(0.43236732) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(0.50306815) q[3];
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
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844834) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(0.094141468) q[0];
rz(2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(-2.24618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5603134) q[0];
sx q[0];
rz(-0.34438294) q[0];
sx q[0];
rz(0.11238213) q[0];
rz(1.3307829) q[2];
sx q[2];
rz(-2.1905106) q[2];
sx q[2];
rz(1.0735219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4737282) q[1];
sx q[1];
rz(-1.5267936) q[1];
sx q[1];
rz(1.4021137) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.496109) q[3];
sx q[3];
rz(-1.5731305) q[3];
sx q[3];
rz(-2.1722349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-2.3383979) q[2];
rz(1.1903654) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(-0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068709277) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(-0.51914006) q[0];
rz(2.5601162) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(-1.8849467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089766895) q[0];
sx q[0];
rz(-1.8221812) q[0];
sx q[0];
rz(-3.0999523) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9728327) q[2];
sx q[2];
rz(-2.6425344) q[2];
sx q[2];
rz(-1.4025584) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9650967) q[1];
sx q[1];
rz(-2.8526222) q[1];
sx q[1];
rz(-0.22775905) q[1];
rz(-pi) q[2];
rz(-1.2506966) q[3];
sx q[3];
rz(-2.2543636) q[3];
sx q[3];
rz(0.86910955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98823035) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(1.3640277) q[2];
rz(-0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0751188) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(0.30817729) q[0];
rz(-3.0691052) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(-2.7546308) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68650866) q[0];
sx q[0];
rz(-0.63502705) q[0];
sx q[0];
rz(-1.1839266) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2690291) q[2];
sx q[2];
rz(-1.9930895) q[2];
sx q[2];
rz(-2.1625105) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-2.4657001) q[1];
sx q[1];
rz(3.0216316) q[1];
rz(1.79027) q[3];
sx q[3];
rz(-1.0843127) q[3];
sx q[3];
rz(2.6284077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62362921) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(-0.44000885) q[2];
rz(2.4258339) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(2.0619152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(-2.0429042) q[0];
rz(-2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-0.049302014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7002174) q[0];
sx q[0];
rz(-0.8964296) q[0];
sx q[0];
rz(-2.3250439) q[0];
x q[1];
rz(2.1575035) q[2];
sx q[2];
rz(-0.77312914) q[2];
sx q[2];
rz(0.099260515) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.94433632) q[1];
sx q[1];
rz(-1.1068871) q[1];
sx q[1];
rz(2.0891561) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0588166) q[3];
sx q[3];
rz(-0.89142311) q[3];
sx q[3];
rz(-2.8454012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.07842841) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(-0.72193974) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14770517) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(1.0797427) q[0];
rz(-1.059277) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(1.7396897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3106874) q[0];
sx q[0];
rz(-1.3584104) q[0];
sx q[0];
rz(3.0124245) q[0];
rz(2.9548168) q[2];
sx q[2];
rz(-1.5439856) q[2];
sx q[2];
rz(1.9865799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2682174) q[1];
sx q[1];
rz(-0.89726071) q[1];
sx q[1];
rz(2.7482277) q[1];
rz(-pi) q[2];
x q[2];
rz(2.510342) q[3];
sx q[3];
rz(-1.4025941) q[3];
sx q[3];
rz(1.9025308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4828651) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(1.6213017) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(-2.4441392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.993492) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(2.2254754) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(-0.82952164) q[2];
sx q[2];
rz(-0.98486949) q[2];
sx q[2];
rz(-2.9688901) q[2];
rz(-1.0536853) q[3];
sx q[3];
rz(-1.5862982) q[3];
sx q[3];
rz(1.3261212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];