OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.779939) q[0];
sx q[0];
rz(-0.10804636) q[0];
sx q[0];
rz(1.894423) q[0];
rz(-0.78020686) q[1];
sx q[1];
rz(8.297774) q[1];
sx q[1];
rz(10.170593) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0184263) q[0];
sx q[0];
rz(-1.9115051) q[0];
sx q[0];
rz(-0.19827224) q[0];
rz(1.4100513) q[2];
sx q[2];
rz(-1.5622592) q[2];
sx q[2];
rz(1.05913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31065658) q[1];
sx q[1];
rz(-1.6864221) q[1];
sx q[1];
rz(-1.7831037) q[1];
rz(-pi) q[2];
rz(-2.8549544) q[3];
sx q[3];
rz(-2.6034546) q[3];
sx q[3];
rz(1.2179483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1625533) q[2];
sx q[2];
rz(-2.5472842) q[2];
sx q[2];
rz(1.7621367) q[2];
rz(2.4123794) q[3];
sx q[3];
rz(-0.76494923) q[3];
sx q[3];
rz(-2.2123857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58296975) q[0];
sx q[0];
rz(-2.0780777) q[0];
sx q[0];
rz(2.9003918) q[0];
rz(-0.25602117) q[1];
sx q[1];
rz(-2.119901) q[1];
sx q[1];
rz(2.3604438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22750073) q[0];
sx q[0];
rz(-2.1315398) q[0];
sx q[0];
rz(0.95277159) q[0];
rz(-2.5767869) q[2];
sx q[2];
rz(-1.6132659) q[2];
sx q[2];
rz(0.04893411) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1835412) q[1];
sx q[1];
rz(-2.1528917) q[1];
sx q[1];
rz(2.2954825) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2595176) q[3];
sx q[3];
rz(-0.14434179) q[3];
sx q[3];
rz(2.8022769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.29391924) q[2];
sx q[2];
rz(-1.5519698) q[2];
sx q[2];
rz(1.5576564) q[2];
rz(3.1368351) q[3];
sx q[3];
rz(-0.59499756) q[3];
sx q[3];
rz(0.34576542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3874338) q[0];
sx q[0];
rz(-1.4492946) q[0];
sx q[0];
rz(-2.9505728) q[0];
rz(0.35107958) q[1];
sx q[1];
rz(-2.7103238) q[1];
sx q[1];
rz(-1.692159) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96286839) q[0];
sx q[0];
rz(-2.7984848) q[0];
sx q[0];
rz(2.9897811) q[0];
x q[1];
rz(-2.6963553) q[2];
sx q[2];
rz(-1.9539333) q[2];
sx q[2];
rz(-0.68510011) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.23891029) q[1];
sx q[1];
rz(-2.604831) q[1];
sx q[1];
rz(-2.6031659) q[1];
rz(0.48630096) q[3];
sx q[3];
rz(-0.6502021) q[3];
sx q[3];
rz(-2.7335258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8423975) q[2];
sx q[2];
rz(-1.7123875) q[2];
sx q[2];
rz(0.12796417) q[2];
rz(-1.1658824) q[3];
sx q[3];
rz(-2.2539625) q[3];
sx q[3];
rz(-2.8036346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8562451) q[0];
sx q[0];
rz(-1.0839533) q[0];
sx q[0];
rz(-1.4803084) q[0];
rz(3.0604494) q[1];
sx q[1];
rz(-2.0442043) q[1];
sx q[1];
rz(-2.2611484) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1937516) q[0];
sx q[0];
rz(-1.5955076) q[0];
sx q[0];
rz(1.5773043) q[0];
rz(-pi) q[1];
x q[1];
rz(1.501785) q[2];
sx q[2];
rz(-0.61997783) q[2];
sx q[2];
rz(-0.91199707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2374946) q[1];
sx q[1];
rz(-1.6818871) q[1];
sx q[1];
rz(2.6059125) q[1];
rz(-3.1034558) q[3];
sx q[3];
rz(-0.55805579) q[3];
sx q[3];
rz(-1.7779153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60260281) q[2];
sx q[2];
rz(-2.4675641) q[2];
sx q[2];
rz(-1.7146141) q[2];
rz(-1.9849298) q[3];
sx q[3];
rz(-1.7159228) q[3];
sx q[3];
rz(-0.15779933) q[3];
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
rz(-pi/2) q[0];
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
rz(-0.91931525) q[0];
sx q[0];
rz(-2.3265525) q[0];
sx q[0];
rz(2.2015233) q[0];
rz(-2.6433511) q[1];
sx q[1];
rz(-0.91612852) q[1];
sx q[1];
rz(1.1246276) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4152545) q[0];
sx q[0];
rz(-2.4372516) q[0];
sx q[0];
rz(1.5096774) q[0];
rz(-1.1246936) q[2];
sx q[2];
rz(-0.44844018) q[2];
sx q[2];
rz(-2.5288127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7371205) q[1];
sx q[1];
rz(-1.1386765) q[1];
sx q[1];
rz(2.4678556) q[1];
rz(-pi) q[2];
rz(0.78557555) q[3];
sx q[3];
rz(-1.5557369) q[3];
sx q[3];
rz(0.3846947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70790946) q[2];
sx q[2];
rz(-0.83432546) q[2];
sx q[2];
rz(-1.4617807) q[2];
rz(0.24623571) q[3];
sx q[3];
rz(-2.2610531) q[3];
sx q[3];
rz(-1.0730526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4435302) q[0];
sx q[0];
rz(-1.9140697) q[0];
sx q[0];
rz(-0.45502934) q[0];
rz(-0.21806923) q[1];
sx q[1];
rz(-0.90556216) q[1];
sx q[1];
rz(-2.6630482) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9537331) q[0];
sx q[0];
rz(-1.4728913) q[0];
sx q[0];
rz(-2.4763473) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9265106) q[2];
sx q[2];
rz(-0.96181574) q[2];
sx q[2];
rz(0.019542309) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8763833) q[1];
sx q[1];
rz(-0.57683101) q[1];
sx q[1];
rz(2.2028707) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1699311) q[3];
sx q[3];
rz(-2.1671579) q[3];
sx q[3];
rz(-0.14113035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1856508) q[2];
sx q[2];
rz(-1.6472683) q[2];
sx q[2];
rz(-0.68598023) q[2];
rz(1.3395122) q[3];
sx q[3];
rz(-1.1687665) q[3];
sx q[3];
rz(-1.3919938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5865536) q[0];
sx q[0];
rz(-1.3864484) q[0];
sx q[0];
rz(-2.6626124) q[0];
rz(-0.85887495) q[1];
sx q[1];
rz(-0.54194599) q[1];
sx q[1];
rz(-3.0368793) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0514487) q[0];
sx q[0];
rz(-1.3810754) q[0];
sx q[0];
rz(-2.7499697) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2042047) q[2];
sx q[2];
rz(-1.1127923) q[2];
sx q[2];
rz(1.7228433) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0942451) q[1];
sx q[1];
rz(-2.1434103) q[1];
sx q[1];
rz(0.25783536) q[1];
rz(-pi) q[2];
rz(2.0488304) q[3];
sx q[3];
rz(-1.5198738) q[3];
sx q[3];
rz(1.6855406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0571478) q[2];
sx q[2];
rz(-1.8859325) q[2];
sx q[2];
rz(-2.2596333) q[2];
rz(0.070934892) q[3];
sx q[3];
rz(-1.8788012) q[3];
sx q[3];
rz(-0.96496636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1607745) q[0];
sx q[0];
rz(-0.47176281) q[0];
sx q[0];
rz(-1.8881352) q[0];
rz(-0.71022931) q[1];
sx q[1];
rz(-1.9435147) q[1];
sx q[1];
rz(1.089878) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37140815) q[0];
sx q[0];
rz(-0.84781269) q[0];
sx q[0];
rz(0.56130479) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5023191) q[2];
sx q[2];
rz(-2.6007923) q[2];
sx q[2];
rz(-1.6317473) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0347669) q[1];
sx q[1];
rz(-2.2112101) q[1];
sx q[1];
rz(-0.57042112) q[1];
x q[2];
rz(-2.5651781) q[3];
sx q[3];
rz(-1.4428328) q[3];
sx q[3];
rz(1.1157677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39500427) q[2];
sx q[2];
rz(-2.2825664) q[2];
sx q[2];
rz(-2.11917) q[2];
rz(0.59631452) q[3];
sx q[3];
rz(-0.78323451) q[3];
sx q[3];
rz(0.35308009) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31350598) q[0];
sx q[0];
rz(-2.6441898) q[0];
sx q[0];
rz(-1.5561546) q[0];
rz(-1.4900788) q[1];
sx q[1];
rz(-0.45526344) q[1];
sx q[1];
rz(-0.99686399) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5307025) q[0];
sx q[0];
rz(-0.35323745) q[0];
sx q[0];
rz(0.20521407) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7559986) q[2];
sx q[2];
rz(-1.5723937) q[2];
sx q[2];
rz(-1.8271556) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91040033) q[1];
sx q[1];
rz(-2.2072189) q[1];
sx q[1];
rz(0.6271805) q[1];
x q[2];
rz(-2.1622657) q[3];
sx q[3];
rz(-1.1445657) q[3];
sx q[3];
rz(-1.1104079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8263714) q[2];
sx q[2];
rz(-2.4085277) q[2];
sx q[2];
rz(0.23615393) q[2];
rz(-0.30596966) q[3];
sx q[3];
rz(-1.1924084) q[3];
sx q[3];
rz(-1.6570305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.21208256) q[0];
sx q[0];
rz(-0.65123737) q[0];
sx q[0];
rz(0.65810743) q[0];
rz(-1.2492389) q[1];
sx q[1];
rz(-2.55195) q[1];
sx q[1];
rz(2.6499937) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8458873) q[0];
sx q[0];
rz(-1.0330366) q[0];
sx q[0];
rz(-2.0989492) q[0];
x q[1];
rz(-1.7936321) q[2];
sx q[2];
rz(-2.1249944) q[2];
sx q[2];
rz(-1.8775307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1695849) q[1];
sx q[1];
rz(-1.8465202) q[1];
sx q[1];
rz(-1.5805401) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.059547) q[3];
sx q[3];
rz(-1.4982035) q[3];
sx q[3];
rz(-1.5047764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7813985) q[2];
sx q[2];
rz(-0.34803826) q[2];
sx q[2];
rz(1.4813102) q[2];
rz(-0.71634746) q[3];
sx q[3];
rz(-2.4951388) q[3];
sx q[3];
rz(-0.20814482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79138712) q[0];
sx q[0];
rz(-1.5972142) q[0];
sx q[0];
rz(0.41440339) q[0];
rz(2.1898337) q[1];
sx q[1];
rz(-1.2928243) q[1];
sx q[1];
rz(-1.181319) q[1];
rz(2.4454115) q[2];
sx q[2];
rz(-1.0715108) q[2];
sx q[2];
rz(0.61658695) q[2];
rz(1.7813206) q[3];
sx q[3];
rz(-0.68974907) q[3];
sx q[3];
rz(-0.30221119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
