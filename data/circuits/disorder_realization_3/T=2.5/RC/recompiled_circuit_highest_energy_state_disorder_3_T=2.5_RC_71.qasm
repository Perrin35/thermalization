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
rz(0.74547493) q[0];
sx q[0];
rz(-0.59362721) q[0];
sx q[0];
rz(-2.797085) q[0];
rz(-1.966882) q[1];
sx q[1];
rz(-0.36829683) q[1];
sx q[1];
rz(2.5370497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6942351) q[0];
sx q[0];
rz(-1.6295054) q[0];
sx q[0];
rz(2.319058) q[0];
rz(-pi) q[1];
rz(0.069434631) q[2];
sx q[2];
rz(-1.3926818) q[2];
sx q[2];
rz(1.7770065) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82935059) q[1];
sx q[1];
rz(-1.568746) q[1];
sx q[1];
rz(-2.7228628) q[1];
rz(-0.80896583) q[3];
sx q[3];
rz(-2.3936317) q[3];
sx q[3];
rz(0.35424074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.979226) q[2];
sx q[2];
rz(-1.484885) q[2];
sx q[2];
rz(-1.7171198) q[2];
rz(-3.034397) q[3];
sx q[3];
rz(-0.19459477) q[3];
sx q[3];
rz(-2.181982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87218881) q[0];
sx q[0];
rz(-0.93282455) q[0];
sx q[0];
rz(1.3099571) q[0];
rz(0.45920363) q[1];
sx q[1];
rz(-1.867086) q[1];
sx q[1];
rz(-0.31712636) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91823365) q[0];
sx q[0];
rz(-0.14592136) q[0];
sx q[0];
rz(-2.8097665) q[0];
x q[1];
rz(-0.049811157) q[2];
sx q[2];
rz(-1.5684109) q[2];
sx q[2];
rz(2.002169) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9307946) q[1];
sx q[1];
rz(-2.5902777) q[1];
sx q[1];
rz(3.0153794) q[1];
x q[2];
rz(1.7789654) q[3];
sx q[3];
rz(-0.82566264) q[3];
sx q[3];
rz(-2.7221808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.7121048) q[2];
sx q[2];
rz(-1.2965344) q[2];
sx q[2];
rz(-0.8075766) q[2];
rz(0.56337041) q[3];
sx q[3];
rz(-2.6475776) q[3];
sx q[3];
rz(-1.0529244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5092369) q[0];
sx q[0];
rz(-0.18993264) q[0];
sx q[0];
rz(-2.307039) q[0];
rz(2.6368311) q[1];
sx q[1];
rz(-2.7752462) q[1];
sx q[1];
rz(1.4588446) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6989649) q[0];
sx q[0];
rz(-2.1219398) q[0];
sx q[0];
rz(-1.3359265) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3332852) q[2];
sx q[2];
rz(-0.9918074) q[2];
sx q[2];
rz(2.4579687) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7968081) q[1];
sx q[1];
rz(-2.4151509) q[1];
sx q[1];
rz(-0.11747719) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57679983) q[3];
sx q[3];
rz(-2.24772) q[3];
sx q[3];
rz(-3.0056382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0316281) q[2];
sx q[2];
rz(-1.4723023) q[2];
sx q[2];
rz(-2.7785981) q[2];
rz(1.1686769) q[3];
sx q[3];
rz(-2.3790338) q[3];
sx q[3];
rz(2.3762083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7211683) q[0];
sx q[0];
rz(-0.52614251) q[0];
sx q[0];
rz(-2.5508733) q[0];
rz(-0.72714725) q[1];
sx q[1];
rz(-0.37494451) q[1];
sx q[1];
rz(-0.66853833) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1626133) q[0];
sx q[0];
rz(-0.091446459) q[0];
sx q[0];
rz(1.0981111) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0425778) q[2];
sx q[2];
rz(-2.2792247) q[2];
sx q[2];
rz(2.0076795) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8188386) q[1];
sx q[1];
rz(-1.1024945) q[1];
sx q[1];
rz(0.78264758) q[1];
rz(2.0987857) q[3];
sx q[3];
rz(-0.90195459) q[3];
sx q[3];
rz(-0.5878512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.10355243) q[2];
sx q[2];
rz(-1.870564) q[2];
sx q[2];
rz(1.5070149) q[2];
rz(0.42190894) q[3];
sx q[3];
rz(-0.76015893) q[3];
sx q[3];
rz(0.71934492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5108532) q[0];
sx q[0];
rz(-0.35950867) q[0];
sx q[0];
rz(-2.0810293) q[0];
rz(3.0781436) q[1];
sx q[1];
rz(-2.5109406) q[1];
sx q[1];
rz(-0.55387703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3097274) q[0];
sx q[0];
rz(-2.1500282) q[0];
sx q[0];
rz(2.6088891) q[0];
x q[1];
rz(-0.92775821) q[2];
sx q[2];
rz(-2.0915439) q[2];
sx q[2];
rz(-0.82633229) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.088515048) q[1];
sx q[1];
rz(-1.6106399) q[1];
sx q[1];
rz(0.30327176) q[1];
rz(3.0054566) q[3];
sx q[3];
rz(-1.247974) q[3];
sx q[3];
rz(-0.28188595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.267103) q[2];
sx q[2];
rz(-2.5449982) q[2];
sx q[2];
rz(-0.60274094) q[2];
rz(-0.069325773) q[3];
sx q[3];
rz(-1.5452789) q[3];
sx q[3];
rz(0.16665211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.199274) q[0];
sx q[0];
rz(-2.5397904) q[0];
sx q[0];
rz(-0.45403516) q[0];
rz(1.8858689) q[1];
sx q[1];
rz(-2.1818826) q[1];
sx q[1];
rz(1.4383291) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5471449) q[0];
sx q[0];
rz(-1.1714352) q[0];
sx q[0];
rz(2.2844676) q[0];
rz(-pi) q[1];
rz(-2.8642162) q[2];
sx q[2];
rz(-1.2834594) q[2];
sx q[2];
rz(-0.6544906) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.033129902) q[1];
sx q[1];
rz(-2.3537043) q[1];
sx q[1];
rz(2.2910396) q[1];
rz(-pi) q[2];
rz(1.0694396) q[3];
sx q[3];
rz(-0.92235288) q[3];
sx q[3];
rz(1.7164149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85256514) q[2];
sx q[2];
rz(-1.270741) q[2];
sx q[2];
rz(-1.3555869) q[2];
rz(-1.2191314) q[3];
sx q[3];
rz(-1.6617323) q[3];
sx q[3];
rz(1.5752972) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3179625) q[0];
sx q[0];
rz(-0.83331236) q[0];
sx q[0];
rz(1.3939567) q[0];
rz(-2.6849003) q[1];
sx q[1];
rz(-0.87867457) q[1];
sx q[1];
rz(-1.7535694) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0972892) q[0];
sx q[0];
rz(-0.40709041) q[0];
sx q[0];
rz(1.2319698) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51808234) q[2];
sx q[2];
rz(-0.23031313) q[2];
sx q[2];
rz(1.7948593) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79560515) q[1];
sx q[1];
rz(-1.8966764) q[1];
sx q[1];
rz(-1.4326976) q[1];
rz(2.4782466) q[3];
sx q[3];
rz(-2.4069549) q[3];
sx q[3];
rz(-1.2879077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2726511) q[2];
sx q[2];
rz(-2.36918) q[2];
sx q[2];
rz(-0.19115494) q[2];
rz(-1.8027421) q[3];
sx q[3];
rz(-0.50722417) q[3];
sx q[3];
rz(0.52201456) q[3];
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
rz(-2.1538447) q[0];
sx q[0];
rz(-0.67180434) q[0];
sx q[0];
rz(-1.9195358) q[0];
rz(-0.98474312) q[1];
sx q[1];
rz(-2.1538815) q[1];
sx q[1];
rz(2.9827859) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5903871) q[0];
sx q[0];
rz(-0.90999167) q[0];
sx q[0];
rz(1.4676276) q[0];
rz(0.42694636) q[2];
sx q[2];
rz(-1.5265577) q[2];
sx q[2];
rz(-2.6298863) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1013704) q[1];
sx q[1];
rz(-1.1565546) q[1];
sx q[1];
rz(-0.48419063) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9302904) q[3];
sx q[3];
rz(-1.4332472) q[3];
sx q[3];
rz(0.40776238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.020393546) q[2];
sx q[2];
rz(-0.32200107) q[2];
sx q[2];
rz(-0.89312345) q[2];
rz(2.761306) q[3];
sx q[3];
rz(-1.9392574) q[3];
sx q[3];
rz(-1.7072385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073008386) q[0];
sx q[0];
rz(-0.46982729) q[0];
sx q[0];
rz(-1.481886) q[0];
rz(-2.1549554) q[1];
sx q[1];
rz(-1.4886798) q[1];
sx q[1];
rz(-2.5843487) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38735861) q[0];
sx q[0];
rz(-3.0172303) q[0];
sx q[0];
rz(-2.6919305) q[0];
rz(-1.2115914) q[2];
sx q[2];
rz(-1.85382) q[2];
sx q[2];
rz(-2.46794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3448191) q[1];
sx q[1];
rz(-1.0386931) q[1];
sx q[1];
rz(2.565233) q[1];
rz(-pi) q[2];
rz(-2.557665) q[3];
sx q[3];
rz(-3.1213396) q[3];
sx q[3];
rz(-0.7113061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3096932) q[2];
sx q[2];
rz(-2.5688186) q[2];
sx q[2];
rz(-1.0581623) q[2];
rz(1.7582827) q[3];
sx q[3];
rz(-1.0388831) q[3];
sx q[3];
rz(-0.96341187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6365373) q[0];
sx q[0];
rz(-1.1543244) q[0];
sx q[0];
rz(-2.7099047) q[0];
rz(2.441326) q[1];
sx q[1];
rz(-1.2338748) q[1];
sx q[1];
rz(-2.4468927) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2496172) q[0];
sx q[0];
rz(-1.204899) q[0];
sx q[0];
rz(-2.1137733) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0236597) q[2];
sx q[2];
rz(-1.4768355) q[2];
sx q[2];
rz(0.31459034) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8164639) q[1];
sx q[1];
rz(-2.0184787) q[1];
sx q[1];
rz(-2.1486077) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9630345) q[3];
sx q[3];
rz(-1.6135912) q[3];
sx q[3];
rz(-0.6346441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6557287) q[2];
sx q[2];
rz(-2.8317917) q[2];
sx q[2];
rz(2.7244205) q[2];
rz(0.74472767) q[3];
sx q[3];
rz(-1.797902) q[3];
sx q[3];
rz(-2.1630796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3485296) q[0];
sx q[0];
rz(-1.2611669) q[0];
sx q[0];
rz(1.4766759) q[0];
rz(1.2186125) q[1];
sx q[1];
rz(-1.1398133) q[1];
sx q[1];
rz(-2.555991) q[1];
rz(-2.583859) q[2];
sx q[2];
rz(-0.9716059) q[2];
sx q[2];
rz(-1.5699408) q[2];
rz(2.811609) q[3];
sx q[3];
rz(-2.1654304) q[3];
sx q[3];
rz(2.5684857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
