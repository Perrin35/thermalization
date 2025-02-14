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
rz(-2.0237145) q[0];
sx q[0];
rz(-1.0853381) q[0];
sx q[0];
rz(-0.35882741) q[0];
rz(1.2598414) q[1];
sx q[1];
rz(-1.0955732) q[1];
sx q[1];
rz(2.1020558) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0721999) q[0];
sx q[0];
rz(-0.93749267) q[0];
sx q[0];
rz(1.1266438) q[0];
rz(0.02746901) q[2];
sx q[2];
rz(-1.6131764) q[2];
sx q[2];
rz(-0.24951631) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1778587) q[1];
sx q[1];
rz(-1.2716645) q[1];
sx q[1];
rz(-1.9480716) q[1];
rz(0.12173481) q[3];
sx q[3];
rz(-2.4678708) q[3];
sx q[3];
rz(-2.664543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.871668) q[2];
sx q[2];
rz(-2.2871127) q[2];
sx q[2];
rz(0.16647767) q[2];
rz(0.13655937) q[3];
sx q[3];
rz(-2.2763493) q[3];
sx q[3];
rz(-1.4599919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0110375) q[0];
sx q[0];
rz(-2.8189973) q[0];
sx q[0];
rz(0.36001298) q[0];
rz(-0.75045466) q[1];
sx q[1];
rz(-2.2513335) q[1];
sx q[1];
rz(-3.0976565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1587857) q[0];
sx q[0];
rz(-1.7097397) q[0];
sx q[0];
rz(1.7278633) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0925748) q[2];
sx q[2];
rz(-0.92713461) q[2];
sx q[2];
rz(1.1903534) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3219731) q[1];
sx q[1];
rz(-2.3016046) q[1];
sx q[1];
rz(-1.8614443) q[1];
rz(-pi) q[2];
rz(2.206541) q[3];
sx q[3];
rz(-2.8160444) q[3];
sx q[3];
rz(2.9071484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3540078) q[2];
sx q[2];
rz(-2.6968991) q[2];
sx q[2];
rz(1.0178052) q[2];
rz(-0.20778896) q[3];
sx q[3];
rz(-2.5354247) q[3];
sx q[3];
rz(-0.35473216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6924234) q[0];
sx q[0];
rz(-2.848931) q[0];
sx q[0];
rz(-1.0798651) q[0];
rz(1.0049459) q[1];
sx q[1];
rz(-2.4506073) q[1];
sx q[1];
rz(1.2452004) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5737138) q[0];
sx q[0];
rz(-1.880293) q[0];
sx q[0];
rz(-1.482016) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7342308) q[2];
sx q[2];
rz(-1.15112) q[2];
sx q[2];
rz(-0.070022665) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.019466467) q[1];
sx q[1];
rz(-1.2858741) q[1];
sx q[1];
rz(-1.0128922) q[1];
rz(-pi) q[2];
rz(1.336634) q[3];
sx q[3];
rz(-2.5076205) q[3];
sx q[3];
rz(2.0869107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9006742) q[2];
sx q[2];
rz(-1.2718879) q[2];
sx q[2];
rz(0.34350485) q[2];
rz(-3.1268696) q[3];
sx q[3];
rz(-0.69535178) q[3];
sx q[3];
rz(1.2408313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(3.0111276) q[0];
sx q[0];
rz(-2.0978329) q[0];
sx q[0];
rz(0.68487942) q[0];
rz(-0.20826805) q[1];
sx q[1];
rz(-1.8332278) q[1];
sx q[1];
rz(2.3484255) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7670534) q[0];
sx q[0];
rz(-1.6691795) q[0];
sx q[0];
rz(-1.0258425) q[0];
rz(2.9689413) q[2];
sx q[2];
rz(-1.358629) q[2];
sx q[2];
rz(1.4868975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4186331) q[1];
sx q[1];
rz(-1.5546989) q[1];
sx q[1];
rz(2.1820585) q[1];
rz(-2.0871619) q[3];
sx q[3];
rz(-2.722009) q[3];
sx q[3];
rz(-0.53158599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3785279) q[2];
sx q[2];
rz(-1.6918162) q[2];
sx q[2];
rz(-2.6587963) q[2];
rz(-1.637508) q[3];
sx q[3];
rz(-0.62862325) q[3];
sx q[3];
rz(-0.32046902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046431635) q[0];
sx q[0];
rz(-0.23155364) q[0];
sx q[0];
rz(-0.10145536) q[0];
rz(-0.21508148) q[1];
sx q[1];
rz(-1.2098034) q[1];
sx q[1];
rz(-1.2717517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47574193) q[0];
sx q[0];
rz(-1.8545194) q[0];
sx q[0];
rz(-2.0497149) q[0];
rz(-pi) q[1];
rz(-0.59753363) q[2];
sx q[2];
rz(-1.1206822) q[2];
sx q[2];
rz(2.2240153) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7600619) q[1];
sx q[1];
rz(-0.4222798) q[1];
sx q[1];
rz(1.1997919) q[1];
rz(-pi) q[2];
rz(-2.1121096) q[3];
sx q[3];
rz(-2.7446163) q[3];
sx q[3];
rz(3.027806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35470206) q[2];
sx q[2];
rz(-0.44822794) q[2];
sx q[2];
rz(-2.9917713) q[2];
rz(0.24770728) q[3];
sx q[3];
rz(-0.37023538) q[3];
sx q[3];
rz(-2.3410334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.596452) q[0];
sx q[0];
rz(-2.5189724) q[0];
sx q[0];
rz(-1.2860292) q[0];
rz(2.3226358) q[1];
sx q[1];
rz(-0.30636925) q[1];
sx q[1];
rz(2.9584296) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.055119) q[0];
sx q[0];
rz(-1.3684784) q[0];
sx q[0];
rz(-2.5760465) q[0];
rz(-pi) q[1];
rz(-2.9565565) q[2];
sx q[2];
rz(-1.5079947) q[2];
sx q[2];
rz(-0.35279314) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0135063) q[1];
sx q[1];
rz(-1.9486843) q[1];
sx q[1];
rz(-2.3289447) q[1];
rz(-pi) q[2];
rz(-2.6425903) q[3];
sx q[3];
rz(-1.2899644) q[3];
sx q[3];
rz(-2.4256191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5861107) q[2];
sx q[2];
rz(-1.7038944) q[2];
sx q[2];
rz(-2.6044031) q[2];
rz(-1.7389343) q[3];
sx q[3];
rz(-1.2638998) q[3];
sx q[3];
rz(2.5130443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1527767) q[0];
sx q[0];
rz(-0.39867908) q[0];
sx q[0];
rz(-3.0231349) q[0];
rz(1.1064103) q[1];
sx q[1];
rz(-1.9170008) q[1];
sx q[1];
rz(-0.67682636) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7962435) q[0];
sx q[0];
rz(-1.5156094) q[0];
sx q[0];
rz(0.075166239) q[0];
rz(0.17497988) q[2];
sx q[2];
rz(-0.6616486) q[2];
sx q[2];
rz(-1.511919) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.54296498) q[1];
sx q[1];
rz(-1.7522546) q[1];
sx q[1];
rz(-1.7398029) q[1];
rz(-1.5561701) q[3];
sx q[3];
rz(-1.5717603) q[3];
sx q[3];
rz(0.87982501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7532838) q[2];
sx q[2];
rz(-2.1948185) q[2];
sx q[2];
rz(3.0906265) q[2];
rz(-2.986159) q[3];
sx q[3];
rz(-1.6499465) q[3];
sx q[3];
rz(-1.0203993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8841356) q[0];
sx q[0];
rz(-1.1549042) q[0];
sx q[0];
rz(0.98404032) q[0];
rz(2.6718196) q[1];
sx q[1];
rz(-1.677294) q[1];
sx q[1];
rz(-0.22107548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8399682) q[0];
sx q[0];
rz(-1.589612) q[0];
sx q[0];
rz(-3.0473659) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65024489) q[2];
sx q[2];
rz(-1.2653554) q[2];
sx q[2];
rz(0.98750706) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.03121491) q[1];
sx q[1];
rz(-1.7965914) q[1];
sx q[1];
rz(1.8757191) q[1];
rz(-pi) q[2];
rz(0.89721591) q[3];
sx q[3];
rz(-0.69475896) q[3];
sx q[3];
rz(-1.0396797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2880963) q[2];
sx q[2];
rz(-2.4595538) q[2];
sx q[2];
rz(1.811553) q[2];
rz(2.5174777) q[3];
sx q[3];
rz(-0.95804405) q[3];
sx q[3];
rz(-0.38930711) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57988155) q[0];
sx q[0];
rz(-1.6498673) q[0];
sx q[0];
rz(1.0505744) q[0];
rz(0.9802649) q[1];
sx q[1];
rz(-1.0824208) q[1];
sx q[1];
rz(-1.6260653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4517364) q[0];
sx q[0];
rz(-2.0961267) q[0];
sx q[0];
rz(-1.3888874) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7057436) q[2];
sx q[2];
rz(-0.51700355) q[2];
sx q[2];
rz(-0.85384254) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84534553) q[1];
sx q[1];
rz(-1.913684) q[1];
sx q[1];
rz(1.898349) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6330746) q[3];
sx q[3];
rz(-2.4313723) q[3];
sx q[3];
rz(2.4274277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3773697) q[2];
sx q[2];
rz(-1.919596) q[2];
sx q[2];
rz(0.73966217) q[2];
rz(2.9861279) q[3];
sx q[3];
rz(-0.37853095) q[3];
sx q[3];
rz(-0.19708656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4231606) q[0];
sx q[0];
rz(-2.6078556) q[0];
sx q[0];
rz(2.0045795) q[0];
rz(0.63277376) q[1];
sx q[1];
rz(-2.4686765) q[1];
sx q[1];
rz(0.64579642) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5006913) q[0];
sx q[0];
rz(-2.7703022) q[0];
sx q[0];
rz(-2.6616606) q[0];
rz(-pi) q[1];
rz(2.907173) q[2];
sx q[2];
rz(-0.79371292) q[2];
sx q[2];
rz(-0.91784263) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91612047) q[1];
sx q[1];
rz(-1.6399621) q[1];
sx q[1];
rz(-2.3232077) q[1];
x q[2];
rz(-2.0955576) q[3];
sx q[3];
rz(-1.2933795) q[3];
sx q[3];
rz(0.21464561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48468963) q[2];
sx q[2];
rz(-0.36269665) q[2];
sx q[2];
rz(-3.1032069) q[2];
rz(-0.11135993) q[3];
sx q[3];
rz(-2.7146118) q[3];
sx q[3];
rz(-2.4666983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69042324) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(1.3896598) q[0];
rz(1.2870862) q[1];
sx q[1];
rz(-0.99936395) q[1];
sx q[1];
rz(-2.0157464) q[1];
rz(-0.95697516) q[2];
sx q[2];
rz(-2.0718832) q[2];
sx q[2];
rz(2.8399443) q[2];
rz(1.7467563) q[3];
sx q[3];
rz(-0.8192438) q[3];
sx q[3];
rz(1.4024709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
