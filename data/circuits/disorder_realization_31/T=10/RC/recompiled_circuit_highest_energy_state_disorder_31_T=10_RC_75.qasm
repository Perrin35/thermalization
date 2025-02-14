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
rz(1.0072768) q[0];
sx q[0];
rz(-0.64581031) q[0];
sx q[0];
rz(-2.7752152) q[0];
rz(0.74908787) q[1];
sx q[1];
rz(-1.5146511) q[1];
sx q[1];
rz(-0.14263022) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3245904) q[0];
sx q[0];
rz(-2.3329371) q[0];
sx q[0];
rz(2.3456744) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7103315) q[2];
sx q[2];
rz(-1.7931316) q[2];
sx q[2];
rz(2.2408533) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.2995404) q[1];
sx q[1];
rz(-1.6900542) q[1];
sx q[1];
rz(1.0930608) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47358114) q[3];
sx q[3];
rz(-1.8627121) q[3];
sx q[3];
rz(-1.1137258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35965219) q[2];
sx q[2];
rz(-0.47986978) q[2];
sx q[2];
rz(-1.5894319) q[2];
rz(-1.747067) q[3];
sx q[3];
rz(-2.7701869) q[3];
sx q[3];
rz(2.4594405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29779103) q[0];
sx q[0];
rz(-2.5318635) q[0];
sx q[0];
rz(-0.30526701) q[0];
rz(-1.8318532) q[1];
sx q[1];
rz(-2.7204308) q[1];
sx q[1];
rz(-3.0895244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6604495) q[0];
sx q[0];
rz(-1.7178806) q[0];
sx q[0];
rz(-1.4597465) q[0];
rz(-1.270625) q[2];
sx q[2];
rz(-1.9294293) q[2];
sx q[2];
rz(-1.7110362) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4365719) q[1];
sx q[1];
rz(-1.3051045) q[1];
sx q[1];
rz(-0.17036856) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0092485) q[3];
sx q[3];
rz(-2.2346046) q[3];
sx q[3];
rz(0.19451605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.084595844) q[2];
sx q[2];
rz(-1.8085542) q[2];
sx q[2];
rz(-1.7219273) q[2];
rz(-2.2325884) q[3];
sx q[3];
rz(-2.3125068) q[3];
sx q[3];
rz(-1.4066633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97515714) q[0];
sx q[0];
rz(-1.1844013) q[0];
sx q[0];
rz(-1.7433521) q[0];
rz(-0.94683975) q[1];
sx q[1];
rz(-1.0942752) q[1];
sx q[1];
rz(1.2907226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4104267) q[0];
sx q[0];
rz(-1.7058347) q[0];
sx q[0];
rz(0.076083994) q[0];
rz(-pi) q[1];
rz(0.47934808) q[2];
sx q[2];
rz(-0.84576037) q[2];
sx q[2];
rz(-0.59402639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7056928) q[1];
sx q[1];
rz(-1.5997643) q[1];
sx q[1];
rz(1.6375919) q[1];
rz(-pi) q[2];
rz(-2.6638899) q[3];
sx q[3];
rz(-1.3906428) q[3];
sx q[3];
rz(-2.2733062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0134086) q[2];
sx q[2];
rz(-0.42542294) q[2];
sx q[2];
rz(-2.181633) q[2];
rz(2.8175682) q[3];
sx q[3];
rz(-2.6257381) q[3];
sx q[3];
rz(-2.7014151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0078761) q[0];
sx q[0];
rz(-0.77819264) q[0];
sx q[0];
rz(-2.9414951) q[0];
rz(-2.5307185) q[1];
sx q[1];
rz(-2.4433544) q[1];
sx q[1];
rz(0.79455882) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3171697) q[0];
sx q[0];
rz(-2.9545435) q[0];
sx q[0];
rz(-0.6424128) q[0];
rz(-0.17996338) q[2];
sx q[2];
rz(-0.90562253) q[2];
sx q[2];
rz(-2.6106204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17758372) q[1];
sx q[1];
rz(-1.4191886) q[1];
sx q[1];
rz(-1.8880009) q[1];
x q[2];
rz(0.67312397) q[3];
sx q[3];
rz(-0.84852058) q[3];
sx q[3];
rz(2.9150042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8841298) q[2];
sx q[2];
rz(-0.56371671) q[2];
sx q[2];
rz(1.5719315) q[2];
rz(1.3563159) q[3];
sx q[3];
rz(-1.1607728) q[3];
sx q[3];
rz(-0.13264382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86303478) q[0];
sx q[0];
rz(-0.25535169) q[0];
sx q[0];
rz(2.417946) q[0];
rz(3.0408995) q[1];
sx q[1];
rz(-2.0932525) q[1];
sx q[1];
rz(0.06631276) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50936705) q[0];
sx q[0];
rz(-2.5697587) q[0];
sx q[0];
rz(-1.454005) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91513779) q[2];
sx q[2];
rz(-1.1205289) q[2];
sx q[2];
rz(1.7129667) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4206754) q[1];
sx q[1];
rz(-2.2128792) q[1];
sx q[1];
rz(-2.9319619) q[1];
rz(-2.4661786) q[3];
sx q[3];
rz(-2.6468122) q[3];
sx q[3];
rz(0.29008415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2012607) q[2];
sx q[2];
rz(-1.9278434) q[2];
sx q[2];
rz(0.21855375) q[2];
rz(-0.39441937) q[3];
sx q[3];
rz(-0.68587488) q[3];
sx q[3];
rz(0.39943892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047693096) q[0];
sx q[0];
rz(-0.91966367) q[0];
sx q[0];
rz(0.0041051824) q[0];
rz(-0.63402367) q[1];
sx q[1];
rz(-1.6396294) q[1];
sx q[1];
rz(-2.1852469) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078236899) q[0];
sx q[0];
rz(-0.18785297) q[0];
sx q[0];
rz(-0.22042033) q[0];
rz(-pi) q[1];
rz(1.1460828) q[2];
sx q[2];
rz(-2.0419952) q[2];
sx q[2];
rz(1.4482519) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9117212) q[1];
sx q[1];
rz(-1.2930065) q[1];
sx q[1];
rz(-1.576527) q[1];
rz(0.34849747) q[3];
sx q[3];
rz(-1.6808884) q[3];
sx q[3];
rz(-1.7736721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52028209) q[2];
sx q[2];
rz(-1.4352398) q[2];
sx q[2];
rz(0.83462805) q[2];
rz(-0.37072119) q[3];
sx q[3];
rz(-0.70086896) q[3];
sx q[3];
rz(-2.4247775) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9653559) q[0];
sx q[0];
rz(-3.1270202) q[0];
sx q[0];
rz(-2.6915349) q[0];
rz(-0.4484446) q[1];
sx q[1];
rz(-0.83201718) q[1];
sx q[1];
rz(0.73743302) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5784604) q[0];
sx q[0];
rz(-1.9942584) q[0];
sx q[0];
rz(2.9388301) q[0];
rz(1.4294523) q[2];
sx q[2];
rz(-1.6108043) q[2];
sx q[2];
rz(1.227898) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0417252) q[1];
sx q[1];
rz(-1.2507572) q[1];
sx q[1];
rz(2.7766224) q[1];
rz(-2.8705863) q[3];
sx q[3];
rz(-0.79870634) q[3];
sx q[3];
rz(-2.9661111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8488778) q[2];
sx q[2];
rz(-2.541031) q[2];
sx q[2];
rz(-2.4265477) q[2];
rz(-0.50825608) q[3];
sx q[3];
rz(-1.6786989) q[3];
sx q[3];
rz(1.3983294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7364863) q[0];
sx q[0];
rz(-2.5953601) q[0];
sx q[0];
rz(2.7832094) q[0];
rz(-0.37250039) q[1];
sx q[1];
rz(-1.7820216) q[1];
sx q[1];
rz(0.43982664) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7067703) q[0];
sx q[0];
rz(-2.3447038) q[0];
sx q[0];
rz(-0.057698542) q[0];
rz(-pi) q[1];
x q[1];
rz(1.110564) q[2];
sx q[2];
rz(-1.2046736) q[2];
sx q[2];
rz(1.5535959) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3258626) q[1];
sx q[1];
rz(-1.7472032) q[1];
sx q[1];
rz(2.8730064) q[1];
rz(-pi) q[2];
rz(1.5164916) q[3];
sx q[3];
rz(-1.2479932) q[3];
sx q[3];
rz(-1.0564547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0899352) q[2];
sx q[2];
rz(-0.097044162) q[2];
sx q[2];
rz(0.70866054) q[2];
rz(0.25157252) q[3];
sx q[3];
rz(-2.1422062) q[3];
sx q[3];
rz(-0.034041762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2947023) q[0];
sx q[0];
rz(-2.1069694) q[0];
sx q[0];
rz(-2.5501472) q[0];
rz(-0.017631831) q[1];
sx q[1];
rz(-1.0523187) q[1];
sx q[1];
rz(-0.10094053) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.712698) q[0];
sx q[0];
rz(-0.34358812) q[0];
sx q[0];
rz(-1.982292) q[0];
rz(-pi) q[1];
rz(3.0885126) q[2];
sx q[2];
rz(-2.5042222) q[2];
sx q[2];
rz(-0.87634477) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7314607) q[1];
sx q[1];
rz(-1.4515775) q[1];
sx q[1];
rz(1.340926) q[1];
x q[2];
rz(-1.798257) q[3];
sx q[3];
rz(-0.69256562) q[3];
sx q[3];
rz(-0.3161968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4883604) q[2];
sx q[2];
rz(-2.4182352) q[2];
sx q[2];
rz(-1.6548659) q[2];
rz(-0.15445736) q[3];
sx q[3];
rz(-1.1052701) q[3];
sx q[3];
rz(-0.36623335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79009295) q[0];
sx q[0];
rz(-0.14227754) q[0];
sx q[0];
rz(-2.067814) q[0];
rz(-1.0723266) q[1];
sx q[1];
rz(-2.8876979) q[1];
sx q[1];
rz(0.059018746) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.430141) q[0];
sx q[0];
rz(-1.3716328) q[0];
sx q[0];
rz(1.5368622) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8779936) q[2];
sx q[2];
rz(-1.0754183) q[2];
sx q[2];
rz(1.6558629) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5328588) q[1];
sx q[1];
rz(-2.1648832) q[1];
sx q[1];
rz(0.913175) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2667472) q[3];
sx q[3];
rz(-2.2922462) q[3];
sx q[3];
rz(1.9375999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8129639) q[2];
sx q[2];
rz(-1.2920516) q[2];
sx q[2];
rz(-2.9126419) q[2];
rz(-1.9479343) q[3];
sx q[3];
rz(-0.3242068) q[3];
sx q[3];
rz(1.4540023) q[3];
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
rz(3.1345632) q[0];
sx q[0];
rz(-0.80637359) q[0];
sx q[0];
rz(-0.73643186) q[0];
rz(-1.7560584) q[1];
sx q[1];
rz(-1.7693188) q[1];
sx q[1];
rz(1.7973695) q[1];
rz(-0.028111721) q[2];
sx q[2];
rz(-2.1098469) q[2];
sx q[2];
rz(-2.5218365) q[2];
rz(1.368461) q[3];
sx q[3];
rz(-0.96755316) q[3];
sx q[3];
rz(0.62300976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
