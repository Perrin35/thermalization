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
rz(1.6588563) q[0];
sx q[0];
rz(-0.98134494) q[0];
sx q[0];
rz(1.097783) q[0];
rz(-0.78805796) q[1];
sx q[1];
rz(-2.0907953) q[1];
sx q[1];
rz(3.1346336) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5612457) q[0];
sx q[0];
rz(-0.92134199) q[0];
sx q[0];
rz(0.070493079) q[0];
rz(-pi) q[1];
rz(0.2351951) q[2];
sx q[2];
rz(-0.93738745) q[2];
sx q[2];
rz(-0.051284479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4930058) q[1];
sx q[1];
rz(-1.3351591) q[1];
sx q[1];
rz(3.1070903) q[1];
x q[2];
rz(-1.4640635) q[3];
sx q[3];
rz(-1.2327223) q[3];
sx q[3];
rz(2.6298912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2385345) q[2];
sx q[2];
rz(-1.3471341) q[2];
sx q[2];
rz(-1.6713589) q[2];
rz(0.12601958) q[3];
sx q[3];
rz(-2.6991548) q[3];
sx q[3];
rz(0.68388763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11494342) q[0];
sx q[0];
rz(-2.6925955) q[0];
sx q[0];
rz(-2.5057416) q[0];
rz(2.3682829) q[1];
sx q[1];
rz(-0.67716235) q[1];
sx q[1];
rz(-1.4453452) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6757386) q[0];
sx q[0];
rz(-0.93198085) q[0];
sx q[0];
rz(-0.56446989) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0436297) q[2];
sx q[2];
rz(-0.84542984) q[2];
sx q[2];
rz(0.54976094) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8938039) q[1];
sx q[1];
rz(-0.87444982) q[1];
sx q[1];
rz(1.3651834) q[1];
x q[2];
rz(-1.7023193) q[3];
sx q[3];
rz(-1.6386541) q[3];
sx q[3];
rz(-1.1179641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.14629743) q[2];
sx q[2];
rz(-1.7701021) q[2];
sx q[2];
rz(0.13776097) q[2];
rz(2.6855101) q[3];
sx q[3];
rz(-0.59499732) q[3];
sx q[3];
rz(2.6661787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59738961) q[0];
sx q[0];
rz(-1.5578288) q[0];
sx q[0];
rz(-0.57918817) q[0];
rz(-0.10313615) q[1];
sx q[1];
rz(-1.3214279) q[1];
sx q[1];
rz(-0.78027049) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73585498) q[0];
sx q[0];
rz(-3.0445703) q[0];
sx q[0];
rz(0.19543044) q[0];
x q[1];
rz(-2.4603526) q[2];
sx q[2];
rz(-1.4461755) q[2];
sx q[2];
rz(-1.7523426) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6578778) q[1];
sx q[1];
rz(-0.45770744) q[1];
sx q[1];
rz(1.0545516) q[1];
rz(2.2663651) q[3];
sx q[3];
rz(-1.4000633) q[3];
sx q[3];
rz(3.1306992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66037336) q[2];
sx q[2];
rz(-1.4619091) q[2];
sx q[2];
rz(-0.090864651) q[2];
rz(-1.4800492) q[3];
sx q[3];
rz(-1.816498) q[3];
sx q[3];
rz(-2.3265694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12260967) q[0];
sx q[0];
rz(-2.9538587) q[0];
sx q[0];
rz(-2.6222099) q[0];
rz(-1.9620365) q[1];
sx q[1];
rz(-0.68125454) q[1];
sx q[1];
rz(0.048351668) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78084842) q[0];
sx q[0];
rz(-2.0791302) q[0];
sx q[0];
rz(-2.2717649) q[0];
rz(0.91421336) q[2];
sx q[2];
rz(-2.5541452) q[2];
sx q[2];
rz(0.32876164) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.4614233) q[1];
sx q[1];
rz(-1.5986414) q[1];
sx q[1];
rz(0.48770406) q[1];
rz(-pi) q[2];
rz(0.79508852) q[3];
sx q[3];
rz(-1.3833191) q[3];
sx q[3];
rz(0.80611967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68253303) q[2];
sx q[2];
rz(-0.30632633) q[2];
sx q[2];
rz(-1.691386) q[2];
rz(-2.0356483) q[3];
sx q[3];
rz(-1.148843) q[3];
sx q[3];
rz(-2.2753184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3806216) q[0];
sx q[0];
rz(-2.3411317) q[0];
sx q[0];
rz(-0.77734429) q[0];
rz(0.49579534) q[1];
sx q[1];
rz(-0.40758857) q[1];
sx q[1];
rz(0.19283238) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042943311) q[0];
sx q[0];
rz(-2.1060364) q[0];
sx q[0];
rz(-1.9283867) q[0];
rz(-2.8933346) q[2];
sx q[2];
rz(-0.78310668) q[2];
sx q[2];
rz(-3.1241036) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3461756) q[1];
sx q[1];
rz(-0.3095135) q[1];
sx q[1];
rz(-1.2796106) q[1];
x q[2];
rz(1.2019346) q[3];
sx q[3];
rz(-0.4274803) q[3];
sx q[3];
rz(-0.99586801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5567646) q[2];
sx q[2];
rz(-0.77797055) q[2];
sx q[2];
rz(0.83621109) q[2];
rz(-0.95544514) q[3];
sx q[3];
rz(-1.7183869) q[3];
sx q[3];
rz(0.53171617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.29109508) q[0];
sx q[0];
rz(-1.9170772) q[0];
sx q[0];
rz(0.033705458) q[0];
rz(-0.91824245) q[1];
sx q[1];
rz(-2.4372209) q[1];
sx q[1];
rz(-2.4423626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540772) q[0];
sx q[0];
rz(-1.0654385) q[0];
sx q[0];
rz(2.2127241) q[0];
x q[1];
rz(2.9543058) q[2];
sx q[2];
rz(-1.5178871) q[2];
sx q[2];
rz(-2.1250181) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.110975) q[1];
sx q[1];
rz(-0.4781107) q[1];
sx q[1];
rz(0.39076372) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27856234) q[3];
sx q[3];
rz(-1.8980366) q[3];
sx q[3];
rz(-0.95765169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.08192) q[2];
sx q[2];
rz(-0.6531738) q[2];
sx q[2];
rz(0.095452249) q[2];
rz(-2.0898315) q[3];
sx q[3];
rz(-1.8973408) q[3];
sx q[3];
rz(-2.2473647) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8555701) q[0];
sx q[0];
rz(-0.48208553) q[0];
sx q[0];
rz(1.1676189) q[0];
rz(-2.7883912) q[1];
sx q[1];
rz(-2.628852) q[1];
sx q[1];
rz(2.8599427) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7324556) q[0];
sx q[0];
rz(-2.0444217) q[0];
sx q[0];
rz(0.098129674) q[0];
rz(3.0007576) q[2];
sx q[2];
rz(-0.65208921) q[2];
sx q[2];
rz(-1.2229133) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6503164) q[1];
sx q[1];
rz(-1.7185083) q[1];
sx q[1];
rz(0.13593918) q[1];
x q[2];
rz(1.0135039) q[3];
sx q[3];
rz(-1.9907328) q[3];
sx q[3];
rz(2.1518663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1436651) q[2];
sx q[2];
rz(-1.2219656) q[2];
sx q[2];
rz(1.3227051) q[2];
rz(1.6392684) q[3];
sx q[3];
rz(-2.0570677) q[3];
sx q[3];
rz(3.0433906) q[3];
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
rz(-2.1439576) q[0];
sx q[0];
rz(-1.2459545) q[0];
sx q[0];
rz(-2.8676046) q[0];
rz(-2.5471089) q[1];
sx q[1];
rz(-1.7285873) q[1];
sx q[1];
rz(-0.53057539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45121058) q[0];
sx q[0];
rz(-1.5669826) q[0];
sx q[0];
rz(0.03937748) q[0];
x q[1];
rz(0.12315788) q[2];
sx q[2];
rz(-2.5854857) q[2];
sx q[2];
rz(1.3431988) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1889607) q[1];
sx q[1];
rz(-0.41302339) q[1];
sx q[1];
rz(-2.7160591) q[1];
rz(-pi) q[2];
x q[2];
rz(1.382393) q[3];
sx q[3];
rz(-2.419988) q[3];
sx q[3];
rz(2.5327275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8990367) q[2];
sx q[2];
rz(-1.2071004) q[2];
sx q[2];
rz(0.25699082) q[2];
rz(-1.9422003) q[3];
sx q[3];
rz(-1.2446087) q[3];
sx q[3];
rz(1.6283584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.5892107) q[0];
sx q[0];
rz(-2.1871545) q[0];
sx q[0];
rz(2.6115665) q[0];
rz(-1.5026622) q[1];
sx q[1];
rz(-1.9866147) q[1];
sx q[1];
rz(-2.2030305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5230704) q[0];
sx q[0];
rz(-2.0393199) q[0];
sx q[0];
rz(0.65798385) q[0];
x q[1];
rz(-2.4468719) q[2];
sx q[2];
rz(-0.18604569) q[2];
sx q[2];
rz(-0.11878601) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5512498) q[1];
sx q[1];
rz(-1.9971202) q[1];
sx q[1];
rz(0.23386441) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0051232) q[3];
sx q[3];
rz(-2.1519063) q[3];
sx q[3];
rz(-2.7032397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75371257) q[2];
sx q[2];
rz(-2.2215999) q[2];
sx q[2];
rz(-2.06125) q[2];
rz(-2.3267817) q[3];
sx q[3];
rz(-1.9459414) q[3];
sx q[3];
rz(-1.7241534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1431047) q[0];
sx q[0];
rz(-1.4838706) q[0];
sx q[0];
rz(1.0888354) q[0];
rz(3.0740652) q[1];
sx q[1];
rz(-1.2779526) q[1];
sx q[1];
rz(0.75928226) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7831206) q[0];
sx q[0];
rz(-1.4519339) q[0];
sx q[0];
rz(-1.8274183) q[0];
rz(-0.49726059) q[2];
sx q[2];
rz(-2.8858375) q[2];
sx q[2];
rz(-1.349468) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2226505) q[1];
sx q[1];
rz(-0.83178751) q[1];
sx q[1];
rz(1.2779425) q[1];
rz(-pi) q[2];
rz(-0.57709007) q[3];
sx q[3];
rz(-0.40501696) q[3];
sx q[3];
rz(1.4228203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8074983) q[2];
sx q[2];
rz(-1.0958025) q[2];
sx q[2];
rz(-2.5661772) q[2];
rz(2.5790162) q[3];
sx q[3];
rz(-0.96499363) q[3];
sx q[3];
rz(-1.3728728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0811049) q[0];
sx q[0];
rz(-1.8555547) q[0];
sx q[0];
rz(-0.9077358) q[0];
rz(1.0870712) q[1];
sx q[1];
rz(-1.1320976) q[1];
sx q[1];
rz(-0.017398106) q[1];
rz(0.74093735) q[2];
sx q[2];
rz(-2.9334684) q[2];
sx q[2];
rz(2.8192782) q[2];
rz(-1.9418874) q[3];
sx q[3];
rz(-2.1368847) q[3];
sx q[3];
rz(0.18659244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
