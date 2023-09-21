OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0032089631) q[0];
sx q[0];
rz(-0.15455833) q[0];
sx q[0];
rz(0.69252339) q[0];
rz(5.0737557) q[1];
sx q[1];
rz(4.3901246) q[1];
sx q[1];
rz(7.6683383) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9557578) q[0];
sx q[0];
rz(-1.251936) q[0];
sx q[0];
rz(2.3715109) q[0];
rz(-pi) q[1];
rz(2.5897883) q[2];
sx q[2];
rz(-1.8066415) q[2];
sx q[2];
rz(0.15197309) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9618122) q[1];
sx q[1];
rz(-1.4084067) q[1];
sx q[1];
rz(-0.23602545) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0632283) q[3];
sx q[3];
rz(-0.97798097) q[3];
sx q[3];
rz(1.4584695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2549071) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(-0.20516667) q[2];
rz(2.3702879) q[3];
sx q[3];
rz(-2.3588534) q[3];
sx q[3];
rz(-1.1024968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40760621) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(-2.6876887) q[0];
rz(-2.1167963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(-1.9143547) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95943806) q[0];
sx q[0];
rz(-1.6998788) q[0];
sx q[0];
rz(0.076230787) q[0];
rz(-pi) q[1];
rz(-0.71364673) q[2];
sx q[2];
rz(-2.4948641) q[2];
sx q[2];
rz(1.5734067) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7770575) q[1];
sx q[1];
rz(-2.747066) q[1];
sx q[1];
rz(-2.415034) q[1];
rz(-0.65900461) q[3];
sx q[3];
rz(-2.8895624) q[3];
sx q[3];
rz(-0.27740955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1002905) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(2.5741637) q[2];
rz(2.7764017) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(2.173483) q[3];
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
rz(-2.658618) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(-2.2429402) q[0];
rz(-2.1458416) q[1];
sx q[1];
rz(-1.5834705) q[1];
sx q[1];
rz(-0.333289) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8329187) q[0];
sx q[0];
rz(-2.0325639) q[0];
sx q[0];
rz(-2.4426016) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7477481) q[2];
sx q[2];
rz(-1.4743917) q[2];
sx q[2];
rz(-2.5071438) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.052735141) q[1];
sx q[1];
rz(-0.6328859) q[1];
sx q[1];
rz(1.8295374) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7379052) q[3];
sx q[3];
rz(-1.0567769) q[3];
sx q[3];
rz(-0.46883632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68625346) q[2];
sx q[2];
rz(-1.3505961) q[2];
sx q[2];
rz(-1.1509482) q[2];
rz(2.3006556) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(1.9870728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44160098) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(-2.4568795) q[0];
rz(-2.1060064) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(1.205014) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9158463) q[0];
sx q[0];
rz(-1.5374743) q[0];
sx q[0];
rz(0.87945625) q[0];
x q[1];
rz(0.20385216) q[2];
sx q[2];
rz(-0.9365558) q[2];
sx q[2];
rz(-0.39433345) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.64441427) q[1];
sx q[1];
rz(-2.2481611) q[1];
sx q[1];
rz(-0.35269423) q[1];
x q[2];
rz(-1.6550001) q[3];
sx q[3];
rz(-2.4349672) q[3];
sx q[3];
rz(0.018761793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.24923199) q[2];
sx q[2];
rz(-1.7148596) q[2];
sx q[2];
rz(0.37115804) q[2];
rz(-1.7403729) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(1.1192809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-3.086833) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(3.0084685) q[0];
rz(0.99331028) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(-2.5865119) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9585882) q[0];
sx q[0];
rz(-2.8259813) q[0];
sx q[0];
rz(1.5297699) q[0];
rz(-1.416989) q[2];
sx q[2];
rz(-2.7222735) q[2];
sx q[2];
rz(-2.9749982) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.811378) q[1];
sx q[1];
rz(-2.0159617) q[1];
sx q[1];
rz(0.021275714) q[1];
rz(-pi) q[2];
rz(-3.0606907) q[3];
sx q[3];
rz(-1.5187851) q[3];
sx q[3];
rz(0.73976433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.83539) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(0.13892697) q[2];
rz(-2.1991918) q[3];
sx q[3];
rz(-1.5709632) q[3];
sx q[3];
rz(0.55148235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979364) q[0];
sx q[0];
rz(-2.5718226) q[0];
sx q[0];
rz(0.57975769) q[0];
rz(3.014091) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(-1.6019843) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.416076) q[0];
sx q[0];
rz(-0.86891876) q[0];
sx q[0];
rz(0.26728018) q[0];
rz(-pi) q[1];
rz(1.5251665) q[2];
sx q[2];
rz(-1.7012193) q[2];
sx q[2];
rz(-1.4554731) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.786799) q[1];
sx q[1];
rz(-1.2878294) q[1];
sx q[1];
rz(1.1250886) q[1];
x q[2];
rz(-1.9220011) q[3];
sx q[3];
rz(-0.71577365) q[3];
sx q[3];
rz(0.90482611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55398983) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(2.8721151) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(-0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44678974) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(-2.457298) q[0];
rz(0.11958312) q[1];
sx q[1];
rz(-1.8493098) q[1];
sx q[1];
rz(-0.51876846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6249734) q[0];
sx q[0];
rz(-0.74665753) q[0];
sx q[0];
rz(-1.3730511) q[0];
x q[1];
rz(1.7905551) q[2];
sx q[2];
rz(-2.3292654) q[2];
sx q[2];
rz(-1.0202368) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7873951) q[1];
sx q[1];
rz(-2.0320315) q[1];
sx q[1];
rz(-0.72871491) q[1];
rz(-pi) q[2];
rz(1.9714868) q[3];
sx q[3];
rz(-1.5466994) q[3];
sx q[3];
rz(-1.8322242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2542904) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(-0.56345144) q[2];
rz(-0.051579483) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(-2.1896867) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96034399) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(-1.3990336) q[0];
rz(-0.78701204) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(0.74434892) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9104011) q[0];
sx q[0];
rz(-1.4230799) q[0];
sx q[0];
rz(-1.6258679) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8366351) q[2];
sx q[2];
rz(-2.611534) q[2];
sx q[2];
rz(2.838138) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5163755) q[1];
sx q[1];
rz(-2.5853734) q[1];
sx q[1];
rz(-0.43362995) q[1];
rz(2.6870071) q[3];
sx q[3];
rz(-1.5919519) q[3];
sx q[3];
rz(1.3190312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7156334) q[2];
sx q[2];
rz(-1.3954433) q[2];
sx q[2];
rz(-2.5320833) q[2];
rz(2.4842747) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(-2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1881926) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(1.6375861) q[0];
rz(1.9001182) q[1];
sx q[1];
rz(-1.991661) q[1];
sx q[1];
rz(0.77493587) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0477714) q[0];
sx q[0];
rz(-1.1822961) q[0];
sx q[0];
rz(2.2872778) q[0];
rz(-2.0140892) q[2];
sx q[2];
rz(-1.7121676) q[2];
sx q[2];
rz(3.1183426) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7271125) q[1];
sx q[1];
rz(-2.6467122) q[1];
sx q[1];
rz(-0.40463573) q[1];
rz(2.1433467) q[3];
sx q[3];
rz(-0.25902723) q[3];
sx q[3];
rz(-1.6365901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9541786) q[2];
sx q[2];
rz(-0.22970197) q[2];
sx q[2];
rz(2.9837218) q[2];
rz(1.212451) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(1.8635748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.071844) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(2.9272595) q[0];
rz(0.65746039) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(2.0956031) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56453088) q[0];
sx q[0];
rz(-2.7663167) q[0];
sx q[0];
rz(2.8481086) q[0];
x q[1];
rz(0.043141456) q[2];
sx q[2];
rz(-0.59141814) q[2];
sx q[2];
rz(-2.6864348) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6007081) q[1];
sx q[1];
rz(-0.67544671) q[1];
sx q[1];
rz(2.1727968) q[1];
x q[2];
rz(3.0797327) q[3];
sx q[3];
rz(-2.7354771) q[3];
sx q[3];
rz(-1.5861685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7252698) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(-1.0894758) q[2];
rz(-1.5661092) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(-2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2789223) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(-1.5325585) q[1];
sx q[1];
rz(-1.4668203) q[1];
sx q[1];
rz(-1.1046881) q[1];
rz(2.7325148) q[2];
sx q[2];
rz(-1.8553875) q[2];
sx q[2];
rz(2.6864048) q[2];
rz(3.0388721) q[3];
sx q[3];
rz(-2.5892047) q[3];
sx q[3];
rz(1.8451167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];