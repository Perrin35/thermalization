OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(-0.50146377) q[0];
rz(1.4986829) q[1];
sx q[1];
rz(-2.745435) q[1];
sx q[1];
rz(-0.3224386) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3381391) q[0];
sx q[0];
rz(-2.8087466) q[0];
sx q[0];
rz(-1.8344318) q[0];
rz(-2.5311004) q[2];
sx q[2];
rz(-2.1571026) q[2];
sx q[2];
rz(2.5551978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7079561) q[1];
sx q[1];
rz(-0.70322137) q[1];
sx q[1];
rz(-1.0183079) q[1];
rz(-0.93482165) q[3];
sx q[3];
rz(-1.1879731) q[3];
sx q[3];
rz(2.8179907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(-2.5855529) q[2];
rz(2.3089144) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(0.26113025) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(3.0325586) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.581668) q[0];
sx q[0];
rz(-1.7413057) q[0];
sx q[0];
rz(-1.8018363) q[0];
x q[1];
rz(-2.878506) q[2];
sx q[2];
rz(-1.7619942) q[2];
sx q[2];
rz(-2.4441602) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8728767) q[1];
sx q[1];
rz(-2.0375588) q[1];
sx q[1];
rz(0.66653911) q[1];
rz(-0.29626366) q[3];
sx q[3];
rz(-0.54013541) q[3];
sx q[3];
rz(1.4512645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2382425) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(-1.2634574) q[2];
rz(2.8144828) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(-1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0771714) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.3431312) q[0];
rz(-2.893977) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(-2.6599191) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96268481) q[0];
sx q[0];
rz(-0.6324397) q[0];
sx q[0];
rz(-2.2326367) q[0];
rz(-pi) q[1];
rz(-1.1142251) q[2];
sx q[2];
rz(-2.4764428) q[2];
sx q[2];
rz(0.84012023) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3232705) q[1];
sx q[1];
rz(-0.96410492) q[1];
sx q[1];
rz(-2.2491749) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2504911) q[3];
sx q[3];
rz(-0.84935969) q[3];
sx q[3];
rz(0.47618714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3383011) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(2.6417007) q[2];
rz(0.56097427) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.4612173) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9445779) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(0.552185) q[0];
rz(1.5532956) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(-1.8968556) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9024076) q[0];
sx q[0];
rz(-0.33948487) q[0];
sx q[0];
rz(-3.1367338) q[0];
rz(-pi) q[1];
rz(-0.130529) q[2];
sx q[2];
rz(-1.7703238) q[2];
sx q[2];
rz(-0.56632698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4501805) q[1];
sx q[1];
rz(-1.0346864) q[1];
sx q[1];
rz(-0.24713534) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5707385) q[3];
sx q[3];
rz(-1.3802765) q[3];
sx q[3];
rz(-1.966147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2923979) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(-1.4771279) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(-0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078159049) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(3.0601236) q[0];
rz(-3.07913) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(1.5030456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7446639) q[0];
sx q[0];
rz(-1.443112) q[0];
sx q[0];
rz(-0.98608195) q[0];
x q[1];
rz(3.0585074) q[2];
sx q[2];
rz(-1.462888) q[2];
sx q[2];
rz(-1.8748869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.11322184) q[1];
sx q[1];
rz(-1.8824667) q[1];
sx q[1];
rz(-3.034364) q[1];
x q[2];
rz(0.50932192) q[3];
sx q[3];
rz(-1.29042) q[3];
sx q[3];
rz(2.544739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6333255) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(-1.2379237) q[2];
rz(-1.1226908) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6102819) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(2.8748728) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(-2.3430603) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5312936) q[0];
sx q[0];
rz(-1.4302505) q[0];
sx q[0];
rz(1.2929686) q[0];
rz(-pi) q[1];
rz(1.919826) q[2];
sx q[2];
rz(-1.0997699) q[2];
sx q[2];
rz(0.40086056) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.861607) q[1];
sx q[1];
rz(-1.0000739) q[1];
sx q[1];
rz(-0.062203783) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1875601) q[3];
sx q[3];
rz(-0.87152374) q[3];
sx q[3];
rz(-1.1773393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(0.36671656) q[2];
rz(-1.2612873) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40903184) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(-0.60638705) q[0];
rz(-0.19730332) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(0.46404776) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460829) q[0];
sx q[0];
rz(-1.1469139) q[0];
sx q[0];
rz(0.8830107) q[0];
rz(-1.3299994) q[2];
sx q[2];
rz(-1.5511998) q[2];
sx q[2];
rz(-0.400825) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5210919) q[1];
sx q[1];
rz(-1.2695841) q[1];
sx q[1];
rz(-2.0786933) q[1];
rz(-pi) q[2];
rz(-1.4887605) q[3];
sx q[3];
rz(-1.335252) q[3];
sx q[3];
rz(1.1788648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.93418926) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(-0.25804538) q[2];
rz(1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(0.38129693) q[0];
rz(3.0463468) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.7000748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1822752) q[0];
sx q[0];
rz(-1.5858985) q[0];
sx q[0];
rz(3.0781151) q[0];
rz(-pi) q[1];
rz(1.6091789) q[2];
sx q[2];
rz(-2.1830609) q[2];
sx q[2];
rz(1.6910764) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4691094) q[1];
sx q[1];
rz(-1.4454578) q[1];
sx q[1];
rz(3.010716) q[1];
rz(-pi) q[2];
rz(-3.0184047) q[3];
sx q[3];
rz(-0.33629575) q[3];
sx q[3];
rz(1.5598701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(0.44858027) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(0.8297689) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.7425849) q[0];
rz(-0.31708583) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(0.98659602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50842677) q[0];
sx q[0];
rz(-1.0796483) q[0];
sx q[0];
rz(-1.7500061) q[0];
x q[1];
rz(0.95790205) q[2];
sx q[2];
rz(-2.3447678) q[2];
sx q[2];
rz(0.56607841) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36531891) q[1];
sx q[1];
rz(-2.1734108) q[1];
sx q[1];
rz(-1.7425294) q[1];
x q[2];
rz(-2.968623) q[3];
sx q[3];
rz(-1.9860455) q[3];
sx q[3];
rz(-3.1239307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9265147) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(-0.26327291) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(1.4051399) q[0];
rz(-2.3204904) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(-1.4155037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6380499) q[0];
sx q[0];
rz(-1.3347795) q[0];
sx q[0];
rz(-0.34061265) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1591522) q[2];
sx q[2];
rz(-2.7124565) q[2];
sx q[2];
rz(0.8796126) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0701329) q[1];
sx q[1];
rz(-1.3576641) q[1];
sx q[1];
rz(0.15503426) q[1];
rz(1.2763001) q[3];
sx q[3];
rz(-0.53139549) q[3];
sx q[3];
rz(-2.6607799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71904174) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(3.0174875) q[2];
rz(2.1758046) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(-2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.538095) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(-0.30766906) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(-0.60097354) q[2];
sx q[2];
rz(-2.8291611) q[2];
sx q[2];
rz(-0.57401382) q[2];
rz(-0.55660558) q[3];
sx q[3];
rz(-1.442933) q[3];
sx q[3];
rz(1.8419151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
