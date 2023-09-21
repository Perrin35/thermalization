OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(0.00014076509) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(1.948184) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084273987) q[0];
sx q[0];
rz(-0.3398551) q[0];
sx q[0];
rz(1.0124595) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2725856) q[2];
sx q[2];
rz(-1.0422921) q[2];
sx q[2];
rz(2.3117711) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6686033) q[1];
sx q[1];
rz(-1.1057292) q[1];
sx q[1];
rz(-0.71778645) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0334929) q[3];
sx q[3];
rz(-0.24216147) q[3];
sx q[3];
rz(2.7778181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6821735) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(1.9127282) q[2];
rz(-1.4131644) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(-1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380149) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(1.0128101) q[0];
rz(-0.027659841) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-2.0181296) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7682122) q[0];
sx q[0];
rz(-1.1102311) q[0];
sx q[0];
rz(0.064231355) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6337109) q[2];
sx q[2];
rz(-0.79195576) q[2];
sx q[2];
rz(-0.5069678) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23501913) q[1];
sx q[1];
rz(-2.0779013) q[1];
sx q[1];
rz(-2.1706235) q[1];
x q[2];
rz(-2.1334322) q[3];
sx q[3];
rz(-0.99346549) q[3];
sx q[3];
rz(1.2853704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79364395) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(-0.91903764) q[2];
rz(-0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(-1.8664237) q[0];
rz(0.69349849) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(-2.0085874) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64273524) q[0];
sx q[0];
rz(-2.9992636) q[0];
sx q[0];
rz(1.5940773) q[0];
rz(2.0793545) q[2];
sx q[2];
rz(-2.2937751) q[2];
sx q[2];
rz(1.522097) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6269835) q[1];
sx q[1];
rz(-1.2830462) q[1];
sx q[1];
rz(2.1327553) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95136178) q[3];
sx q[3];
rz(-1.0361443) q[3];
sx q[3];
rz(-1.149328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2514078) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(-1.8481002) q[2];
rz(0.039316468) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.2599729) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(-1.1608634) q[0];
rz(-0.89598957) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(-0.13555759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40977749) q[0];
sx q[0];
rz(-1.8827794) q[0];
sx q[0];
rz(0.72809763) q[0];
x q[1];
rz(0.074684871) q[2];
sx q[2];
rz(-1.1876145) q[2];
sx q[2];
rz(0.89786868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5038824) q[1];
sx q[1];
rz(-1.4117068) q[1];
sx q[1];
rz(1.7590982) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1578219) q[3];
sx q[3];
rz(-1.2994248) q[3];
sx q[3];
rz(-1.8872758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23665145) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(-3.0974292) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(-2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(-2.1039466) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(2.1283545) q[0];
rz(3.0918616) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(2.0577046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0439107) q[0];
sx q[0];
rz(-1.8143166) q[0];
sx q[0];
rz(-2.9527412) q[0];
x q[1];
rz(0.82917825) q[2];
sx q[2];
rz(-2.7784756) q[2];
sx q[2];
rz(-1.5311637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10478445) q[1];
sx q[1];
rz(-1.0462865) q[1];
sx q[1];
rz(0.12496897) q[1];
rz(-pi) q[2];
rz(-1.5203939) q[3];
sx q[3];
rz(-2.084123) q[3];
sx q[3];
rz(0.36171519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9087387) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(-0.24442913) q[2];
rz(-2.7092253) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(-2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(-3.0474512) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(-0.89541268) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0462414) q[0];
sx q[0];
rz(-1.5329251) q[0];
sx q[0];
rz(-0.34237679) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63352872) q[2];
sx q[2];
rz(-1.765536) q[2];
sx q[2];
rz(2.5031236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.10455924) q[1];
sx q[1];
rz(-1.7393141) q[1];
sx q[1];
rz(3.0969572) q[1];
rz(-pi) q[2];
rz(1.496109) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(0.96935779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0078997) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(-2.3383979) q[2];
rz(1.1903654) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0728834) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(-0.51914006) q[0];
rz(2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(1.8849467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0518258) q[0];
sx q[0];
rz(-1.3194114) q[0];
sx q[0];
rz(-3.0999523) q[0];
x q[1];
rz(-1.16876) q[2];
sx q[2];
rz(-0.49905825) q[2];
sx q[2];
rz(1.7390342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9650967) q[1];
sx q[1];
rz(-0.28897044) q[1];
sx q[1];
rz(-0.22775905) q[1];
x q[2];
rz(0.36862936) q[3];
sx q[3];
rz(-0.7437403) q[3];
sx q[3];
rz(1.3524692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1533623) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68650866) q[0];
sx q[0];
rz(-0.63502705) q[0];
sx q[0];
rz(-1.1839266) q[0];
rz(-pi) q[1];
rz(2.2690291) q[2];
sx q[2];
rz(-1.1485032) q[2];
sx q[2];
rz(2.1625105) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-2.4657001) q[1];
sx q[1];
rz(3.0216316) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49658637) q[3];
sx q[3];
rz(-1.7644617) q[3];
sx q[3];
rz(-0.95369875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(-0.44000885) q[2];
rz(-2.4258339) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(-2.0619152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0004262) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(2.0429042) q[0];
rz(-0.72775841) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(3.0922906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44137529) q[0];
sx q[0];
rz(-0.8964296) q[0];
sx q[0];
rz(2.3250439) q[0];
rz(0.98408913) q[2];
sx q[2];
rz(-2.3684635) q[2];
sx q[2];
rz(-3.0423321) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1972563) q[1];
sx q[1];
rz(-1.1068871) q[1];
sx q[1];
rz(1.0524366) q[1];
rz(2.0588166) q[3];
sx q[3];
rz(-0.89142311) q[3];
sx q[3];
rz(0.29619141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0631642) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(-0.72193974) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(-2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9938875) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(-2.06185) q[0];
rz(-1.059277) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(1.7396897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.429075) q[0];
sx q[0];
rz(-1.6970465) q[0];
sx q[0];
rz(-1.7849126) q[0];
rz(0.14342587) q[2];
sx q[2];
rz(-0.18866814) q[2];
sx q[2];
rz(-0.27486899) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5857196) q[1];
sx q[1];
rz(-1.266529) q[1];
sx q[1];
rz(0.85822206) q[1];
rz(-pi) q[2];
rz(1.7781156) q[3];
sx q[3];
rz(-0.94982409) q[3];
sx q[3];
rz(0.45351115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(1.6213017) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.993492) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(-2.2254754) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(-0.73268391) q[2];
sx q[2];
rz(-0.97326836) q[2];
sx q[2];
rz(1.2748981) q[2];
rz(2.0879073) q[3];
sx q[3];
rz(-1.5862982) q[3];
sx q[3];
rz(1.3261212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
