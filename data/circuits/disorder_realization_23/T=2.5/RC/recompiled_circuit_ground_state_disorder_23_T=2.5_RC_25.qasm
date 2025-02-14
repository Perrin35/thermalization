OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5601114) q[0];
sx q[0];
rz(-1.2553296) q[0];
sx q[0];
rz(-2.6502967) q[0];
rz(-3.4499912) q[1];
sx q[1];
rz(1.8412794) q[1];
sx q[1];
rz(12.507764) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97293905) q[0];
sx q[0];
rz(-1.1298475) q[0];
sx q[0];
rz(-2.1827374) q[0];
rz(-pi) q[1];
rz(-2.0255857) q[2];
sx q[2];
rz(-0.84013772) q[2];
sx q[2];
rz(-2.0065713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.20929259) q[1];
sx q[1];
rz(-1.8718616) q[1];
sx q[1];
rz(1.73897) q[1];
x q[2];
rz(-3.0618151) q[3];
sx q[3];
rz(-1.8892131) q[3];
sx q[3];
rz(-1.7169881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0356174) q[2];
sx q[2];
rz(-2.3507037) q[2];
sx q[2];
rz(-0.38113511) q[2];
rz(-3.1086339) q[3];
sx q[3];
rz(-1.4573174) q[3];
sx q[3];
rz(-1.0491252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9829262) q[0];
sx q[0];
rz(-2.6148112) q[0];
sx q[0];
rz(0.43989936) q[0];
rz(3.0797709) q[1];
sx q[1];
rz(-1.6114176) q[1];
sx q[1];
rz(-0.27438146) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6950521) q[0];
sx q[0];
rz(-2.1529004) q[0];
sx q[0];
rz(-2.5869408) q[0];
x q[1];
rz(2.0437061) q[2];
sx q[2];
rz(-2.6691872) q[2];
sx q[2];
rz(-2.1771569) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3267131) q[1];
sx q[1];
rz(-2.68703) q[1];
sx q[1];
rz(1.3978895) q[1];
x q[2];
rz(-0.37526826) q[3];
sx q[3];
rz(-0.57884848) q[3];
sx q[3];
rz(-1.8414094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8476094) q[2];
sx q[2];
rz(-0.01394883) q[2];
sx q[2];
rz(3.0713165) q[2];
rz(0.62486068) q[3];
sx q[3];
rz(-1.3288386) q[3];
sx q[3];
rz(0.074987324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4570419) q[0];
sx q[0];
rz(-1.5032737) q[0];
sx q[0];
rz(2.7497838) q[0];
rz(-2.1458972) q[1];
sx q[1];
rz(-2.3052146) q[1];
sx q[1];
rz(0.044274274) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05737722) q[0];
sx q[0];
rz(-0.65859933) q[0];
sx q[0];
rz(2.088571) q[0];
x q[1];
rz(-2.2703553) q[2];
sx q[2];
rz(-0.92513621) q[2];
sx q[2];
rz(1.439656) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7024732) q[1];
sx q[1];
rz(-2.8608192) q[1];
sx q[1];
rz(-2.4149815) q[1];
x q[2];
rz(2.9999773) q[3];
sx q[3];
rz(-0.60185233) q[3];
sx q[3];
rz(-2.600649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5135045) q[2];
sx q[2];
rz(-2.2254483) q[2];
sx q[2];
rz(1.2543031) q[2];
rz(0.11145505) q[3];
sx q[3];
rz(-2.1696551) q[3];
sx q[3];
rz(0.85533065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8348389) q[0];
sx q[0];
rz(-2.4618537) q[0];
sx q[0];
rz(-0.87598959) q[0];
rz(-0.39729473) q[1];
sx q[1];
rz(-1.5715503) q[1];
sx q[1];
rz(1.3321336) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5529163) q[0];
sx q[0];
rz(-1.1015) q[0];
sx q[0];
rz(-1.7847212) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2418873) q[2];
sx q[2];
rz(-1.2061276) q[2];
sx q[2];
rz(-2.784035) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1173199) q[1];
sx q[1];
rz(-1.7242715) q[1];
sx q[1];
rz(1.491311) q[1];
x q[2];
rz(-1.187866) q[3];
sx q[3];
rz(-0.57102244) q[3];
sx q[3];
rz(1.3298705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1490271) q[2];
sx q[2];
rz(-1.488648) q[2];
sx q[2];
rz(-1.7472902) q[2];
rz(-2.1483138) q[3];
sx q[3];
rz(-1.1634049) q[3];
sx q[3];
rz(-1.8786028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.29436) q[0];
sx q[0];
rz(-0.17437409) q[0];
sx q[0];
rz(2.2827523) q[0];
rz(-0.42293388) q[1];
sx q[1];
rz(-2.6507288) q[1];
sx q[1];
rz(1.4609969) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2160667) q[0];
sx q[0];
rz(-0.51617814) q[0];
sx q[0];
rz(-0.25052089) q[0];
x q[1];
rz(1.5761574) q[2];
sx q[2];
rz(-2.6963384) q[2];
sx q[2];
rz(1.7536947) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36375052) q[1];
sx q[1];
rz(-1.7763174) q[1];
sx q[1];
rz(-1.3745893) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8247174) q[3];
sx q[3];
rz(-0.51957209) q[3];
sx q[3];
rz(-2.7412753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.6964263) q[2];
sx q[2];
rz(-0.40881613) q[2];
sx q[2];
rz(1.4985296) q[2];
rz(-2.0131352) q[3];
sx q[3];
rz(-1.107629) q[3];
sx q[3];
rz(-0.64277738) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4437101) q[0];
sx q[0];
rz(-1.3297798) q[0];
sx q[0];
rz(-0.71769303) q[0];
rz(-2.7806661) q[1];
sx q[1];
rz(-0.91054994) q[1];
sx q[1];
rz(-0.31401971) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234834) q[0];
sx q[0];
rz(-0.49719329) q[0];
sx q[0];
rz(3.0548993) q[0];
rz(-pi) q[1];
rz(1.5271565) q[2];
sx q[2];
rz(-2.4167762) q[2];
sx q[2];
rz(-0.61312719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2870658) q[1];
sx q[1];
rz(-2.5662133) q[1];
sx q[1];
rz(2.1772312) q[1];
x q[2];
rz(1.3889503) q[3];
sx q[3];
rz(-0.80917984) q[3];
sx q[3];
rz(-1.941118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6641984) q[2];
sx q[2];
rz(-1.7919284) q[2];
sx q[2];
rz(-2.9273709) q[2];
rz(-0.91841206) q[3];
sx q[3];
rz(-2.1954506) q[3];
sx q[3];
rz(-1.4001747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4825831) q[0];
sx q[0];
rz(-2.5795586) q[0];
sx q[0];
rz(2.5157978) q[0];
rz(1.2311426) q[1];
sx q[1];
rz(-1.8742671) q[1];
sx q[1];
rz(2.321718) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5441262) q[0];
sx q[0];
rz(-1.0080999) q[0];
sx q[0];
rz(1.94348) q[0];
x q[1];
rz(2.0076114) q[2];
sx q[2];
rz(-2.2333418) q[2];
sx q[2];
rz(0.32811859) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.956683) q[1];
sx q[1];
rz(-1.239983) q[1];
sx q[1];
rz(-0.079903101) q[1];
x q[2];
rz(0.19015892) q[3];
sx q[3];
rz(-0.09626046) q[3];
sx q[3];
rz(3.0031693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.95255533) q[2];
sx q[2];
rz(-1.8841212) q[2];
sx q[2];
rz(2.7541568) q[2];
rz(2.6881325) q[3];
sx q[3];
rz(-2.267434) q[3];
sx q[3];
rz(-0.15997729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.74359918) q[0];
sx q[0];
rz(-1.1389808) q[0];
sx q[0];
rz(-1.1238264) q[0];
rz(-0.75346142) q[1];
sx q[1];
rz(-1.9517784) q[1];
sx q[1];
rz(-1.168728) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2126287) q[0];
sx q[0];
rz(-0.88200404) q[0];
sx q[0];
rz(-3.0228011) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5587135) q[2];
sx q[2];
rz(-0.45736499) q[2];
sx q[2];
rz(-1.7810019) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.469731) q[1];
sx q[1];
rz(-2.1657092) q[1];
sx q[1];
rz(-0.39284006) q[1];
rz(-pi) q[2];
rz(-0.99645741) q[3];
sx q[3];
rz(-0.74241168) q[3];
sx q[3];
rz(-1.3848828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0231102) q[2];
sx q[2];
rz(-1.0175984) q[2];
sx q[2];
rz(0.47038308) q[2];
rz(-1.4605626) q[3];
sx q[3];
rz(-1.8813671) q[3];
sx q[3];
rz(-2.1605261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7505782) q[0];
sx q[0];
rz(-1.3704726) q[0];
sx q[0];
rz(1.6455261) q[0];
rz(-1.1356614) q[1];
sx q[1];
rz(-1.2874425) q[1];
sx q[1];
rz(-2.6527203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5565692) q[0];
sx q[0];
rz(-2.9980066) q[0];
sx q[0];
rz(-1.7973336) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0592986) q[2];
sx q[2];
rz(-1.6548205) q[2];
sx q[2];
rz(2.3531599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1355537) q[1];
sx q[1];
rz(-0.054243739) q[1];
sx q[1];
rz(-0.42363055) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5034204) q[3];
sx q[3];
rz(-0.95219288) q[3];
sx q[3];
rz(-2.2139795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6966072) q[2];
sx q[2];
rz(-2.5312436) q[2];
sx q[2];
rz(-1.102591) q[2];
rz(1.6164814) q[3];
sx q[3];
rz(-2.3307255) q[3];
sx q[3];
rz(1.6703687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1333604) q[0];
sx q[0];
rz(-0.52953774) q[0];
sx q[0];
rz(1.8089005) q[0];
rz(-0.38707271) q[1];
sx q[1];
rz(-1.8862855) q[1];
sx q[1];
rz(1.056384) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1176193) q[0];
sx q[0];
rz(-1.4314382) q[0];
sx q[0];
rz(-1.3246714) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6956639) q[2];
sx q[2];
rz(-1.1556632) q[2];
sx q[2];
rz(1.268569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3247128) q[1];
sx q[1];
rz(-1.9512647) q[1];
sx q[1];
rz(-0.76075424) q[1];
x q[2];
rz(-0.71223082) q[3];
sx q[3];
rz(-2.3281186) q[3];
sx q[3];
rz(-1.939919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.010783823) q[2];
sx q[2];
rz(-2.2937412) q[2];
sx q[2];
rz(1.7424142) q[2];
rz(2.0684659) q[3];
sx q[3];
rz(-1.9173887) q[3];
sx q[3];
rz(-2.6789902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90049967) q[0];
sx q[0];
rz(-1.4971965) q[0];
sx q[0];
rz(1.5370488) q[0];
rz(0.72147876) q[1];
sx q[1];
rz(-0.57054467) q[1];
sx q[1];
rz(1.9346938) q[1];
rz(0.73312326) q[2];
sx q[2];
rz(-2.5903268) q[2];
sx q[2];
rz(-3.0916284) q[2];
rz(1.4720451) q[3];
sx q[3];
rz(-0.80174123) q[3];
sx q[3];
rz(1.5855736) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
