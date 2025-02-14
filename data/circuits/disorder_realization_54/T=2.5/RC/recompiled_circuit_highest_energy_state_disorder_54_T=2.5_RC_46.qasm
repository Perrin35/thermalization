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
rz(2.8132791) q[0];
sx q[0];
rz(-2.0067196) q[0];
sx q[0];
rz(-0.56587926) q[0];
rz(-2.8973051) q[1];
sx q[1];
rz(-1.7841508) q[1];
sx q[1];
rz(2.607333) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8847647) q[0];
sx q[0];
rz(-2.3951911) q[0];
sx q[0];
rz(1.166626) q[0];
x q[1];
rz(-1.8741605) q[2];
sx q[2];
rz(-2.144226) q[2];
sx q[2];
rz(-2.0801185) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54142648) q[1];
sx q[1];
rz(-1.8712338) q[1];
sx q[1];
rz(-1.5460127) q[1];
rz(-pi) q[2];
rz(-2.665042) q[3];
sx q[3];
rz(-0.50560564) q[3];
sx q[3];
rz(-1.4534392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1285706) q[2];
sx q[2];
rz(-0.2505005) q[2];
sx q[2];
rz(-2.8924083) q[2];
rz(-1.3439882) q[3];
sx q[3];
rz(-1.598571) q[3];
sx q[3];
rz(0.070076076) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4942112) q[0];
sx q[0];
rz(-1.2263612) q[0];
sx q[0];
rz(2.1517854) q[0];
rz(-0.79065943) q[1];
sx q[1];
rz(-2.374687) q[1];
sx q[1];
rz(2.8299423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30325952) q[0];
sx q[0];
rz(-1.0619535) q[0];
sx q[0];
rz(0.80755393) q[0];
x q[1];
rz(-2.0925518) q[2];
sx q[2];
rz(-2.9269896) q[2];
sx q[2];
rz(1.3133446) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5144044) q[1];
sx q[1];
rz(-1.7429105) q[1];
sx q[1];
rz(-1.3283511) q[1];
x q[2];
rz(-0.34087025) q[3];
sx q[3];
rz(-2.1580527) q[3];
sx q[3];
rz(2.1275008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1404169) q[2];
sx q[2];
rz(-1.1602465) q[2];
sx q[2];
rz(-0.95138335) q[2];
rz(0.22826711) q[3];
sx q[3];
rz(-2.164866) q[3];
sx q[3];
rz(1.2736646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839639) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(-0.11446318) q[0];
rz(1.0022256) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(-0.1618298) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069050633) q[0];
sx q[0];
rz(-2.4009573) q[0];
sx q[0];
rz(-0.97461755) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21185565) q[2];
sx q[2];
rz(-0.39688928) q[2];
sx q[2];
rz(0.8241764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56282957) q[1];
sx q[1];
rz(-2.17318) q[1];
sx q[1];
rz(1.2788601) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3409445) q[3];
sx q[3];
rz(-1.8654239) q[3];
sx q[3];
rz(1.1931452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3942922) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(-0.14075819) q[2];
rz(-0.55177871) q[3];
sx q[3];
rz(-1.2740302) q[3];
sx q[3];
rz(-0.75132918) q[3];
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
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20525876) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(0.67489135) q[0];
rz(2.9871125) q[1];
sx q[1];
rz(-1.4090425) q[1];
sx q[1];
rz(-2.6836269) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3727399) q[0];
sx q[0];
rz(-0.64592664) q[0];
sx q[0];
rz(-2.4619308) q[0];
x q[1];
rz(2.0003597) q[2];
sx q[2];
rz(-1.8515966) q[2];
sx q[2];
rz(0.28102885) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.030589015) q[1];
sx q[1];
rz(-1.2007502) q[1];
sx q[1];
rz(-1.7683328) q[1];
x q[2];
rz(1.4639888) q[3];
sx q[3];
rz(-1.065101) q[3];
sx q[3];
rz(0.13733521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0080879) q[2];
sx q[2];
rz(-2.4613481) q[2];
sx q[2];
rz(-2.8221455) q[2];
rz(1.3301298) q[3];
sx q[3];
rz(-1.0165241) q[3];
sx q[3];
rz(0.5602079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9820246) q[0];
sx q[0];
rz(-1.1363131) q[0];
sx q[0];
rz(-1.2215479) q[0];
rz(-2.9087861) q[1];
sx q[1];
rz(-0.56929749) q[1];
sx q[1];
rz(0.98308841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.074551) q[0];
sx q[0];
rz(-2.3877151) q[0];
sx q[0];
rz(0.75410944) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4120502) q[2];
sx q[2];
rz(-2.5149143) q[2];
sx q[2];
rz(1.9959297) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32199643) q[1];
sx q[1];
rz(-0.23769544) q[1];
sx q[1];
rz(-0.6657712) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9931204) q[3];
sx q[3];
rz(-1.7230801) q[3];
sx q[3];
rz(1.1642758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7588707) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(2.7663084) q[2];
rz(-2.0659857) q[3];
sx q[3];
rz(-0.38638249) q[3];
sx q[3];
rz(0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7406834) q[0];
sx q[0];
rz(-1.704819) q[0];
sx q[0];
rz(-1.7042473) q[0];
rz(0.078723343) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(0.54814235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.948471) q[0];
sx q[0];
rz(-1.3594432) q[0];
sx q[0];
rz(2.4708492) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0847381) q[2];
sx q[2];
rz(-2.0693739) q[2];
sx q[2];
rz(-0.14190292) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3972612) q[1];
sx q[1];
rz(-1.9588699) q[1];
sx q[1];
rz(-0.71340386) q[1];
rz(-pi) q[2];
rz(2.4924566) q[3];
sx q[3];
rz(-2.2723291) q[3];
sx q[3];
rz(2.8809406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4569725) q[2];
sx q[2];
rz(-2.520884) q[2];
sx q[2];
rz(-0.6753298) q[2];
rz(1.5572549) q[3];
sx q[3];
rz(-1.8341589) q[3];
sx q[3];
rz(0.61711446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5439344) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(-2.6455998) q[0];
rz(1.687382) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(-2.0248263) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3163213) q[0];
sx q[0];
rz(-2.6351235) q[0];
sx q[0];
rz(2.3536918) q[0];
rz(0.26979763) q[2];
sx q[2];
rz(-1.0493663) q[2];
sx q[2];
rz(-2.9228589) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8161113) q[1];
sx q[1];
rz(-1.2524091) q[1];
sx q[1];
rz(0.35811425) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6326866) q[3];
sx q[3];
rz(-1.1155432) q[3];
sx q[3];
rz(0.52857182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3108959) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(2.7867219) q[2];
rz(-0.76534671) q[3];
sx q[3];
rz(-0.86330515) q[3];
sx q[3];
rz(-3.0520458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.68882051) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(-2.3759957) q[0];
rz(2.2701524) q[1];
sx q[1];
rz(-2.3759418) q[1];
sx q[1];
rz(1.8188459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2521556) q[0];
sx q[0];
rz(-1.3103292) q[0];
sx q[0];
rz(-1.5171264) q[0];
x q[1];
rz(-0.89143388) q[2];
sx q[2];
rz(-1.3855033) q[2];
sx q[2];
rz(-0.44033716) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8606931) q[1];
sx q[1];
rz(-2.6329663) q[1];
sx q[1];
rz(-0.033837576) q[1];
x q[2];
rz(1.7962097) q[3];
sx q[3];
rz(-1.3691878) q[3];
sx q[3];
rz(-2.9359093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65779296) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(0.22171177) q[2];
rz(-2.0486369) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(-2.4634585) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073864989) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(-2.8644323) q[0];
rz(2.3932338) q[1];
sx q[1];
rz(-0.44735083) q[1];
sx q[1];
rz(-1.998273) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7171913) q[0];
sx q[0];
rz(-2.7729308) q[0];
sx q[0];
rz(-2.1564846) q[0];
x q[1];
rz(1.3242993) q[2];
sx q[2];
rz(-0.72028226) q[2];
sx q[2];
rz(1.5697073) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1705679) q[1];
sx q[1];
rz(-1.2538365) q[1];
sx q[1];
rz(-0.16522878) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42871126) q[3];
sx q[3];
rz(-2.747251) q[3];
sx q[3];
rz(-1.1129781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9144168) q[2];
sx q[2];
rz(-2.1180426) q[2];
sx q[2];
rz(-3.0926256) q[2];
rz(2.07552) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(2.9858885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12298909) q[0];
sx q[0];
rz(-1.5276696) q[0];
sx q[0];
rz(0.019512026) q[0];
rz(2.3628269) q[1];
sx q[1];
rz(-2.4486783) q[1];
sx q[1];
rz(-1.5787517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21694788) q[0];
sx q[0];
rz(-2.9836982) q[0];
sx q[0];
rz(-0.96952207) q[0];
x q[1];
rz(-2.1807272) q[2];
sx q[2];
rz(-1.4205407) q[2];
sx q[2];
rz(-0.021266887) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.73247468) q[1];
sx q[1];
rz(-2.1555087) q[1];
sx q[1];
rz(-0.52701108) q[1];
rz(-pi) q[2];
rz(-1.7354129) q[3];
sx q[3];
rz(-2.4009628) q[3];
sx q[3];
rz(0.1333789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19266263) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(-0.16555244) q[2];
rz(2.5340269) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(-2.6732388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.7623357) q[0];
sx q[0];
rz(-0.48185928) q[0];
sx q[0];
rz(3.0961105) q[0];
rz(-1.5070076) q[1];
sx q[1];
rz(-1.6804809) q[1];
sx q[1];
rz(1.5483821) q[1];
rz(2.9956837) q[2];
sx q[2];
rz(-0.78073954) q[2];
sx q[2];
rz(-2.2013046) q[2];
rz(-0.28260091) q[3];
sx q[3];
rz(-0.62487124) q[3];
sx q[3];
rz(1.665653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
