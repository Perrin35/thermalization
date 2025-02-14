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
rz(1.16601) q[0];
sx q[0];
rz(2.2628885) q[0];
sx q[0];
rz(9.2287697) q[0];
rz(1.6485543) q[1];
sx q[1];
rz(-2.2033446) q[1];
sx q[1];
rz(-0.87747639) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24841079) q[0];
sx q[0];
rz(-2.0201004) q[0];
sx q[0];
rz(0.23907678) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9196366) q[2];
sx q[2];
rz(-1.9720805) q[2];
sx q[2];
rz(0.87043412) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2171195) q[1];
sx q[1];
rz(-1.3018601) q[1];
sx q[1];
rz(2.0321789) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0935173) q[3];
sx q[3];
rz(-1.358489) q[3];
sx q[3];
rz(-2.4119854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8087372) q[2];
sx q[2];
rz(-2.4337807) q[2];
sx q[2];
rz(3.1044002) q[2];
rz(2.228179) q[3];
sx q[3];
rz(-0.98228684) q[3];
sx q[3];
rz(-1.5090401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7971802) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(-0.2221701) q[0];
rz(-0.19368681) q[1];
sx q[1];
rz(-1.8733459) q[1];
sx q[1];
rz(2.8817315) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7714789) q[0];
sx q[0];
rz(-0.51379097) q[0];
sx q[0];
rz(-2.0822078) q[0];
rz(-pi) q[1];
rz(-0.1078306) q[2];
sx q[2];
rz(-1.7485776) q[2];
sx q[2];
rz(-1.9974634) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9495595) q[1];
sx q[1];
rz(-1.7550526) q[1];
sx q[1];
rz(-2.9263587) q[1];
rz(-pi) q[2];
rz(-2.0570898) q[3];
sx q[3];
rz(-2.2688365) q[3];
sx q[3];
rz(-1.2459038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3746609) q[2];
sx q[2];
rz(-2.5587475) q[2];
sx q[2];
rz(0.15677162) q[2];
rz(2.3356656) q[3];
sx q[3];
rz(-1.0904049) q[3];
sx q[3];
rz(2.6167615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90461516) q[0];
sx q[0];
rz(-0.1122864) q[0];
sx q[0];
rz(1.1861381) q[0];
rz(0.063749464) q[1];
sx q[1];
rz(-1.6155764) q[1];
sx q[1];
rz(-2.729111) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97226364) q[0];
sx q[0];
rz(-1.0522075) q[0];
sx q[0];
rz(-0.56904582) q[0];
rz(-0.15699129) q[2];
sx q[2];
rz(-2.2447033) q[2];
sx q[2];
rz(1.5012596) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3220249) q[1];
sx q[1];
rz(-1.7734999) q[1];
sx q[1];
rz(2.8611819) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2687421) q[3];
sx q[3];
rz(-2.344909) q[3];
sx q[3];
rz(2.7484675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.56846642) q[2];
sx q[2];
rz(-0.44034475) q[2];
sx q[2];
rz(1.9754515) q[2];
rz(0.79683534) q[3];
sx q[3];
rz(-1.7723869) q[3];
sx q[3];
rz(0.78127512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9154938) q[0];
sx q[0];
rz(-1.6330566) q[0];
sx q[0];
rz(1.4133806) q[0];
rz(2.6426897) q[1];
sx q[1];
rz(-1.3830802) q[1];
sx q[1];
rz(1.2938719) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.897936) q[0];
sx q[0];
rz(-1.5439646) q[0];
sx q[0];
rz(-1.6559589) q[0];
x q[1];
rz(-1.8982688) q[2];
sx q[2];
rz(-1.4988384) q[2];
sx q[2];
rz(-1.5639992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41999757) q[1];
sx q[1];
rz(-1.5087869) q[1];
sx q[1];
rz(-1.5676143) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7440552) q[3];
sx q[3];
rz(-3.0500055) q[3];
sx q[3];
rz(-2.0360144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52900195) q[2];
sx q[2];
rz(-2.4606073) q[2];
sx q[2];
rz(-2.9717818) q[2];
rz(-2.7134907) q[3];
sx q[3];
rz(-1.4295108) q[3];
sx q[3];
rz(2.9812109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47996461) q[0];
sx q[0];
rz(-2.7492838) q[0];
sx q[0];
rz(2.3062134) q[0];
rz(-2.5021878) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(1.0909874) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3314455) q[0];
sx q[0];
rz(-1.5584599) q[0];
sx q[0];
rz(-2.4452565) q[0];
x q[1];
rz(2.8418737) q[2];
sx q[2];
rz(-1.4106361) q[2];
sx q[2];
rz(2.0998139) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9613232) q[1];
sx q[1];
rz(-1.1774447) q[1];
sx q[1];
rz(2.250227) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32650487) q[3];
sx q[3];
rz(-2.0707664) q[3];
sx q[3];
rz(-0.96159354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7731446) q[2];
sx q[2];
rz(-1.1870563) q[2];
sx q[2];
rz(-3.1375569) q[2];
rz(-2.0648547) q[3];
sx q[3];
rz(-2.4439947) q[3];
sx q[3];
rz(2.9546886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3183597) q[0];
sx q[0];
rz(-1.7802745) q[0];
sx q[0];
rz(2.5079492) q[0];
rz(2.7404495) q[1];
sx q[1];
rz(-0.70924962) q[1];
sx q[1];
rz(0.69222442) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61718231) q[0];
sx q[0];
rz(-1.9815191) q[0];
sx q[0];
rz(2.3401712) q[0];
rz(0.76898076) q[2];
sx q[2];
rz(-1.7959463) q[2];
sx q[2];
rz(1.0234635) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8521523) q[1];
sx q[1];
rz(-1.0045144) q[1];
sx q[1];
rz(-0.29636611) q[1];
rz(-pi) q[2];
rz(-2.4603087) q[3];
sx q[3];
rz(-2.9046106) q[3];
sx q[3];
rz(-2.9221688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.20574084) q[2];
sx q[2];
rz(-2.3571099) q[2];
sx q[2];
rz(1.240823) q[2];
rz(1.1848909) q[3];
sx q[3];
rz(-1.840206) q[3];
sx q[3];
rz(-1.544781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.009509) q[0];
sx q[0];
rz(-1.8385831) q[0];
sx q[0];
rz(-1.7708154) q[0];
rz(2.7235183) q[1];
sx q[1];
rz(-1.4825753) q[1];
sx q[1];
rz(0.48666993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3330227) q[0];
sx q[0];
rz(-1.6679224) q[0];
sx q[0];
rz(1.0104695) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7667747) q[2];
sx q[2];
rz(-2.0897802) q[2];
sx q[2];
rz(-1.5216375) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6981737) q[1];
sx q[1];
rz(-1.3935745) q[1];
sx q[1];
rz(1.7755356) q[1];
rz(-pi) q[2];
rz(1.5237996) q[3];
sx q[3];
rz(-2.040427) q[3];
sx q[3];
rz(2.3835973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.439094) q[2];
sx q[2];
rz(-2.0926026) q[2];
sx q[2];
rz(0.19719633) q[2];
rz(2.6321453) q[3];
sx q[3];
rz(-2.8809437) q[3];
sx q[3];
rz(-2.313224) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8615123) q[0];
sx q[0];
rz(-1.0741638) q[0];
sx q[0];
rz(1.7751088) q[0];
rz(-2.3249783) q[1];
sx q[1];
rz(-1.0895224) q[1];
sx q[1];
rz(2.8588967) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0007083) q[0];
sx q[0];
rz(-1.3907897) q[0];
sx q[0];
rz(-0.96830826) q[0];
rz(-pi) q[1];
rz(0.28699283) q[2];
sx q[2];
rz(-0.32397917) q[2];
sx q[2];
rz(-2.2215077) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1360838) q[1];
sx q[1];
rz(-0.89673954) q[1];
sx q[1];
rz(0.23917178) q[1];
rz(-1.181284) q[3];
sx q[3];
rz(-2.0593025) q[3];
sx q[3];
rz(2.8581885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.035159811) q[2];
sx q[2];
rz(-2.1356434) q[2];
sx q[2];
rz(-0.98313588) q[2];
rz(1.2784917) q[3];
sx q[3];
rz(-2.216414) q[3];
sx q[3];
rz(-2.2127051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4814257) q[0];
sx q[0];
rz(-1.3496512) q[0];
sx q[0];
rz(-2.8283258) q[0];
rz(-2.616864) q[1];
sx q[1];
rz(-0.26291651) q[1];
sx q[1];
rz(-1.252334) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1838054) q[0];
sx q[0];
rz(-2.5982214) q[0];
sx q[0];
rz(-2.9068391) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5400794) q[2];
sx q[2];
rz(-0.85962112) q[2];
sx q[2];
rz(1.892923) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5460427) q[1];
sx q[1];
rz(-0.79985207) q[1];
sx q[1];
rz(-2.355666) q[1];
x q[2];
rz(1.657293) q[3];
sx q[3];
rz(-0.75631006) q[3];
sx q[3];
rz(-2.5027517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9423882) q[2];
sx q[2];
rz(-0.62493268) q[2];
sx q[2];
rz(-1.3814629) q[2];
rz(2.840461) q[3];
sx q[3];
rz(-0.98265177) q[3];
sx q[3];
rz(-2.7605831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(-2.9075266) q[0];
sx q[0];
rz(-2.1387687) q[0];
sx q[0];
rz(1.7940849) q[0];
rz(-0.8849591) q[1];
sx q[1];
rz(-0.62565175) q[1];
sx q[1];
rz(-1.398483) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2725582) q[0];
sx q[0];
rz(-1.5601275) q[0];
sx q[0];
rz(1.9452626) q[0];
x q[1];
rz(-0.18357205) q[2];
sx q[2];
rz(-2.2677448) q[2];
sx q[2];
rz(-1.1810034) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81512886) q[1];
sx q[1];
rz(-0.088610709) q[1];
sx q[1];
rz(-2.7068702) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22866727) q[3];
sx q[3];
rz(-2.1211984) q[3];
sx q[3];
rz(-2.500734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.053190319) q[2];
sx q[2];
rz(-1.3206427) q[2];
sx q[2];
rz(2.1957446) q[2];
rz(2.451402) q[3];
sx q[3];
rz(-1.918957) q[3];
sx q[3];
rz(1.3010196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.239554) q[0];
sx q[0];
rz(-1.5382465) q[0];
sx q[0];
rz(2.4334346) q[0];
rz(3.0802849) q[1];
sx q[1];
rz(-2.5262482) q[1];
sx q[1];
rz(1.6940438) q[1];
rz(-1.4881143) q[2];
sx q[2];
rz(-1.8607685) q[2];
sx q[2];
rz(0.49989732) q[2];
rz(1.7508909) q[3];
sx q[3];
rz(-1.4994499) q[3];
sx q[3];
rz(-2.3105619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
