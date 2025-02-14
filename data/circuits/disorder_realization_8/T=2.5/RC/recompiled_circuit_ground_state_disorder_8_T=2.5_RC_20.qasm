OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72579223) q[0];
sx q[0];
rz(1.9771201) q[0];
sx q[0];
rz(10.676072) q[0];
rz(2.8407821) q[1];
sx q[1];
rz(-1.7164813) q[1];
sx q[1];
rz(-3.0256174) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.011256) q[0];
sx q[0];
rz(-2.7948871) q[0];
sx q[0];
rz(1.1954855) q[0];
rz(-2.9094522) q[2];
sx q[2];
rz(-2.0439185) q[2];
sx q[2];
rz(-2.6363475) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.310461) q[1];
sx q[1];
rz(-0.77378213) q[1];
sx q[1];
rz(-1.6602064) q[1];
rz(-pi) q[2];
rz(-1.7099781) q[3];
sx q[3];
rz(-2.109189) q[3];
sx q[3];
rz(2.0288426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3731602) q[2];
sx q[2];
rz(-3.0765522) q[2];
sx q[2];
rz(1.6981286) q[2];
rz(2.4453435) q[3];
sx q[3];
rz(-1.5158451) q[3];
sx q[3];
rz(-0.35370091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8732805) q[0];
sx q[0];
rz(-0.074957632) q[0];
sx q[0];
rz(0.36175501) q[0];
rz(-1.3842281) q[1];
sx q[1];
rz(-0.39159602) q[1];
sx q[1];
rz(2.1143544) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8105668) q[0];
sx q[0];
rz(-0.89052708) q[0];
sx q[0];
rz(-2.9085338) q[0];
x q[1];
rz(1.7387668) q[2];
sx q[2];
rz(-0.78064519) q[2];
sx q[2];
rz(-1.3813684) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4467795) q[1];
sx q[1];
rz(-2.4449592) q[1];
sx q[1];
rz(-1.9637061) q[1];
x q[2];
rz(0.35137041) q[3];
sx q[3];
rz(-2.9020269) q[3];
sx q[3];
rz(-0.75557652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.1343753) q[2];
sx q[2];
rz(-0.5642887) q[2];
sx q[2];
rz(-1.0312274) q[2];
rz(-0.68650308) q[3];
sx q[3];
rz(-1.406823) q[3];
sx q[3];
rz(-2.8197207) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6512076) q[0];
sx q[0];
rz(-1.4116766) q[0];
sx q[0];
rz(-1.7899293) q[0];
rz(-1.6199813) q[1];
sx q[1];
rz(-0.23623513) q[1];
sx q[1];
rz(1.9550386) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4410858) q[0];
sx q[0];
rz(-1.1916627) q[0];
sx q[0];
rz(0.4613045) q[0];
rz(-pi) q[1];
rz(2.988838) q[2];
sx q[2];
rz(-0.4604333) q[2];
sx q[2];
rz(0.82376119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12315519) q[1];
sx q[1];
rz(-1.5644524) q[1];
sx q[1];
rz(-0.0090339418) q[1];
x q[2];
rz(0.29680829) q[3];
sx q[3];
rz(-2.2389484) q[3];
sx q[3];
rz(0.53199024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97570193) q[2];
sx q[2];
rz(-1.2624792) q[2];
sx q[2];
rz(-0.28124896) q[2];
rz(2.7665372) q[3];
sx q[3];
rz(-2.227759) q[3];
sx q[3];
rz(-2.8121172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220045) q[0];
sx q[0];
rz(-0.81398886) q[0];
sx q[0];
rz(0.32459146) q[0];
rz(-1.5931256) q[1];
sx q[1];
rz(-2.8209768) q[1];
sx q[1];
rz(-2.2494907) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2587713) q[0];
sx q[0];
rz(-3.1207001) q[0];
sx q[0];
rz(-2.5740795) q[0];
x q[1];
rz(1.6692144) q[2];
sx q[2];
rz(-0.8618671) q[2];
sx q[2];
rz(0.44385982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2617787) q[1];
sx q[1];
rz(-0.38912762) q[1];
sx q[1];
rz(-2.3107982) q[1];
x q[2];
rz(1.5209429) q[3];
sx q[3];
rz(-0.59403235) q[3];
sx q[3];
rz(-1.0083053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0804245) q[2];
sx q[2];
rz(-1.3865546) q[2];
sx q[2];
rz(-0.57677734) q[2];
rz(3.0171677) q[3];
sx q[3];
rz(-0.7500698) q[3];
sx q[3];
rz(0.09593825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6782003) q[0];
sx q[0];
rz(-0.41825774) q[0];
sx q[0];
rz(0.28045714) q[0];
rz(-2.8434143) q[1];
sx q[1];
rz(-3.0715946) q[1];
sx q[1];
rz(-0.15204522) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66528377) q[0];
sx q[0];
rz(-1.5179993) q[0];
sx q[0];
rz(-3.1140181) q[0];
rz(-1.0481669) q[2];
sx q[2];
rz(-1.657147) q[2];
sx q[2];
rz(-0.62355838) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8271395) q[1];
sx q[1];
rz(-1.9681853) q[1];
sx q[1];
rz(1.1177776) q[1];
x q[2];
rz(-1.1880615) q[3];
sx q[3];
rz(-1.5627341) q[3];
sx q[3];
rz(-2.1915367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2966447) q[2];
sx q[2];
rz(-1.7903906) q[2];
sx q[2];
rz(-3.1015934) q[2];
rz(-3.0039039) q[3];
sx q[3];
rz(-0.44860336) q[3];
sx q[3];
rz(2.3562446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8765151) q[0];
sx q[0];
rz(-3.0410933) q[0];
sx q[0];
rz(-0.58078372) q[0];
rz(2.4526217) q[1];
sx q[1];
rz(-0.13801485) q[1];
sx q[1];
rz(-1.2977915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0397206) q[0];
sx q[0];
rz(-2.5301109) q[0];
sx q[0];
rz(-2.3518556) q[0];
rz(-pi) q[1];
rz(1.0270732) q[2];
sx q[2];
rz(-1.5158733) q[2];
sx q[2];
rz(1.3872272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7766469) q[1];
sx q[1];
rz(-0.45393911) q[1];
sx q[1];
rz(2.620082) q[1];
rz(2.5491015) q[3];
sx q[3];
rz(-1.1969229) q[3];
sx q[3];
rz(0.44523525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40655228) q[2];
sx q[2];
rz(-1.3383957) q[2];
sx q[2];
rz(1.5470541) q[2];
rz(-0.43153396) q[3];
sx q[3];
rz(-2.0412692) q[3];
sx q[3];
rz(0.6507473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4887017) q[0];
sx q[0];
rz(-3.1259013) q[0];
sx q[0];
rz(-2.5456862) q[0];
rz(-1.7911004) q[1];
sx q[1];
rz(-0.48521438) q[1];
sx q[1];
rz(-2.5686666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112537) q[0];
sx q[0];
rz(-1.5736394) q[0];
sx q[0];
rz(-0.009441998) q[0];
rz(-2.7155034) q[2];
sx q[2];
rz(-1.3182148) q[2];
sx q[2];
rz(-1.9468228) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0966037) q[1];
sx q[1];
rz(-1.2536238) q[1];
sx q[1];
rz(0.60896341) q[1];
x q[2];
rz(0.75667824) q[3];
sx q[3];
rz(-1.8493493) q[3];
sx q[3];
rz(-3.0708282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4493745) q[2];
sx q[2];
rz(-0.91392475) q[2];
sx q[2];
rz(-2.2268028) q[2];
rz(-1.5335013) q[3];
sx q[3];
rz(-2.0017109) q[3];
sx q[3];
rz(-0.99005121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6221301) q[0];
sx q[0];
rz(-1.7161481) q[0];
sx q[0];
rz(-3.07716) q[0];
rz(-0.26647767) q[1];
sx q[1];
rz(-0.12424145) q[1];
sx q[1];
rz(-1.2640094) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2867133) q[0];
sx q[0];
rz(-2.7164796) q[0];
sx q[0];
rz(-0.69010587) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6904628) q[2];
sx q[2];
rz(-2.3887803) q[2];
sx q[2];
rz(3.0835763) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19169584) q[1];
sx q[1];
rz(-2.765871) q[1];
sx q[1];
rz(-2.5671178) q[1];
rz(-pi) q[2];
rz(-1.0059935) q[3];
sx q[3];
rz(-2.0801615) q[3];
sx q[3];
rz(-2.6107174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5233351) q[2];
sx q[2];
rz(-1.9674415) q[2];
sx q[2];
rz(1.0443643) q[2];
rz(2.4339645) q[3];
sx q[3];
rz(-2.7510567) q[3];
sx q[3];
rz(-1.0467168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7758961) q[0];
sx q[0];
rz(-1.5713659) q[0];
sx q[0];
rz(2.7112992) q[0];
rz(-2.4039092) q[1];
sx q[1];
rz(-3.0569515) q[1];
sx q[1];
rz(0.35668361) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9707613) q[0];
sx q[0];
rz(-2.8342842) q[0];
sx q[0];
rz(-0.31913646) q[0];
x q[1];
rz(-1.4595152) q[2];
sx q[2];
rz(-0.26646895) q[2];
sx q[2];
rz(-2.1499718) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0903527) q[1];
sx q[1];
rz(-1.151671) q[1];
sx q[1];
rz(2.5263949) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0218871) q[3];
sx q[3];
rz(-1.6448391) q[3];
sx q[3];
rz(0.56379217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3714527) q[2];
sx q[2];
rz(-2.3944103) q[2];
sx q[2];
rz(2.746197) q[2];
rz(-1.4622408) q[3];
sx q[3];
rz(-1.7489) q[3];
sx q[3];
rz(-0.028954884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8545561) q[0];
sx q[0];
rz(-0.77571464) q[0];
sx q[0];
rz(-2.4391644) q[0];
rz(-2.9632945) q[1];
sx q[1];
rz(-2.4874004) q[1];
sx q[1];
rz(1.5231232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.130647) q[0];
sx q[0];
rz(-0.48640448) q[0];
sx q[0];
rz(-2.0964699) q[0];
x q[1];
rz(-1.7179397) q[2];
sx q[2];
rz(-3.0426263) q[2];
sx q[2];
rz(0.56848016) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9257912) q[1];
sx q[1];
rz(-0.58432441) q[1];
sx q[1];
rz(2.8846442) q[1];
x q[2];
rz(-2.5036292) q[3];
sx q[3];
rz(-2.5416982) q[3];
sx q[3];
rz(-2.0323582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6452597) q[2];
sx q[2];
rz(-0.96480227) q[2];
sx q[2];
rz(0.94950914) q[2];
rz(1.1358787) q[3];
sx q[3];
rz(-2.1954326) q[3];
sx q[3];
rz(-2.5778263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62354325) q[0];
sx q[0];
rz(-1.8907056) q[0];
sx q[0];
rz(-0.84778669) q[0];
rz(1.4354979) q[1];
sx q[1];
rz(-1.2789627) q[1];
sx q[1];
rz(0.36437558) q[1];
rz(2.2709104) q[2];
sx q[2];
rz(-1.489594) q[2];
sx q[2];
rz(1.7832047) q[2];
rz(-0.061312231) q[3];
sx q[3];
rz(-1.8777962) q[3];
sx q[3];
rz(2.5971436) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
