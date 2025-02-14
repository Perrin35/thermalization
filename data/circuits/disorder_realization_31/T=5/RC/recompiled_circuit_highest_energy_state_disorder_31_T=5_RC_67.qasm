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
rz(1.2827058) q[0];
sx q[0];
rz(-2.1252706) q[0];
sx q[0];
rz(1.8860201) q[0];
rz(-1.3610871) q[1];
sx q[1];
rz(-0.95870107) q[1];
sx q[1];
rz(0.60671848) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17598303) q[0];
sx q[0];
rz(-2.2874766) q[0];
sx q[0];
rz(-0.59881439) q[0];
rz(-3.1090281) q[2];
sx q[2];
rz(-1.3059907) q[2];
sx q[2];
rz(2.3015568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1944988) q[1];
sx q[1];
rz(-1.441212) q[1];
sx q[1];
rz(1.2034077) q[1];
rz(-pi) q[2];
rz(2.9806271) q[3];
sx q[3];
rz(-1.4518629) q[3];
sx q[3];
rz(1.87784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2029734) q[2];
sx q[2];
rz(-1.9846658) q[2];
sx q[2];
rz(2.970676) q[2];
rz(-1.1275229) q[3];
sx q[3];
rz(-2.5850962) q[3];
sx q[3];
rz(-0.053430406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.663986) q[0];
sx q[0];
rz(-0.32653102) q[0];
sx q[0];
rz(-2.4531181) q[0];
rz(-1.7636048) q[1];
sx q[1];
rz(-1.0208027) q[1];
sx q[1];
rz(3.0000906) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3162075) q[0];
sx q[0];
rz(-2.0658464) q[0];
sx q[0];
rz(-2.0284257) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42948855) q[2];
sx q[2];
rz(-1.6613243) q[2];
sx q[2];
rz(-1.2538101) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0857755) q[1];
sx q[1];
rz(-1.2165099) q[1];
sx q[1];
rz(-1.2694759) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3745295) q[3];
sx q[3];
rz(-2.7500023) q[3];
sx q[3];
rz(-2.151769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3839533) q[2];
sx q[2];
rz(-0.44718224) q[2];
sx q[2];
rz(-3.1128856) q[2];
rz(-0.38255295) q[3];
sx q[3];
rz(-1.6304723) q[3];
sx q[3];
rz(-2.9401275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2460227) q[0];
sx q[0];
rz(-2.0992794) q[0];
sx q[0];
rz(2.7699455) q[0];
rz(2.9314575) q[1];
sx q[1];
rz(-0.34966436) q[1];
sx q[1];
rz(-1.9336611) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4233404) q[0];
sx q[0];
rz(-2.1499942) q[0];
sx q[0];
rz(0.96286003) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0591427) q[2];
sx q[2];
rz(-1.5918613) q[2];
sx q[2];
rz(2.9631071) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74151285) q[1];
sx q[1];
rz(-2.1472405) q[1];
sx q[1];
rz(-0.6735354) q[1];
x q[2];
rz(-1.6116834) q[3];
sx q[3];
rz(-1.2256116) q[3];
sx q[3];
rz(-0.57905771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5189884) q[2];
sx q[2];
rz(-0.092978803) q[2];
sx q[2];
rz(0.63817111) q[2];
rz(-0.41247621) q[3];
sx q[3];
rz(-1.6715489) q[3];
sx q[3];
rz(1.1023869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07688044) q[0];
sx q[0];
rz(-2.2314254) q[0];
sx q[0];
rz(1.3077211) q[0];
rz(-2.4619596) q[1];
sx q[1];
rz(-2.4145587) q[1];
sx q[1];
rz(1.1839428) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1550876) q[0];
sx q[0];
rz(-2.3851624) q[0];
sx q[0];
rz(2.6593701) q[0];
rz(2.9058564) q[2];
sx q[2];
rz(-1.4582468) q[2];
sx q[2];
rz(-1.7692009) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1278504) q[1];
sx q[1];
rz(-2.7847185) q[1];
sx q[1];
rz(2.9452717) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3700962) q[3];
sx q[3];
rz(-2.960254) q[3];
sx q[3];
rz(-0.69252959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1188941) q[2];
sx q[2];
rz(-0.1476295) q[2];
sx q[2];
rz(2.0856608) q[2];
rz(2.8583156) q[3];
sx q[3];
rz(-1.2669468) q[3];
sx q[3];
rz(1.7578846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0530171) q[0];
sx q[0];
rz(-2.179189) q[0];
sx q[0];
rz(1.4085294) q[0];
rz(-2.9119496) q[1];
sx q[1];
rz(-2.2848928) q[1];
sx q[1];
rz(-1.1913258) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2076805) q[0];
sx q[0];
rz(-1.4290853) q[0];
sx q[0];
rz(-1.1674387) q[0];
rz(-1.2651839) q[2];
sx q[2];
rz(-1.4508045) q[2];
sx q[2];
rz(-3.0141399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7745668) q[1];
sx q[1];
rz(-1.3210249) q[1];
sx q[1];
rz(-0.43065177) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91893371) q[3];
sx q[3];
rz(-1.100538) q[3];
sx q[3];
rz(1.7277499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64629897) q[2];
sx q[2];
rz(-2.7806492) q[2];
sx q[2];
rz(1.6634644) q[2];
rz(2.9804001) q[3];
sx q[3];
rz(-1.5321923) q[3];
sx q[3];
rz(0.78981367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5021055) q[0];
sx q[0];
rz(-2.7171071) q[0];
sx q[0];
rz(2.2232527) q[0];
rz(-2.8847252) q[1];
sx q[1];
rz(-0.99595064) q[1];
sx q[1];
rz(1.1161944) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6300121) q[0];
sx q[0];
rz(-1.0718379) q[0];
sx q[0];
rz(-0.60006882) q[0];
rz(-2.8681197) q[2];
sx q[2];
rz(-1.6533706) q[2];
sx q[2];
rz(1.7805748) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5815365) q[1];
sx q[1];
rz(-1.3764842) q[1];
sx q[1];
rz(0.24340731) q[1];
rz(2.8982387) q[3];
sx q[3];
rz(-1.4019308) q[3];
sx q[3];
rz(2.4143743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18672289) q[2];
sx q[2];
rz(-2.335151) q[2];
sx q[2];
rz(0.40697971) q[2];
rz(-0.61378971) q[3];
sx q[3];
rz(-1.5301306) q[3];
sx q[3];
rz(2.6500402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75519049) q[0];
sx q[0];
rz(-0.84119216) q[0];
sx q[0];
rz(0.92591539) q[0];
rz(2.6689802) q[1];
sx q[1];
rz(-2.0874529) q[1];
sx q[1];
rz(2.6714163) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7177954) q[0];
sx q[0];
rz(-1.6542302) q[0];
sx q[0];
rz(-0.83441894) q[0];
rz(-pi) q[1];
rz(2.2165259) q[2];
sx q[2];
rz(-1.888141) q[2];
sx q[2];
rz(0.44337526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8308816) q[1];
sx q[1];
rz(-2.355) q[1];
sx q[1];
rz(1.5194233) q[1];
x q[2];
rz(-1.6117134) q[3];
sx q[3];
rz(-2.4114263) q[3];
sx q[3];
rz(-2.1388291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4818695) q[2];
sx q[2];
rz(-2.0707776) q[2];
sx q[2];
rz(-2.3107963) q[2];
rz(0.080549084) q[3];
sx q[3];
rz(-2.060067) q[3];
sx q[3];
rz(2.1848333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.047121) q[0];
sx q[0];
rz(-1.6166649) q[0];
sx q[0];
rz(2.9834874) q[0];
rz(0.72599167) q[1];
sx q[1];
rz(-2.7114365) q[1];
sx q[1];
rz(-0.37011883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86240101) q[0];
sx q[0];
rz(-1.4790311) q[0];
sx q[0];
rz(0.33359762) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49741715) q[2];
sx q[2];
rz(-0.82225613) q[2];
sx q[2];
rz(-0.17360877) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.41851317) q[1];
sx q[1];
rz(-1.1795768) q[1];
sx q[1];
rz(1.9926975) q[1];
x q[2];
rz(-2.1158989) q[3];
sx q[3];
rz(-1.2152367) q[3];
sx q[3];
rz(-0.73425435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98081723) q[2];
sx q[2];
rz(-1.2951916) q[2];
sx q[2];
rz(3.004461) q[2];
rz(1.6454227) q[3];
sx q[3];
rz(-2.4479595) q[3];
sx q[3];
rz(-2.6579198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5982323) q[0];
sx q[0];
rz(-1.0022663) q[0];
sx q[0];
rz(-3.015236) q[0];
rz(-0.57451808) q[1];
sx q[1];
rz(-2.2210329) q[1];
sx q[1];
rz(-0.7472907) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68796038) q[0];
sx q[0];
rz(-1.9575141) q[0];
sx q[0];
rz(-0.20129119) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4805675) q[2];
sx q[2];
rz(-2.3213561) q[2];
sx q[2];
rz(-2.837452) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0382749) q[1];
sx q[1];
rz(-0.59477287) q[1];
sx q[1];
rz(-0.86869356) q[1];
x q[2];
rz(-2.6797557) q[3];
sx q[3];
rz(-0.29134068) q[3];
sx q[3];
rz(-0.1801462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1025461) q[2];
sx q[2];
rz(-2.0206082) q[2];
sx q[2];
rz(-0.54753629) q[2];
rz(-1.840379) q[3];
sx q[3];
rz(-2.0042714) q[3];
sx q[3];
rz(-0.12714061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8319156) q[0];
sx q[0];
rz(-1.7597821) q[0];
sx q[0];
rz(0.58050138) q[0];
rz(2.8864587) q[1];
sx q[1];
rz(-2.3105123) q[1];
sx q[1];
rz(-1.9560122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1398753) q[0];
sx q[0];
rz(-0.96065646) q[0];
sx q[0];
rz(-2.3185179) q[0];
x q[1];
rz(-1.5747936) q[2];
sx q[2];
rz(-2.7786715) q[2];
sx q[2];
rz(-1.6238348) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.45764582) q[1];
sx q[1];
rz(-2.6630029) q[1];
sx q[1];
rz(-2.9475888) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85082741) q[3];
sx q[3];
rz(-1.730207) q[3];
sx q[3];
rz(3.0992143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4127976) q[2];
sx q[2];
rz(-0.86270657) q[2];
sx q[2];
rz(-2.1346788) q[2];
rz(1.6179196) q[3];
sx q[3];
rz(-0.98098743) q[3];
sx q[3];
rz(1.1716918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11378743) q[0];
sx q[0];
rz(-1.2600949) q[0];
sx q[0];
rz(-2.8366198) q[0];
rz(1.2420568) q[1];
sx q[1];
rz(-1.2868953) q[1];
sx q[1];
rz(-2.3931265) q[1];
rz(1.1390252) q[2];
sx q[2];
rz(-1.4130269) q[2];
sx q[2];
rz(-2.4561981) q[2];
rz(2.2119568) q[3];
sx q[3];
rz(-0.84992483) q[3];
sx q[3];
rz(0.58569943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
