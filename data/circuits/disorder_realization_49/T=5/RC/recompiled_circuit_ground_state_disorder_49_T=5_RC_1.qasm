OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.16728084) q[0];
sx q[0];
rz(3.4790223) q[0];
sx q[0];
rz(12.364388) q[0];
rz(2.1972411) q[1];
sx q[1];
rz(-0.75466067) q[1];
sx q[1];
rz(1.4442297) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9044735) q[0];
sx q[0];
rz(-2.5982862) q[0];
sx q[0];
rz(0.41443698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7770257) q[2];
sx q[2];
rz(-1.4650868) q[2];
sx q[2];
rz(-1.4240698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3814427) q[1];
sx q[1];
rz(-2.0408148) q[1];
sx q[1];
rz(-2.2862741) q[1];
rz(-pi) q[2];
rz(1.4524649) q[3];
sx q[3];
rz(-2.7744996) q[3];
sx q[3];
rz(1.5637507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8586388) q[2];
sx q[2];
rz(-3.0573461) q[2];
sx q[2];
rz(-1.6072744) q[2];
rz(-0.18680799) q[3];
sx q[3];
rz(-0.45972937) q[3];
sx q[3];
rz(-1.648858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1780136) q[0];
sx q[0];
rz(-0.78129617) q[0];
sx q[0];
rz(0.70469967) q[0];
rz(1.0164545) q[1];
sx q[1];
rz(-1.1926788) q[1];
sx q[1];
rz(-1.1526398) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1267075) q[0];
sx q[0];
rz(-2.5683157) q[0];
sx q[0];
rz(1.2529729) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75147475) q[2];
sx q[2];
rz(-2.412999) q[2];
sx q[2];
rz(2.1214243) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1981467) q[1];
sx q[1];
rz(-1.0662765) q[1];
sx q[1];
rz(0.90934244) q[1];
x q[2];
rz(-0.14039881) q[3];
sx q[3];
rz(-0.75274668) q[3];
sx q[3];
rz(-1.5206159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6013201) q[2];
sx q[2];
rz(-0.86869621) q[2];
sx q[2];
rz(1.8216088) q[2];
rz(-0.071062239) q[3];
sx q[3];
rz(-0.080852121) q[3];
sx q[3];
rz(-2.0474056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0367301) q[0];
sx q[0];
rz(-2.7920089) q[0];
sx q[0];
rz(-0.87376755) q[0];
rz(-1.8050487) q[1];
sx q[1];
rz(-2.4599894) q[1];
sx q[1];
rz(-0.85493404) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23531518) q[0];
sx q[0];
rz(-1.3585257) q[0];
sx q[0];
rz(-2.2215823) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94742355) q[2];
sx q[2];
rz(-0.50259841) q[2];
sx q[2];
rz(-2.0386774) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1347772) q[1];
sx q[1];
rz(-1.2689277) q[1];
sx q[1];
rz(1.4567489) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.072238) q[3];
sx q[3];
rz(-0.66594687) q[3];
sx q[3];
rz(1.7518821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2130337) q[2];
sx q[2];
rz(-1.1308257) q[2];
sx q[2];
rz(1.0553168) q[2];
rz(-2.1999551) q[3];
sx q[3];
rz(-1.0712653) q[3];
sx q[3];
rz(0.94282237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68614352) q[0];
sx q[0];
rz(-2.2930155) q[0];
sx q[0];
rz(-0.69547478) q[0];
rz(-1.4675568) q[1];
sx q[1];
rz(-1.8753884) q[1];
sx q[1];
rz(-3.0470336) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7692663) q[0];
sx q[0];
rz(-2.8594198) q[0];
sx q[0];
rz(-1.397294) q[0];
x q[1];
rz(-2.5369625) q[2];
sx q[2];
rz(-2.1994091) q[2];
sx q[2];
rz(0.18456799) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2219005) q[1];
sx q[1];
rz(-1.2893493) q[1];
sx q[1];
rz(0.20604659) q[1];
rz(0.88282449) q[3];
sx q[3];
rz(-1.3677639) q[3];
sx q[3];
rz(-0.95850755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62119421) q[2];
sx q[2];
rz(-0.79402557) q[2];
sx q[2];
rz(2.3801079) q[2];
rz(0.17101184) q[3];
sx q[3];
rz(-1.2820425) q[3];
sx q[3];
rz(0.096693501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8347725) q[0];
sx q[0];
rz(-0.67245317) q[0];
sx q[0];
rz(0.49224725) q[0];
rz(1.3193839) q[1];
sx q[1];
rz(-1.4232114) q[1];
sx q[1];
rz(2.6548903) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0224441) q[0];
sx q[0];
rz(-2.7239425) q[0];
sx q[0];
rz(-0.29830427) q[0];
rz(-2.3721061) q[2];
sx q[2];
rz(-1.8459919) q[2];
sx q[2];
rz(-1.9857032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0472187) q[1];
sx q[1];
rz(-1.9969121) q[1];
sx q[1];
rz(-2.7414118) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0858889) q[3];
sx q[3];
rz(-1.3660226) q[3];
sx q[3];
rz(1.0398454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95185602) q[2];
sx q[2];
rz(-0.8437914) q[2];
sx q[2];
rz(-1.0978511) q[2];
rz(2.3073933) q[3];
sx q[3];
rz(-2.4662377) q[3];
sx q[3];
rz(-1.437291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51682353) q[0];
sx q[0];
rz(-0.56684816) q[0];
sx q[0];
rz(-0.045850642) q[0];
rz(0.74371964) q[1];
sx q[1];
rz(-0.36879483) q[1];
sx q[1];
rz(-3.0812982) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9550388) q[0];
sx q[0];
rz(-2.5738724) q[0];
sx q[0];
rz(-0.7132775) q[0];
rz(0.39581516) q[2];
sx q[2];
rz(-0.91378736) q[2];
sx q[2];
rz(2.3558877) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0452779) q[1];
sx q[1];
rz(-1.244925) q[1];
sx q[1];
rz(-0.23479825) q[1];
x q[2];
rz(-1.4521818) q[3];
sx q[3];
rz(-0.68202344) q[3];
sx q[3];
rz(2.8220146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9149949) q[2];
sx q[2];
rz(-1.8151585) q[2];
sx q[2];
rz(0.9787406) q[2];
rz(-1.2616875) q[3];
sx q[3];
rz(-1.0190957) q[3];
sx q[3];
rz(0.88920465) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9642445) q[0];
sx q[0];
rz(-1.6151936) q[0];
sx q[0];
rz(-0.88441315) q[0];
rz(-0.85909596) q[1];
sx q[1];
rz(-2.0180549) q[1];
sx q[1];
rz(0.20800796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73586845) q[0];
sx q[0];
rz(-2.5130773) q[0];
sx q[0];
rz(-2.5121538) q[0];
rz(-2.538717) q[2];
sx q[2];
rz(-1.920421) q[2];
sx q[2];
rz(0.58793758) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1237202) q[1];
sx q[1];
rz(-1.8353288) q[1];
sx q[1];
rz(-2.2295647) q[1];
rz(-1.0674123) q[3];
sx q[3];
rz(-1.5725822) q[3];
sx q[3];
rz(-1.3942277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.38006833) q[2];
sx q[2];
rz(-2.2218349) q[2];
sx q[2];
rz(2.5999787) q[2];
rz(-1.8875061) q[3];
sx q[3];
rz(-1.5897635) q[3];
sx q[3];
rz(0.90887535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33614531) q[0];
sx q[0];
rz(-2.8271524) q[0];
sx q[0];
rz(1.0910777) q[0];
rz(0.82487851) q[1];
sx q[1];
rz(-0.42366091) q[1];
sx q[1];
rz(-1.0728015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7110853) q[0];
sx q[0];
rz(-2.6737464) q[0];
sx q[0];
rz(-2.7070295) q[0];
rz(-pi) q[1];
rz(-0.013117803) q[2];
sx q[2];
rz(-1.5934391) q[2];
sx q[2];
rz(-0.72005872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3923918) q[1];
sx q[1];
rz(-1.8129686) q[1];
sx q[1];
rz(1.7081038) q[1];
x q[2];
rz(-1.9735224) q[3];
sx q[3];
rz(-1.3478876) q[3];
sx q[3];
rz(2.9532578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23952809) q[2];
sx q[2];
rz(-0.54486474) q[2];
sx q[2];
rz(-2.876335) q[2];
rz(0.15051633) q[3];
sx q[3];
rz(-2.033332) q[3];
sx q[3];
rz(-0.36293852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6657669) q[0];
sx q[0];
rz(-0.095223991) q[0];
sx q[0];
rz(-1.4775403) q[0];
rz(-2.3382969) q[1];
sx q[1];
rz(-1.7684312) q[1];
sx q[1];
rz(2.1243748) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8724339) q[0];
sx q[0];
rz(-0.24155051) q[0];
sx q[0];
rz(-2.1369324) q[0];
rz(-pi) q[1];
x q[1];
rz(2.656865) q[2];
sx q[2];
rz(-0.51186168) q[2];
sx q[2];
rz(0.32527015) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0086144) q[1];
sx q[1];
rz(-0.45615754) q[1];
sx q[1];
rz(2.4309979) q[1];
rz(-pi) q[2];
rz(2.3961762) q[3];
sx q[3];
rz(-0.43161296) q[3];
sx q[3];
rz(-3.0109757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21641714) q[2];
sx q[2];
rz(-1.5745796) q[2];
sx q[2];
rz(2.4851921) q[2];
rz(-2.9260855) q[3];
sx q[3];
rz(-1.7301205) q[3];
sx q[3];
rz(1.6925252) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0553174) q[0];
sx q[0];
rz(-1.3441939) q[0];
sx q[0];
rz(2.2917746) q[0];
rz(0.96018106) q[1];
sx q[1];
rz(-0.4147059) q[1];
sx q[1];
rz(-0.93920952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612027) q[0];
sx q[0];
rz(-1.8341258) q[0];
sx q[0];
rz(2.3376746) q[0];
rz(0.0077590436) q[2];
sx q[2];
rz(-1.9853704) q[2];
sx q[2];
rz(1.7637635) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43329217) q[1];
sx q[1];
rz(-0.95656779) q[1];
sx q[1];
rz(-2.0310165) q[1];
rz(-0.47795145) q[3];
sx q[3];
rz(-2.047973) q[3];
sx q[3];
rz(-1.7750224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5592929) q[2];
sx q[2];
rz(-1.7377995) q[2];
sx q[2];
rz(0.55029184) q[2];
rz(-0.127921) q[3];
sx q[3];
rz(-3.0030799) q[3];
sx q[3];
rz(0.96854717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.6561103) q[0];
sx q[0];
rz(-0.67018296) q[0];
sx q[0];
rz(-0.67076587) q[0];
rz(0.35519629) q[1];
sx q[1];
rz(-1.6658446) q[1];
sx q[1];
rz(2.7459941) q[1];
rz(-0.14329362) q[2];
sx q[2];
rz(-1.431965) q[2];
sx q[2];
rz(1.9964249) q[2];
rz(-1.1199748) q[3];
sx q[3];
rz(-0.49334855) q[3];
sx q[3];
rz(1.7176513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
