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
rz(0.1306611) q[0];
sx q[0];
rz(4.9424439) q[0];
sx q[0];
rz(10.688936) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(-1.0882508) q[1];
sx q[1];
rz(0.01297125) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1931827) q[0];
sx q[0];
rz(-1.7275066) q[0];
sx q[0];
rz(2.66314) q[0];
x q[1];
rz(-0.58836909) q[2];
sx q[2];
rz(-0.19489637) q[2];
sx q[2];
rz(-2.4590059) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88718677) q[1];
sx q[1];
rz(-1.2206843) q[1];
sx q[1];
rz(-1.2302047) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.784934) q[3];
sx q[3];
rz(-0.99839568) q[3];
sx q[3];
rz(0.97808394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8791447) q[2];
sx q[2];
rz(-1.4538572) q[2];
sx q[2];
rz(0.095414735) q[2];
rz(0.48424193) q[3];
sx q[3];
rz(-2.3384194) q[3];
sx q[3];
rz(-1.4580844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478304) q[0];
sx q[0];
rz(-1.4365124) q[0];
sx q[0];
rz(-2.4875212) q[0];
rz(-1.7895128) q[1];
sx q[1];
rz(-0.65579954) q[1];
sx q[1];
rz(3.0316614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7470972) q[0];
sx q[0];
rz(-1.5506859) q[0];
sx q[0];
rz(-0.012152541) q[0];
rz(-pi) q[1];
rz(2.201033) q[2];
sx q[2];
rz(-1.9498683) q[2];
sx q[2];
rz(2.4565168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1812352) q[1];
sx q[1];
rz(-1.3981888) q[1];
sx q[1];
rz(-0.45876518) q[1];
x q[2];
rz(-0.55628784) q[3];
sx q[3];
rz(-2.6835052) q[3];
sx q[3];
rz(-2.213495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8249417) q[2];
sx q[2];
rz(-1.9330838) q[2];
sx q[2];
rz(-1.4671154) q[2];
rz(-1.0351099) q[3];
sx q[3];
rz(-0.78101522) q[3];
sx q[3];
rz(-1.8234113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255945) q[0];
sx q[0];
rz(-0.23402973) q[0];
sx q[0];
rz(-0.79063928) q[0];
rz(0.74360338) q[1];
sx q[1];
rz(-1.5433106) q[1];
sx q[1];
rz(-2.5620983) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4369219) q[0];
sx q[0];
rz(-1.0365067) q[0];
sx q[0];
rz(-2.2493258) q[0];
rz(-pi) q[1];
rz(-2.2966301) q[2];
sx q[2];
rz(-2.7562592) q[2];
sx q[2];
rz(-2.8266738) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6404785) q[1];
sx q[1];
rz(-1.6081297) q[1];
sx q[1];
rz(-1.8289315) q[1];
rz(2.6121796) q[3];
sx q[3];
rz(-2.0017481) q[3];
sx q[3];
rz(1.4612518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4641331) q[2];
sx q[2];
rz(-2.3049998) q[2];
sx q[2];
rz(-2.7602688) q[2];
rz(-0.85121202) q[3];
sx q[3];
rz(-0.92654735) q[3];
sx q[3];
rz(0.73494953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6094991) q[0];
sx q[0];
rz(-0.61436009) q[0];
sx q[0];
rz(2.7771948) q[0];
rz(2.9413307) q[1];
sx q[1];
rz(-1.7297144) q[1];
sx q[1];
rz(3.0416378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44463377) q[0];
sx q[0];
rz(-1.4958188) q[0];
sx q[0];
rz(-0.020453756) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9995297) q[2];
sx q[2];
rz(-1.6661281) q[2];
sx q[2];
rz(-2.7470061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8122319) q[1];
sx q[1];
rz(-1.4782923) q[1];
sx q[1];
rz(-1.5483556) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1853946) q[3];
sx q[3];
rz(-1.0805939) q[3];
sx q[3];
rz(-2.4407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0692856) q[2];
sx q[2];
rz(-2.0733209) q[2];
sx q[2];
rz(0.11492534) q[2];
rz(2.0011486) q[3];
sx q[3];
rz(-1.2830257) q[3];
sx q[3];
rz(-0.97525245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2045778) q[0];
sx q[0];
rz(-2.1890722) q[0];
sx q[0];
rz(0.17247795) q[0];
rz(1.0109488) q[1];
sx q[1];
rz(-0.5609678) q[1];
sx q[1];
rz(0.15484658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2456937) q[0];
sx q[0];
rz(-1.2749199) q[0];
sx q[0];
rz(0.22801836) q[0];
rz(1.8932976) q[2];
sx q[2];
rz(-0.66993827) q[2];
sx q[2];
rz(-2.8084286) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.78945827) q[1];
sx q[1];
rz(-2.704014) q[1];
sx q[1];
rz(0.34559135) q[1];
rz(-pi) q[2];
rz(2.5845086) q[3];
sx q[3];
rz(-2.9324813) q[3];
sx q[3];
rz(1.4534674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0871206) q[2];
sx q[2];
rz(-2.8643769) q[2];
sx q[2];
rz(-2.7692914) q[2];
rz(0.39792684) q[3];
sx q[3];
rz(-1.9979265) q[3];
sx q[3];
rz(-0.00042375617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.953124) q[0];
sx q[0];
rz(-0.57250452) q[0];
sx q[0];
rz(-1.5484126) q[0];
rz(-0.36987034) q[1];
sx q[1];
rz(-1.3984171) q[1];
sx q[1];
rz(0.63873783) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2681779) q[0];
sx q[0];
rz(-1.333433) q[0];
sx q[0];
rz(0.15837146) q[0];
rz(-pi) q[1];
rz(-1.0453141) q[2];
sx q[2];
rz(-1.3945701) q[2];
sx q[2];
rz(-2.3228563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.418701) q[1];
sx q[1];
rz(-1.1851553) q[1];
sx q[1];
rz(1.0807178) q[1];
rz(-2.1356167) q[3];
sx q[3];
rz(-2.2788439) q[3];
sx q[3];
rz(-2.4410409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0729735) q[2];
sx q[2];
rz(-2.0899453) q[2];
sx q[2];
rz(2.8886786) q[2];
rz(0.84826338) q[3];
sx q[3];
rz(-1.4784644) q[3];
sx q[3];
rz(0.084913582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0406168) q[0];
sx q[0];
rz(-2.9315797) q[0];
sx q[0];
rz(-2.9926391) q[0];
rz(1.3111929) q[1];
sx q[1];
rz(-1.4102178) q[1];
sx q[1];
rz(3.0874918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9210756) q[0];
sx q[0];
rz(-1.0722677) q[0];
sx q[0];
rz(1.6706628) q[0];
rz(-1.615754) q[2];
sx q[2];
rz(-2.3727594) q[2];
sx q[2];
rz(2.2368778) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91840345) q[1];
sx q[1];
rz(-1.4610054) q[1];
sx q[1];
rz(0.70007433) q[1];
rz(0.0232969) q[3];
sx q[3];
rz(-2.0459649) q[3];
sx q[3];
rz(-0.12926973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3397303) q[2];
sx q[2];
rz(-1.170155) q[2];
sx q[2];
rz(-2.1745963) q[2];
rz(-2.4407834) q[3];
sx q[3];
rz(-0.98027027) q[3];
sx q[3];
rz(0.35082671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0130149) q[0];
sx q[0];
rz(-1.9676493) q[0];
sx q[0];
rz(-1.5482192) q[0];
rz(-0.37823996) q[1];
sx q[1];
rz(-2.389237) q[1];
sx q[1];
rz(-2.7517448) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10791099) q[0];
sx q[0];
rz(-1.0709239) q[0];
sx q[0];
rz(-1.7938016) q[0];
rz(0.016640113) q[2];
sx q[2];
rz(-1.2667873) q[2];
sx q[2];
rz(-1.9250411) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9842062) q[1];
sx q[1];
rz(-2.0573924) q[1];
sx q[1];
rz(1.1729878) q[1];
rz(-pi) q[2];
rz(2.8865783) q[3];
sx q[3];
rz(-1.7863635) q[3];
sx q[3];
rz(0.6620342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0924015) q[2];
sx q[2];
rz(-1.1213877) q[2];
sx q[2];
rz(1.258705) q[2];
rz(-0.24426584) q[3];
sx q[3];
rz(-1.3552856) q[3];
sx q[3];
rz(-2.145483) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8213537) q[0];
sx q[0];
rz(-1.2293674) q[0];
sx q[0];
rz(-1.356333) q[0];
rz(-0.16607302) q[1];
sx q[1];
rz(-0.95172721) q[1];
sx q[1];
rz(0.43620268) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2552196) q[0];
sx q[0];
rz(-0.56111911) q[0];
sx q[0];
rz(-2.5291689) q[0];
rz(-pi) q[1];
rz(2.8032254) q[2];
sx q[2];
rz(-1.1778579) q[2];
sx q[2];
rz(1.5054877) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0862941) q[1];
sx q[1];
rz(-0.82189593) q[1];
sx q[1];
rz(-2.8248252) q[1];
x q[2];
rz(1.5028788) q[3];
sx q[3];
rz(-2.006586) q[3];
sx q[3];
rz(0.21623789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0210375) q[2];
sx q[2];
rz(-1.1641116) q[2];
sx q[2];
rz(-0.54350054) q[2];
rz(3.0920658) q[3];
sx q[3];
rz(-1.6560358) q[3];
sx q[3];
rz(-1.653695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5307584) q[0];
sx q[0];
rz(-3.0093091) q[0];
sx q[0];
rz(0.079205967) q[0];
rz(2.7414956) q[1];
sx q[1];
rz(-0.85405093) q[1];
sx q[1];
rz(-1.0265464) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29591376) q[0];
sx q[0];
rz(-0.85703731) q[0];
sx q[0];
rz(-2.8299324) q[0];
rz(0.24377771) q[2];
sx q[2];
rz(-0.37636435) q[2];
sx q[2];
rz(-1.7610904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.77791926) q[1];
sx q[1];
rz(-1.247974) q[1];
sx q[1];
rz(-1.670091) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21422106) q[3];
sx q[3];
rz(-2.7832) q[3];
sx q[3];
rz(-0.58979366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2058699) q[2];
sx q[2];
rz(-1.3047855) q[2];
sx q[2];
rz(1.4523466) q[2];
rz(-2.3116889) q[3];
sx q[3];
rz(-2.2645576) q[3];
sx q[3];
rz(-0.95480603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5230539) q[0];
sx q[0];
rz(-1.1445615) q[0];
sx q[0];
rz(-1.6339697) q[0];
rz(-1.2670831) q[1];
sx q[1];
rz(-1.0825842) q[1];
sx q[1];
rz(-2.4172197) q[1];
rz(0.83618589) q[2];
sx q[2];
rz(-0.7291353) q[2];
sx q[2];
rz(1.8020204) q[2];
rz(2.5821154) q[3];
sx q[3];
rz(-1.8636129) q[3];
sx q[3];
rz(-1.158798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
