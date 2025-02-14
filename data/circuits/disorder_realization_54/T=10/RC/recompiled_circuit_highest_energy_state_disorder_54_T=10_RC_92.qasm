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
rz(2.430001) q[0];
sx q[0];
rz(-1.1049668) q[0];
sx q[0];
rz(2.0018863) q[0];
rz(2.3501514) q[1];
sx q[1];
rz(-2.1005519) q[1];
sx q[1];
rz(1.1295553) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4576042) q[0];
sx q[0];
rz(-0.86375551) q[0];
sx q[0];
rz(-2.9024505) q[0];
rz(-pi) q[1];
rz(1.5989675) q[2];
sx q[2];
rz(-0.90785387) q[2];
sx q[2];
rz(0.65671317) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0186386) q[1];
sx q[1];
rz(-1.8644511) q[1];
sx q[1];
rz(2.3112626) q[1];
x q[2];
rz(3.0409052) q[3];
sx q[3];
rz(-0.33439454) q[3];
sx q[3];
rz(1.7582126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0396042) q[2];
sx q[2];
rz(-2.3309989) q[2];
sx q[2];
rz(-0.90768901) q[2];
rz(-3.0229819) q[3];
sx q[3];
rz(-1.6018931) q[3];
sx q[3];
rz(2.8906726) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1949961) q[0];
sx q[0];
rz(-2.4793766) q[0];
sx q[0];
rz(-3.0372341) q[0];
rz(1.1907499) q[1];
sx q[1];
rz(-1.933681) q[1];
sx q[1];
rz(0.45713919) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32923082) q[0];
sx q[0];
rz(-2.9456249) q[0];
sx q[0];
rz(1.4120483) q[0];
rz(2.1188583) q[2];
sx q[2];
rz(-0.75621683) q[2];
sx q[2];
rz(0.05847419) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7105661) q[1];
sx q[1];
rz(-1.885864) q[1];
sx q[1];
rz(-3.0418212) q[1];
rz(-pi) q[2];
rz(0.13808226) q[3];
sx q[3];
rz(-0.90681091) q[3];
sx q[3];
rz(-2.9241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5856058) q[2];
sx q[2];
rz(-1.1588691) q[2];
sx q[2];
rz(-0.15698329) q[2];
rz(-2.6202776) q[3];
sx q[3];
rz(-2.3020404) q[3];
sx q[3];
rz(3.1275911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74612015) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(-0.34573063) q[0];
rz(-1.455447) q[1];
sx q[1];
rz(-1.2724178) q[1];
sx q[1];
rz(-1.5706496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5356876) q[0];
sx q[0];
rz(-0.55258026) q[0];
sx q[0];
rz(-1.7406169) q[0];
rz(-0.94812553) q[2];
sx q[2];
rz(-1.5737163) q[2];
sx q[2];
rz(-2.7013283) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8040672) q[1];
sx q[1];
rz(-0.21906549) q[1];
sx q[1];
rz(0.22727025) q[1];
x q[2];
rz(2.9390252) q[3];
sx q[3];
rz(-0.93155471) q[3];
sx q[3];
rz(2.5646184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3674783) q[2];
sx q[2];
rz(-0.228129) q[2];
sx q[2];
rz(0.95743123) q[2];
rz(-1.1227603) q[3];
sx q[3];
rz(-1.2343531) q[3];
sx q[3];
rz(-1.7694337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289309) q[0];
sx q[0];
rz(-1.7191732) q[0];
sx q[0];
rz(-1.5439532) q[0];
rz(-1.7515901) q[1];
sx q[1];
rz(-1.3163687) q[1];
sx q[1];
rz(-1.4768627) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66946736) q[0];
sx q[0];
rz(-2.5335214) q[0];
sx q[0];
rz(0.50699349) q[0];
rz(-pi) q[1];
rz(3.0576287) q[2];
sx q[2];
rz(-1.8083085) q[2];
sx q[2];
rz(2.824914) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.37472607) q[1];
sx q[1];
rz(-1.0488759) q[1];
sx q[1];
rz(2.7609227) q[1];
rz(-1.847193) q[3];
sx q[3];
rz(-1.2528462) q[3];
sx q[3];
rz(-2.2459787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5433898) q[2];
sx q[2];
rz(-1.8566088) q[2];
sx q[2];
rz(0.99679917) q[2];
rz(-1.8761084) q[3];
sx q[3];
rz(-0.71804738) q[3];
sx q[3];
rz(0.27749458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.0604414) q[0];
sx q[0];
rz(-1.5716946) q[0];
sx q[0];
rz(1.0720217) q[0];
rz(-0.95442665) q[1];
sx q[1];
rz(-1.1773033) q[1];
sx q[1];
rz(1.6486453) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7537751) q[0];
sx q[0];
rz(-0.79083453) q[0];
sx q[0];
rz(2.8641939) q[0];
x q[1];
rz(0.49333879) q[2];
sx q[2];
rz(-0.78560053) q[2];
sx q[2];
rz(-2.0480905) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0860436) q[1];
sx q[1];
rz(-1.8978658) q[1];
sx q[1];
rz(-1.4202761) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29376438) q[3];
sx q[3];
rz(-0.33393327) q[3];
sx q[3];
rz(-1.5042083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.91168779) q[2];
sx q[2];
rz(-2.7547084) q[2];
sx q[2];
rz(1.9913199) q[2];
rz(-0.041042717) q[3];
sx q[3];
rz(-1.5951472) q[3];
sx q[3];
rz(0.40157792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3732442) q[0];
sx q[0];
rz(-2.0337489) q[0];
sx q[0];
rz(-1.9819697) q[0];
rz(1.4609963) q[1];
sx q[1];
rz(-0.20080876) q[1];
sx q[1];
rz(1.2332835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4077743) q[0];
sx q[0];
rz(-0.81627699) q[0];
sx q[0];
rz(0.4305779) q[0];
rz(-pi) q[1];
rz(-0.057282863) q[2];
sx q[2];
rz(-1.5174688) q[2];
sx q[2];
rz(1.48207) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94168866) q[1];
sx q[1];
rz(-2.2268128) q[1];
sx q[1];
rz(1.0628878) q[1];
rz(2.7090577) q[3];
sx q[3];
rz(-0.56946856) q[3];
sx q[3];
rz(-1.7908975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.636574) q[2];
sx q[2];
rz(-2.7754112) q[2];
sx q[2];
rz(1.150307) q[2];
rz(0.73927528) q[3];
sx q[3];
rz(-1.8940247) q[3];
sx q[3];
rz(-0.44481835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4742541) q[0];
sx q[0];
rz(-2.4485454) q[0];
sx q[0];
rz(0.44878238) q[0];
rz(-1.601864) q[1];
sx q[1];
rz(-1.1069143) q[1];
sx q[1];
rz(1.9805699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0036573) q[0];
sx q[0];
rz(-1.1627534) q[0];
sx q[0];
rz(-1.28027) q[0];
rz(-pi) q[1];
rz(0.49893219) q[2];
sx q[2];
rz(-2.009444) q[2];
sx q[2];
rz(0.89887757) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5827477) q[1];
sx q[1];
rz(-0.44932355) q[1];
sx q[1];
rz(-0.72968633) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0169898) q[3];
sx q[3];
rz(-1.3447059) q[3];
sx q[3];
rz(-2.6718321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1686958) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(1.1807582) q[2];
rz(-1.0817889) q[3];
sx q[3];
rz(-1.6094145) q[3];
sx q[3];
rz(-2.7150174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.54843724) q[0];
sx q[0];
rz(-1.3189545) q[0];
sx q[0];
rz(0.68761188) q[0];
rz(1.9718735) q[1];
sx q[1];
rz(-0.97427383) q[1];
sx q[1];
rz(2.7489472) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1448185) q[0];
sx q[0];
rz(-0.4985362) q[0];
sx q[0];
rz(1.2720435) q[0];
rz(1.7101426) q[2];
sx q[2];
rz(-1.2616488) q[2];
sx q[2];
rz(1.6935371) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.709014) q[1];
sx q[1];
rz(-0.48332941) q[1];
sx q[1];
rz(-2.9523115) q[1];
rz(-pi) q[2];
rz(1.6527376) q[3];
sx q[3];
rz(-1.382516) q[3];
sx q[3];
rz(-0.57037607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.16200599) q[2];
sx q[2];
rz(-1.0065099) q[2];
sx q[2];
rz(1.6458192) q[2];
rz(1.5227854) q[3];
sx q[3];
rz(-2.3706172) q[3];
sx q[3];
rz(0.60105598) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9252121) q[0];
sx q[0];
rz(-2.7583211) q[0];
sx q[0];
rz(1.8076757) q[0];
rz(-1.1434309) q[1];
sx q[1];
rz(-0.83088487) q[1];
sx q[1];
rz(-0.7935895) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0045354768) q[0];
sx q[0];
rz(-1.5001778) q[0];
sx q[0];
rz(2.9831992) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0751444) q[2];
sx q[2];
rz(-0.96198231) q[2];
sx q[2];
rz(0.43283909) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2379505) q[1];
sx q[1];
rz(-1.7667084) q[1];
sx q[1];
rz(-0.51301051) q[1];
x q[2];
rz(3.034044) q[3];
sx q[3];
rz(-2.0189773) q[3];
sx q[3];
rz(-1.5761689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6952343) q[2];
sx q[2];
rz(-2.0909205) q[2];
sx q[2];
rz(2.0900772) q[2];
rz(0.71600437) q[3];
sx q[3];
rz(-2.1456238) q[3];
sx q[3];
rz(-2.9314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7150772) q[0];
sx q[0];
rz(-2.2861013) q[0];
sx q[0];
rz(-2.9606384) q[0];
rz(-1.9521693) q[1];
sx q[1];
rz(-1.9633429) q[1];
sx q[1];
rz(2.1612371) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10889036) q[0];
sx q[0];
rz(-2.0712896) q[0];
sx q[0];
rz(-1.6653401) q[0];
x q[1];
rz(-2.3848214) q[2];
sx q[2];
rz(-2.5483661) q[2];
sx q[2];
rz(2.7293918) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7804523) q[1];
sx q[1];
rz(-1.3899511) q[1];
sx q[1];
rz(0.7493345) q[1];
x q[2];
rz(2.7465024) q[3];
sx q[3];
rz(-1.7655585) q[3];
sx q[3];
rz(-0.89941809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.02562) q[2];
sx q[2];
rz(-2.9774234) q[2];
sx q[2];
rz(1.8170961) q[2];
rz(0.65274158) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(-3.1033707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7272335) q[0];
sx q[0];
rz(-1.7200732) q[0];
sx q[0];
rz(2.1678069) q[0];
rz(-0.7242135) q[1];
sx q[1];
rz(-2.0661294) q[1];
sx q[1];
rz(-2.9758458) q[1];
rz(-1.235511) q[2];
sx q[2];
rz(-1.4852471) q[2];
sx q[2];
rz(0.5813364) q[2];
rz(-0.11304819) q[3];
sx q[3];
rz(-0.77533508) q[3];
sx q[3];
rz(-1.2067457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
