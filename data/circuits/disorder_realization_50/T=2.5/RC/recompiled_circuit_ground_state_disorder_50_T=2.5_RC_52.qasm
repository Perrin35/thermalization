OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49837056) q[0];
sx q[0];
rz(4.8084365) q[0];
sx q[0];
rz(9.2118535) q[0];
rz(-2.9990745) q[1];
sx q[1];
rz(-0.63894874) q[1];
sx q[1];
rz(-2.6899333) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8925326) q[0];
sx q[0];
rz(-2.0304168) q[0];
sx q[0];
rz(-1.100774) q[0];
x q[1];
rz(-0.75889905) q[2];
sx q[2];
rz(-0.90803185) q[2];
sx q[2];
rz(0.77292216) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2131552) q[1];
sx q[1];
rz(-0.31107956) q[1];
sx q[1];
rz(-2.3401287) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0805243) q[3];
sx q[3];
rz(-2.3754915) q[3];
sx q[3];
rz(1.776772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7529922) q[2];
sx q[2];
rz(-1.708932) q[2];
sx q[2];
rz(2.5373552) q[2];
rz(-2.9271017) q[3];
sx q[3];
rz(-2.4356804) q[3];
sx q[3];
rz(0.96116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251505) q[0];
sx q[0];
rz(-2.3264628) q[0];
sx q[0];
rz(2.2863638) q[0];
rz(1.1280355) q[1];
sx q[1];
rz(-2.3927092) q[1];
sx q[1];
rz(1.6181207) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.010941) q[0];
sx q[0];
rz(-2.391577) q[0];
sx q[0];
rz(-0.86096835) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5775213) q[2];
sx q[2];
rz(-0.84779948) q[2];
sx q[2];
rz(-1.3644219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2115626) q[1];
sx q[1];
rz(-2.6324144) q[1];
sx q[1];
rz(0.58873621) q[1];
rz(-1.4030966) q[3];
sx q[3];
rz(-2.0032231) q[3];
sx q[3];
rz(0.54673115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.35260674) q[2];
sx q[2];
rz(-1.2615729) q[2];
sx q[2];
rz(-1.2471586) q[2];
rz(-0.17036197) q[3];
sx q[3];
rz(-0.43694654) q[3];
sx q[3];
rz(-0.40444571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929297) q[0];
sx q[0];
rz(-0.12032838) q[0];
sx q[0];
rz(-0.89068252) q[0];
rz(2.0590674) q[1];
sx q[1];
rz(-0.71290103) q[1];
sx q[1];
rz(-3.069186) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88806498) q[0];
sx q[0];
rz(-1.0529336) q[0];
sx q[0];
rz(-0.47674556) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35066013) q[2];
sx q[2];
rz(-2.0267916) q[2];
sx q[2];
rz(-0.58619546) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5634656) q[1];
sx q[1];
rz(-1.153128) q[1];
sx q[1];
rz(2.595957) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4506574) q[3];
sx q[3];
rz(-1.9124219) q[3];
sx q[3];
rz(-1.6610314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.34042865) q[2];
sx q[2];
rz(-1.195793) q[2];
sx q[2];
rz(2.8677531) q[2];
rz(-1.0770816) q[3];
sx q[3];
rz(-1.2696973) q[3];
sx q[3];
rz(0.70335189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7766137) q[0];
sx q[0];
rz(-0.7319428) q[0];
sx q[0];
rz(-1.8044385) q[0];
rz(2.8058167) q[1];
sx q[1];
rz(-1.2715205) q[1];
sx q[1];
rz(1.550386) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44061892) q[0];
sx q[0];
rz(-1.0583726) q[0];
sx q[0];
rz(-0.62823624) q[0];
x q[1];
rz(-0.50652047) q[2];
sx q[2];
rz(-1.532785) q[2];
sx q[2];
rz(-0.7061298) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76564255) q[1];
sx q[1];
rz(-0.11184622) q[1];
sx q[1];
rz(1.2601529) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1071572) q[3];
sx q[3];
rz(-2.2175466) q[3];
sx q[3];
rz(1.3440901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35440272) q[2];
sx q[2];
rz(-1.7305817) q[2];
sx q[2];
rz(1.2350941) q[2];
rz(2.6835594) q[3];
sx q[3];
rz(-2.1221275) q[3];
sx q[3];
rz(2.4052446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6422727) q[0];
sx q[0];
rz(-0.90773931) q[0];
sx q[0];
rz(0.15550144) q[0];
rz(0.78089619) q[1];
sx q[1];
rz(-0.28762329) q[1];
sx q[1];
rz(-0.23316613) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1564395) q[0];
sx q[0];
rz(-1.4228983) q[0];
sx q[0];
rz(-0.83507706) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67552394) q[2];
sx q[2];
rz(-1.6758989) q[2];
sx q[2];
rz(1.6750248) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0079769) q[1];
sx q[1];
rz(-1.4583499) q[1];
sx q[1];
rz(1.8419209) q[1];
rz(-2.1114757) q[3];
sx q[3];
rz(-1.0305163) q[3];
sx q[3];
rz(-2.212611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95722014) q[2];
sx q[2];
rz(-1.5868712) q[2];
sx q[2];
rz(-2.8389944) q[2];
rz(-0.042081632) q[3];
sx q[3];
rz(-2.849597) q[3];
sx q[3];
rz(0.80872768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0724723) q[0];
sx q[0];
rz(-1.1061677) q[0];
sx q[0];
rz(-2.1492667) q[0];
rz(0.26556695) q[1];
sx q[1];
rz(-2.3873603) q[1];
sx q[1];
rz(-3.1297562) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965746) q[0];
sx q[0];
rz(-1.2906605) q[0];
sx q[0];
rz(-0.046124553) q[0];
rz(-pi) q[1];
rz(1.0008079) q[2];
sx q[2];
rz(-1.3234183) q[2];
sx q[2];
rz(-0.22817366) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.34764987) q[1];
sx q[1];
rz(-1.2041438) q[1];
sx q[1];
rz(-0.071392752) q[1];
rz(-pi) q[2];
rz(-1.8675667) q[3];
sx q[3];
rz(-2.9057876) q[3];
sx q[3];
rz(-2.2267226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5806233) q[2];
sx q[2];
rz(-2.2443266) q[2];
sx q[2];
rz(0.90274367) q[2];
rz(0.48815253) q[3];
sx q[3];
rz(-1.3267696) q[3];
sx q[3];
rz(-1.3443525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.97486597) q[0];
sx q[0];
rz(-0.45477295) q[0];
sx q[0];
rz(-1.3653261) q[0];
rz(-2.4873554) q[1];
sx q[1];
rz(-0.74834329) q[1];
sx q[1];
rz(-1.0985589) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7099534) q[0];
sx q[0];
rz(-2.0281784) q[0];
sx q[0];
rz(-2.4541992) q[0];
x q[1];
rz(-2.9113414) q[2];
sx q[2];
rz(-2.8881713) q[2];
sx q[2];
rz(2.3524323) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0579018) q[1];
sx q[1];
rz(-2.5595287) q[1];
sx q[1];
rz(-2.9681724) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1992616) q[3];
sx q[3];
rz(-1.7876995) q[3];
sx q[3];
rz(-2.1700117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72413033) q[2];
sx q[2];
rz(-1.3516358) q[2];
sx q[2];
rz(-0.89428085) q[2];
rz(-1.2096679) q[3];
sx q[3];
rz(-3.0318048) q[3];
sx q[3];
rz(-0.82718682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2073479) q[0];
sx q[0];
rz(-1.3149911) q[0];
sx q[0];
rz(0.19499245) q[0];
rz(0.64942819) q[1];
sx q[1];
rz(-2.1959627) q[1];
sx q[1];
rz(1.4494928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60842321) q[0];
sx q[0];
rz(-2.4618645) q[0];
sx q[0];
rz(2.5191368) q[0];
rz(-pi) q[1];
rz(-2.605841) q[2];
sx q[2];
rz(-1.881534) q[2];
sx q[2];
rz(-2.9536332) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7342775) q[1];
sx q[1];
rz(-2.4859634) q[1];
sx q[1];
rz(1.6796965) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9417848) q[3];
sx q[3];
rz(-2.2358353) q[3];
sx q[3];
rz(-0.68913078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0608757) q[2];
sx q[2];
rz(-2.1970811) q[2];
sx q[2];
rz(0.057849217) q[2];
rz(-0.079023376) q[3];
sx q[3];
rz(-1.9130324) q[3];
sx q[3];
rz(1.2270989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(2.3651315) q[0];
sx q[0];
rz(-1.6420028) q[0];
sx q[0];
rz(2.814433) q[0];
rz(0.51045927) q[1];
sx q[1];
rz(-1.3182498) q[1];
sx q[1];
rz(3.0697451) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7450388) q[0];
sx q[0];
rz(-1.7522194) q[0];
sx q[0];
rz(-1.4722162) q[0];
x q[1];
rz(-2.2141404) q[2];
sx q[2];
rz(-2.1610552) q[2];
sx q[2];
rz(0.037742712) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.10061) q[1];
sx q[1];
rz(-2.1400053) q[1];
sx q[1];
rz(-0.74633523) q[1];
rz(-pi) q[2];
rz(-1.184578) q[3];
sx q[3];
rz(-0.90208331) q[3];
sx q[3];
rz(0.45763256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0923126) q[2];
sx q[2];
rz(-1.6593554) q[2];
sx q[2];
rz(-1.2266889) q[2];
rz(-2.6022794) q[3];
sx q[3];
rz(-1.156811) q[3];
sx q[3];
rz(-1.4815319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6698089) q[0];
sx q[0];
rz(-1.1265378) q[0];
sx q[0];
rz(-0.54927611) q[0];
rz(1.1726941) q[1];
sx q[1];
rz(-1.6212308) q[1];
sx q[1];
rz(-1.8424013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8090813) q[0];
sx q[0];
rz(-1.375688) q[0];
sx q[0];
rz(1.5089351) q[0];
rz(-pi) q[1];
rz(-2.4626715) q[2];
sx q[2];
rz(-1.2679865) q[2];
sx q[2];
rz(-0.34679373) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0944984) q[1];
sx q[1];
rz(-0.21378042) q[1];
sx q[1];
rz(0.70014145) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2659508) q[3];
sx q[3];
rz(-1.048363) q[3];
sx q[3];
rz(-1.6506139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4486763) q[2];
sx q[2];
rz(-2.0557949) q[2];
sx q[2];
rz(0.14796251) q[2];
rz(-0.97370094) q[3];
sx q[3];
rz(-2.2805043) q[3];
sx q[3];
rz(2.2073943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9964518) q[0];
sx q[0];
rz(-1.8358163) q[0];
sx q[0];
rz(-1.3124574) q[0];
rz(1.3120069) q[1];
sx q[1];
rz(-0.94999718) q[1];
sx q[1];
rz(-2.2012262) q[1];
rz(1.1205705) q[2];
sx q[2];
rz(-0.2926338) q[2];
sx q[2];
rz(-0.50526239) q[2];
rz(2.9964126) q[3];
sx q[3];
rz(-2.1812781) q[3];
sx q[3];
rz(-1.3727544) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
