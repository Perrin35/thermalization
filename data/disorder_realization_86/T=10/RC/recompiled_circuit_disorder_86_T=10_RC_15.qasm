OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(1.4120742) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(1.6593978) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5858551) q[0];
sx q[0];
rz(-1.1475539) q[0];
sx q[0];
rz(1.4002355) q[0];
x q[1];
rz(-2.8843845) q[2];
sx q[2];
rz(-1.4423941) q[2];
sx q[2];
rz(-0.53127015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9020302) q[1];
sx q[1];
rz(-0.67443496) q[1];
sx q[1];
rz(-2.2873138) q[1];
rz(-pi) q[2];
rz(-1.0079908) q[3];
sx q[3];
rz(-1.4853012) q[3];
sx q[3];
rz(-2.8443955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1564864) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(-2.1872897) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1433379) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(0.026219333) q[0];
rz(-1.6014618) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(2.1781133) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022097691) q[0];
sx q[0];
rz(-0.95675981) q[0];
sx q[0];
rz(-3.1387781) q[0];
rz(-pi) q[1];
rz(-2.1396779) q[2];
sx q[2];
rz(-1.2895673) q[2];
sx q[2];
rz(2.1525454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3669489) q[1];
sx q[1];
rz(-2.0683214) q[1];
sx q[1];
rz(2.7760387) q[1];
rz(1.9217334) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(1.2082781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5144689) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(0.13452402) q[2];
rz(-0.7450122) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298252) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(-2.0939317) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(2.581596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6277498) q[0];
sx q[0];
rz(-1.7729513) q[0];
sx q[0];
rz(-1.1802799) q[0];
x q[1];
rz(1.5494924) q[2];
sx q[2];
rz(-1.2676123) q[2];
sx q[2];
rz(-1.4959178) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8078976) q[1];
sx q[1];
rz(-1.9296608) q[1];
sx q[1];
rz(1.3586504) q[1];
x q[2];
rz(-1.0873763) q[3];
sx q[3];
rz(-1.1421575) q[3];
sx q[3];
rz(-0.40070686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3893163) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(2.1595188) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383706) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(0.28451434) q[0];
rz(-2.8248887) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(-1.8428615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7553058) q[0];
sx q[0];
rz(-0.55314976) q[0];
sx q[0];
rz(1.132071) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3372075) q[2];
sx q[2];
rz(-1.569869) q[2];
sx q[2];
rz(1.550012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.270077) q[1];
sx q[1];
rz(-1.3844826) q[1];
sx q[1];
rz(-0.99888505) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11233791) q[3];
sx q[3];
rz(-1.2445407) q[3];
sx q[3];
rz(2.5374967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(-2.7569125) q[2];
rz(-0.7540594) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383179) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(2.9105913) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(-2.8447661) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9223698) q[0];
sx q[0];
rz(-2.5693359) q[0];
sx q[0];
rz(-0.072364307) q[0];
rz(-pi) q[1];
rz(-1.7237687) q[2];
sx q[2];
rz(-2.3077871) q[2];
sx q[2];
rz(-1.0066777) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.27767147) q[1];
sx q[1];
rz(-1.1754981) q[1];
sx q[1];
rz(-2.0780095) q[1];
rz(-1.6597219) q[3];
sx q[3];
rz(-0.85907912) q[3];
sx q[3];
rz(2.5951648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(-0.39247593) q[2];
rz(1.9893507) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257989) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(0.75138599) q[0];
rz(-1.3279351) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(0.60633916) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7327001) q[0];
sx q[0];
rz(-3.0928287) q[0];
sx q[0];
rz(-0.34838895) q[0];
x q[1];
rz(-2.2491127) q[2];
sx q[2];
rz(-1.9024897) q[2];
sx q[2];
rz(1.734037) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50748435) q[1];
sx q[1];
rz(-1.023523) q[1];
sx q[1];
rz(2.9299111) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.074999853) q[3];
sx q[3];
rz(-1.3538176) q[3];
sx q[3];
rz(2.0892339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63885826) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(2.0022557) q[2];
rz(-1.6566488) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.5320324) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(-1.5628901) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(-0.79024822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3944053) q[0];
sx q[0];
rz(-2.1571775) q[0];
sx q[0];
rz(-1.2140843) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9279518) q[2];
sx q[2];
rz(-1.5593312) q[2];
sx q[2];
rz(-1.849888) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.33752791) q[1];
sx q[1];
rz(-1.3009326) q[1];
sx q[1];
rz(2.7359664) q[1];
x q[2];
rz(2.6584133) q[3];
sx q[3];
rz(-0.70671591) q[3];
sx q[3];
rz(0.28373517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.053085176) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(-0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4247894) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(2.0943663) q[0];
rz(0.60910243) q[1];
sx q[1];
rz(-1.4139688) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6435218) q[0];
sx q[0];
rz(-2.5812491) q[0];
sx q[0];
rz(-1.6269006) q[0];
x q[1];
rz(-1.1461166) q[2];
sx q[2];
rz(-1.0815902) q[2];
sx q[2];
rz(-2.3236772) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70506239) q[1];
sx q[1];
rz(-0.092518004) q[1];
sx q[1];
rz(0.1235048) q[1];
rz(0.29626131) q[3];
sx q[3];
rz(-0.98516609) q[3];
sx q[3];
rz(-2.0010009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9528815) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(2.2593373) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(-1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700478) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(1.4260938) q[0];
rz(-3.0601314) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-2.5833599) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4687913) q[0];
sx q[0];
rz(-1.1835915) q[0];
sx q[0];
rz(2.0331435) q[0];
rz(1.105537) q[2];
sx q[2];
rz(-3.0214587) q[2];
sx q[2];
rz(0.58819729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0135632) q[1];
sx q[1];
rz(-2.3131144) q[1];
sx q[1];
rz(-3.0548884) q[1];
rz(-pi) q[2];
rz(1.7562808) q[3];
sx q[3];
rz(-1.0914601) q[3];
sx q[3];
rz(2.6244147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2925064) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(-0.21197453) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(-1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.50487173) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.6037534) q[0];
rz(2.3161855) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(2.6182981) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9742045) q[0];
sx q[0];
rz(-2.5295969) q[0];
sx q[0];
rz(0.1083072) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9988437) q[2];
sx q[2];
rz(-1.5339282) q[2];
sx q[2];
rz(1.8429304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0303505) q[1];
sx q[1];
rz(-2.6551464) q[1];
sx q[1];
rz(-2.4187947) q[1];
rz(-0.2449805) q[3];
sx q[3];
rz(-1.3462726) q[3];
sx q[3];
rz(0.45104879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.837073) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(2.8722897) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(1.0860898) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3257278) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(-1.6124484) q[2];
sx q[2];
rz(-2.8786082) q[2];
sx q[2];
rz(-1.0515121) q[2];
rz(-0.23444093) q[3];
sx q[3];
rz(-2.3307878) q[3];
sx q[3];
rz(2.5651991) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
