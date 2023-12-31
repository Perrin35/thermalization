OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7944613) q[0];
sx q[0];
rz(-2.1262655) q[0];
sx q[0];
rz(2.6740958) q[0];
rz(2.6214018) q[1];
sx q[1];
rz(-1.7953035) q[1];
sx q[1];
rz(-2.4612114) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677833) q[0];
sx q[0];
rz(-0.73497226) q[0];
sx q[0];
rz(1.3761139) q[0];
rz(-2.1490287) q[2];
sx q[2];
rz(-0.34878584) q[2];
sx q[2];
rz(2.0860096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2441794) q[1];
sx q[1];
rz(-1.3379828) q[1];
sx q[1];
rz(-2.0884872) q[1];
rz(-0.47547961) q[3];
sx q[3];
rz(-2.7694422) q[3];
sx q[3];
rz(1.1572157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15158571) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(0.087466784) q[2];
rz(0.72922373) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.4366348) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(1.0013642) q[0];
rz(0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(-0.52406812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090416106) q[0];
sx q[0];
rz(-0.91958445) q[0];
sx q[0];
rz(-0.37474664) q[0];
rz(-pi) q[1];
rz(-0.87191138) q[2];
sx q[2];
rz(-1.9053835) q[2];
sx q[2];
rz(1.3427693) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3327649) q[1];
sx q[1];
rz(-1.1953224) q[1];
sx q[1];
rz(-1.8217702) q[1];
x q[2];
rz(0.25604053) q[3];
sx q[3];
rz(-2.5898993) q[3];
sx q[3];
rz(1.8191169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9841763) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(-1.3400419) q[2];
rz(0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(-3.1047344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3443417) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(0.62477338) q[0];
rz(-2.1773188) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(0.18951167) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4288085) q[0];
sx q[0];
rz(-1.841294) q[0];
sx q[0];
rz(2.4012647) q[0];
x q[1];
rz(1.6763716) q[2];
sx q[2];
rz(-1.9441248) q[2];
sx q[2];
rz(2.2764652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3526926) q[1];
sx q[1];
rz(-2.8955724) q[1];
sx q[1];
rz(0.23060631) q[1];
rz(-pi) q[2];
rz(0.17574163) q[3];
sx q[3];
rz(-0.97676859) q[3];
sx q[3];
rz(0.16080805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0345962) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(1.754388) q[2];
rz(-2.710279) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(-1.0881902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4819734) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(-1.5455998) q[0];
rz(2.003147) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(-1.0901573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2290303) q[0];
sx q[0];
rz(-1.1704485) q[0];
sx q[0];
rz(1.9489261) q[0];
rz(-pi) q[1];
rz(0.81972576) q[2];
sx q[2];
rz(-1.1851386) q[2];
sx q[2];
rz(-2.0476598) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62413299) q[1];
sx q[1];
rz(-1.0003261) q[1];
sx q[1];
rz(-1.073451) q[1];
rz(-pi) q[2];
rz(-1.3642717) q[3];
sx q[3];
rz(-1.1518475) q[3];
sx q[3];
rz(1.2057613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99889341) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(1.1936197) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(2.2284171) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(-2.1504452) q[1];
sx q[1];
rz(-1.9350355) q[1];
sx q[1];
rz(-1.9715462) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8613375) q[0];
sx q[0];
rz(-0.72685188) q[0];
sx q[0];
rz(1.0759541) q[0];
rz(2.2511475) q[2];
sx q[2];
rz(-1.8257942) q[2];
sx q[2];
rz(1.6971708) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3586732) q[1];
sx q[1];
rz(-1.1766608) q[1];
sx q[1];
rz(-2.6637117) q[1];
rz(0.53797651) q[3];
sx q[3];
rz(-1.7378983) q[3];
sx q[3];
rz(2.3218384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2698007) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(-0.13723792) q[2];
rz(-1.8042701) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18734922) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(1.7911918) q[0];
rz(-2.2987135) q[1];
sx q[1];
rz(-2.402669) q[1];
sx q[1];
rz(0.67289105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85149375) q[0];
sx q[0];
rz(-2.5585472) q[0];
sx q[0];
rz(-1.3821938) q[0];
x q[1];
rz(1.5422836) q[2];
sx q[2];
rz(-2.4270504) q[2];
sx q[2];
rz(1.5039832) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91025483) q[1];
sx q[1];
rz(-2.5006223) q[1];
sx q[1];
rz(-2.7892968) q[1];
rz(-pi) q[2];
rz(-1.6478496) q[3];
sx q[3];
rz(-1.1601845) q[3];
sx q[3];
rz(-1.3239469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.792753) q[2];
sx q[2];
rz(-1.2092084) q[2];
sx q[2];
rz(0.43506452) q[2];
rz(-1.7815636) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.9687987) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(2.0794179) q[0];
rz(-2.0299714) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(1.4020845) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5324425) q[0];
sx q[0];
rz(-2.6770868) q[0];
sx q[0];
rz(0.64446489) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24537556) q[2];
sx q[2];
rz(-2.3615395) q[2];
sx q[2];
rz(1.9666372) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84343108) q[1];
sx q[1];
rz(-1.2097675) q[1];
sx q[1];
rz(1.4069188) q[1];
rz(-pi) q[2];
rz(-1.7756895) q[3];
sx q[3];
rz(-2.5435102) q[3];
sx q[3];
rz(0.3860592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7193675) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(-0.68022234) q[2];
rz(-2.7133572) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(-1.359882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.632804) q[0];
sx q[0];
rz(-1.4161685) q[0];
sx q[0];
rz(-1.592214) q[0];
rz(-0.24958615) q[1];
sx q[1];
rz(-1.1523749) q[1];
sx q[1];
rz(-2.6002398) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2490059) q[0];
sx q[0];
rz(-1.5595946) q[0];
sx q[0];
rz(0.5641933) q[0];
rz(-pi) q[1];
rz(-0.45221046) q[2];
sx q[2];
rz(-1.8794249) q[2];
sx q[2];
rz(-1.2436777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7975446) q[1];
sx q[1];
rz(-0.86595264) q[1];
sx q[1];
rz(-2.123358) q[1];
x q[2];
rz(-1.5726611) q[3];
sx q[3];
rz(-1.1904753) q[3];
sx q[3];
rz(0.62455458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0970739) q[2];
sx q[2];
rz(-0.35733435) q[2];
sx q[2];
rz(-0.77735916) q[2];
rz(2.2903531) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.3210993) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38354307) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(2.126157) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(-1.6962956) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8953067) q[0];
sx q[0];
rz(-2.200897) q[0];
sx q[0];
rz(-0.90755983) q[0];
x q[1];
rz(2.5601013) q[2];
sx q[2];
rz(-0.98052374) q[2];
sx q[2];
rz(-1.3587388) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6362793) q[1];
sx q[1];
rz(-1.5552102) q[1];
sx q[1];
rz(-0.63025766) q[1];
x q[2];
rz(-2.2441838) q[3];
sx q[3];
rz(-1.0703147) q[3];
sx q[3];
rz(-2.9874779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69892591) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(0.60738579) q[2];
rz(1.7025042) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33567515) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(-0.39224958) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0948254) q[0];
sx q[0];
rz(-1.1724768) q[0];
sx q[0];
rz(0.70478583) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62551542) q[2];
sx q[2];
rz(-2.054347) q[2];
sx q[2];
rz(-3.0160883) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1019198) q[1];
sx q[1];
rz(-0.79636785) q[1];
sx q[1];
rz(2.2069195) q[1];
rz(1.5418105) q[3];
sx q[3];
rz(-1.3658804) q[3];
sx q[3];
rz(-0.25587413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.067651) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(-2.5411141) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(1.0935121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(-2.1451163) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-0.36841064) q[2];
sx q[2];
rz(-2.340292) q[2];
sx q[2];
rz(3.1030263) q[2];
rz(-0.92580933) q[3];
sx q[3];
rz(-2.3498597) q[3];
sx q[3];
rz(1.3510977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
