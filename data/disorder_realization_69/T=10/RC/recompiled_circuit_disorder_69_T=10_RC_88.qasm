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
rz(-0.52019083) q[1];
sx q[1];
rz(1.7953035) q[1];
sx q[1];
rz(10.105159) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2517393) q[0];
sx q[0];
rz(-1.4407053) q[0];
sx q[0];
rz(-2.2962909) q[0];
x q[1];
rz(-2.9453968) q[2];
sx q[2];
rz(-1.2805403) q[2];
sx q[2];
rz(-0.44858518) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2441794) q[1];
sx q[1];
rz(-1.3379828) q[1];
sx q[1];
rz(-1.0531055) q[1];
rz(-pi) q[2];
rz(-1.3939788) q[3];
sx q[3];
rz(-1.2415981) q[3];
sx q[3];
rz(-0.65229177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15158571) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(0.087466784) q[2];
rz(-0.72922373) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-2.1402284) q[0];
rz(0.17240605) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4755718) q[0];
sx q[0];
rz(-0.7374987) q[0];
sx q[0];
rz(1.1230418) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42627724) q[2];
sx q[2];
rz(-0.91765109) q[2];
sx q[2];
rz(-0.041235812) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9973035) q[1];
sx q[1];
rz(-1.3376437) q[1];
sx q[1];
rz(-0.38645978) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25604053) q[3];
sx q[3];
rz(-0.55169332) q[3];
sx q[3];
rz(1.3224758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9841763) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(-1.8015507) q[2];
rz(2.346758) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(-0.036858233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3443417) q[0];
sx q[0];
rz(-0.69673711) q[0];
sx q[0];
rz(0.62477338) q[0];
rz(2.1773188) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(-0.18951167) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42650578) q[0];
sx q[0];
rz(-2.3622892) q[0];
sx q[0];
rz(-2.7515609) q[0];
rz(1.465221) q[2];
sx q[2];
rz(-1.1974679) q[2];
sx q[2];
rz(-0.8651274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0057919766) q[1];
sx q[1];
rz(-1.5151007) q[1];
sx q[1];
rz(2.9018351) q[1];
rz(-pi) q[2];
rz(1.3174921) q[3];
sx q[3];
rz(-2.5251303) q[3];
sx q[3];
rz(-2.9951819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0345962) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(-1.3872046) q[2];
rz(-0.43131367) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(-2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4819734) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(1.5455998) q[0];
rz(-1.1384456) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(1.0901573) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2290303) q[0];
sx q[0];
rz(-1.1704485) q[0];
sx q[0];
rz(-1.1926665) q[0];
rz(0.81972576) q[2];
sx q[2];
rz(-1.956454) q[2];
sx q[2];
rz(-1.0939329) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.480099) q[1];
sx q[1];
rz(-1.9839994) q[1];
sx q[1];
rz(0.6306298) q[1];
x q[2];
rz(1.3642717) q[3];
sx q[3];
rz(-1.1518475) q[3];
sx q[3];
rz(-1.2057613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(-2.2857655) q[2];
rz(-1.1936197) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49044931) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-2.9918616) q[0];
rz(2.1504452) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.1700464) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67384185) q[0];
sx q[0];
rz(-1.2497328) q[0];
sx q[0];
rz(2.2348316) q[0];
rz(-pi) q[1];
rz(1.9636376) q[2];
sx q[2];
rz(-0.71937865) q[2];
sx q[2];
rz(0.17578416) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7829195) q[1];
sx q[1];
rz(-1.9649319) q[1];
sx q[1];
rz(0.47788099) q[1];
rz(-pi) q[2];
rz(0.53797651) q[3];
sx q[3];
rz(-1.7378983) q[3];
sx q[3];
rz(-0.81975421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(-0.13723792) q[2];
rz(1.3373226) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(-1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18734922) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(1.7911918) q[0];
rz(2.2987135) q[1];
sx q[1];
rz(-2.402669) q[1];
sx q[1];
rz(2.4687016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.580299) q[0];
sx q[0];
rz(-1.467388) q[0];
sx q[0];
rz(-0.99594492) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2851373) q[2];
sx q[2];
rz(-1.5894784) q[2];
sx q[2];
rz(3.0532388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3402965) q[1];
sx q[1];
rz(-0.9749037) q[1];
sx q[1];
rz(1.3188386) q[1];
x q[2];
rz(-1.6478496) q[3];
sx q[3];
rz(-1.9814081) q[3];
sx q[3];
rz(1.3239469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3488397) q[2];
sx q[2];
rz(-1.2092084) q[2];
sx q[2];
rz(0.43506452) q[2];
rz(-1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9687987) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(2.0794179) q[0];
rz(-1.1116213) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(1.7395082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37004334) q[0];
sx q[0];
rz(-1.84329) q[0];
sx q[0];
rz(-0.38099654) q[0];
rz(-pi) q[1];
rz(2.3767396) q[2];
sx q[2];
rz(-1.7424889) q[2];
sx q[2];
rz(0.21966759) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84343108) q[1];
sx q[1];
rz(-1.9318252) q[1];
sx q[1];
rz(1.7346738) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7756895) q[3];
sx q[3];
rz(-0.59808245) q[3];
sx q[3];
rz(-0.3860592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4222251) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(-0.68022234) q[2];
rz(2.7133572) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(1.359882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.632804) q[0];
sx q[0];
rz(-1.4161685) q[0];
sx q[0];
rz(1.592214) q[0];
rz(-0.24958615) q[1];
sx q[1];
rz(-1.1523749) q[1];
sx q[1];
rz(-2.6002398) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2490059) q[0];
sx q[0];
rz(-1.5595946) q[0];
sx q[0];
rz(-0.5641933) q[0];
rz(-pi) q[1];
rz(1.9114242) q[2];
sx q[2];
rz(-2.0001786) q[2];
sx q[2];
rz(0.18061772) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3440481) q[1];
sx q[1];
rz(-0.86595264) q[1];
sx q[1];
rz(1.0182347) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0046644966) q[3];
sx q[3];
rz(-0.38032535) q[3];
sx q[3];
rz(2.5120146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.044518746) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-0.77735916) q[2];
rz(-2.2903531) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(-1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38354307) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(3.1316277) q[0];
rz(1.0154356) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(-1.6962956) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89307071) q[0];
sx q[0];
rz(-2.0914441) q[0];
sx q[0];
rz(0.74670665) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89500918) q[2];
sx q[2];
rz(-1.0969321) q[2];
sx q[2];
rz(2.5788139) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6362793) q[1];
sx q[1];
rz(-1.5552102) q[1];
sx q[1];
rz(0.63025766) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61049283) q[3];
sx q[3];
rz(-2.1497512) q[3];
sx q[3];
rz(1.0510774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(2.5342069) q[2];
rz(1.4390885) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.8059175) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(-2.5277396) q[0];
rz(-2.0954258) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(2.7493431) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0961571) q[0];
sx q[0];
rz(-2.3491612) q[0];
sx q[0];
rz(-2.5655454) q[0];
x q[1];
rz(0.62551542) q[2];
sx q[2];
rz(-1.0872456) q[2];
sx q[2];
rz(3.0160883) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0543921) q[1];
sx q[1];
rz(-2.0093845) q[1];
sx q[1];
rz(-0.88263504) q[1];
x q[2];
rz(-2.936593) q[3];
sx q[3];
rz(-1.5424171) q[3];
sx q[3];
rz(-1.8207707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0739416) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(2.5411141) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7077211) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(0.99647635) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-2.773182) q[2];
sx q[2];
rz(-0.80130063) q[2];
sx q[2];
rz(-0.038566312) q[2];
rz(-2.2157833) q[3];
sx q[3];
rz(-0.79173294) q[3];
sx q[3];
rz(-1.790495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
