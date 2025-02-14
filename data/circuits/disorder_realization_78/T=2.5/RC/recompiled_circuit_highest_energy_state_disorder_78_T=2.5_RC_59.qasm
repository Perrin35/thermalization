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
rz(0.084963381) q[0];
sx q[0];
rz(-2.8391916) q[0];
sx q[0];
rz(3.1095355) q[0];
rz(3.1337466) q[1];
sx q[1];
rz(3.6454522) q[1];
sx q[1];
rz(7.03581) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40850484) q[0];
sx q[0];
rz(-1.2512387) q[0];
sx q[0];
rz(-1.3590786) q[0];
rz(-1.1967802) q[2];
sx q[2];
rz(-1.2311282) q[2];
sx q[2];
rz(2.938478) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.14559) q[1];
sx q[1];
rz(-1.8217345) q[1];
sx q[1];
rz(-0.35059021) q[1];
x q[2];
rz(1.2129477) q[3];
sx q[3];
rz(-1.7718302) q[3];
sx q[3];
rz(1.7261637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1746615) q[2];
sx q[2];
rz(-1.9661247) q[2];
sx q[2];
rz(2.3215129) q[2];
rz(-0.44102937) q[3];
sx q[3];
rz(-1.1038154) q[3];
sx q[3];
rz(2.5681514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.012861982) q[0];
sx q[0];
rz(-2.2282889) q[0];
sx q[0];
rz(0.41020694) q[0];
rz(-0.38093105) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(1.7832696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5429496) q[0];
sx q[0];
rz(-2.1457167) q[0];
sx q[0];
rz(-0.47358124) q[0];
x q[1];
rz(-2.3097281) q[2];
sx q[2];
rz(-0.68507776) q[2];
sx q[2];
rz(2.5545504) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4315003) q[1];
sx q[1];
rz(-2.9572801) q[1];
sx q[1];
rz(1.1460365) q[1];
rz(-2.210564) q[3];
sx q[3];
rz(-1.4461541) q[3];
sx q[3];
rz(-2.3384936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35839781) q[2];
sx q[2];
rz(-0.40893778) q[2];
sx q[2];
rz(-2.6386293) q[2];
rz(1.8424235) q[3];
sx q[3];
rz(-2.3830569) q[3];
sx q[3];
rz(-2.9097596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4101039) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(-0.9758392) q[0];
rz(-0.93636912) q[1];
sx q[1];
rz(-1.5543289) q[1];
sx q[1];
rz(3.0036614) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16743037) q[0];
sx q[0];
rz(-0.61744055) q[0];
sx q[0];
rz(-0.90888826) q[0];
x q[1];
rz(-1.4549655) q[2];
sx q[2];
rz(-1.9660913) q[2];
sx q[2];
rz(2.9651053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5974849) q[1];
sx q[1];
rz(-0.710809) q[1];
sx q[1];
rz(-1.7729458) q[1];
x q[2];
rz(0.50721747) q[3];
sx q[3];
rz(-2.6196369) q[3];
sx q[3];
rz(-1.0591266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58933538) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(-1.662558) q[2];
rz(0.79834437) q[3];
sx q[3];
rz(-2.2975497) q[3];
sx q[3];
rz(0.020309694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.289157) q[0];
sx q[0];
rz(-3.030179) q[0];
sx q[0];
rz(-1.0668466) q[0];
rz(1.5929219) q[1];
sx q[1];
rz(-1.8916062) q[1];
sx q[1];
rz(0.41935316) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8521312) q[0];
sx q[0];
rz(-1.3416107) q[0];
sx q[0];
rz(2.9782692) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52813645) q[2];
sx q[2];
rz(-1.4924876) q[2];
sx q[2];
rz(0.17158835) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.6904628) q[1];
sx q[1];
rz(-1.9277384) q[1];
sx q[1];
rz(0.068515645) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0556424) q[3];
sx q[3];
rz(-1.2618229) q[3];
sx q[3];
rz(0.072871836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35820094) q[2];
sx q[2];
rz(-0.50938598) q[2];
sx q[2];
rz(1.5979213) q[2];
rz(-1.915043) q[3];
sx q[3];
rz(-1.9251325) q[3];
sx q[3];
rz(1.290087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2655547) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(0.75620404) q[0];
rz(-1.2528231) q[1];
sx q[1];
rz(-0.93106657) q[1];
sx q[1];
rz(1.7255712) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.904793) q[0];
sx q[0];
rz(-2.5052245) q[0];
sx q[0];
rz(2.7704253) q[0];
x q[1];
rz(-0.3438832) q[2];
sx q[2];
rz(-1.6893708) q[2];
sx q[2];
rz(-0.32394513) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.071347) q[1];
sx q[1];
rz(-1.5637378) q[1];
sx q[1];
rz(0.53877212) q[1];
rz(-2.7783423) q[3];
sx q[3];
rz(-1.3923858) q[3];
sx q[3];
rz(1.6476064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8004134) q[2];
sx q[2];
rz(-2.4004553) q[2];
sx q[2];
rz(0.18079147) q[2];
rz(1.4888658) q[3];
sx q[3];
rz(-1.3988262) q[3];
sx q[3];
rz(-1.6256049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4329231) q[0];
sx q[0];
rz(-1.0754508) q[0];
sx q[0];
rz(2.3413626) q[0];
rz(-2.8569787) q[1];
sx q[1];
rz(-1.0601284) q[1];
sx q[1];
rz(-0.12399331) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6033961) q[0];
sx q[0];
rz(-1.5992224) q[0];
sx q[0];
rz(1.70065) q[0];
rz(-pi) q[1];
rz(3.0818865) q[2];
sx q[2];
rz(-1.8054031) q[2];
sx q[2];
rz(-1.6725181) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4912419) q[1];
sx q[1];
rz(-1.6524775) q[1];
sx q[1];
rz(2.8359069) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0923397) q[3];
sx q[3];
rz(-2.4774356) q[3];
sx q[3];
rz(-0.60598999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93512145) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(0.07621152) q[2];
rz(1.8501836) q[3];
sx q[3];
rz(-2.0468057) q[3];
sx q[3];
rz(2.5230303) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718219) q[0];
sx q[0];
rz(-1.5795647) q[0];
sx q[0];
rz(2.3845657) q[0];
rz(-2.3454759) q[1];
sx q[1];
rz(-1.4069822) q[1];
sx q[1];
rz(-2.1597791) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7863203) q[0];
sx q[0];
rz(-2.2844446) q[0];
sx q[0];
rz(2.7631604) q[0];
rz(-pi) q[1];
rz(0.070008833) q[2];
sx q[2];
rz(-1.4638136) q[2];
sx q[2];
rz(1.9221523) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.10985366) q[1];
sx q[1];
rz(-0.64424911) q[1];
sx q[1];
rz(2.2883313) q[1];
rz(-pi) q[2];
x q[2];
rz(2.997274) q[3];
sx q[3];
rz(-1.2929521) q[3];
sx q[3];
rz(1.2793737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7299399) q[2];
sx q[2];
rz(-0.81092683) q[2];
sx q[2];
rz(-2.5977503) q[2];
rz(-3.1206711) q[3];
sx q[3];
rz(-1.436751) q[3];
sx q[3];
rz(2.3232443) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92886096) q[0];
sx q[0];
rz(-2.2305771) q[0];
sx q[0];
rz(2.5182356) q[0];
rz(0.32132545) q[1];
sx q[1];
rz(-2.5265381) q[1];
sx q[1];
rz(0.26434937) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67548114) q[0];
sx q[0];
rz(-0.74680416) q[0];
sx q[0];
rz(1.8854499) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4521763) q[2];
sx q[2];
rz(-0.33782321) q[2];
sx q[2];
rz(1.092697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0496638) q[1];
sx q[1];
rz(-2.290002) q[1];
sx q[1];
rz(-0.37435162) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9569915) q[3];
sx q[3];
rz(-0.34117801) q[3];
sx q[3];
rz(-0.5285078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.136772) q[2];
sx q[2];
rz(-2.4592168) q[2];
sx q[2];
rz(0.47478673) q[2];
rz(0.36561203) q[3];
sx q[3];
rz(-2.6545299) q[3];
sx q[3];
rz(-0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30597618) q[0];
sx q[0];
rz(-0.37537471) q[0];
sx q[0];
rz(2.6557652) q[0];
rz(-1.0659418) q[1];
sx q[1];
rz(-1.7168047) q[1];
sx q[1];
rz(-2.8071383) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850824) q[0];
sx q[0];
rz(-1.872439) q[0];
sx q[0];
rz(1.7443875) q[0];
rz(-0.29269258) q[2];
sx q[2];
rz(-0.45537696) q[2];
sx q[2];
rz(-0.51560452) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85035861) q[1];
sx q[1];
rz(-2.6667074) q[1];
sx q[1];
rz(-0.33109457) q[1];
x q[2];
rz(2.0341088) q[3];
sx q[3];
rz(-2.3114738) q[3];
sx q[3];
rz(-0.23915672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7740384) q[2];
sx q[2];
rz(-0.84711051) q[2];
sx q[2];
rz(-0.96662194) q[2];
rz(1.6537846) q[3];
sx q[3];
rz(-2.8130468) q[3];
sx q[3];
rz(-3.0706792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3286572) q[0];
sx q[0];
rz(-0.85449496) q[0];
sx q[0];
rz(-0.27035126) q[0];
rz(1.6290172) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(-0.62634748) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.644076) q[0];
sx q[0];
rz(-0.77267161) q[0];
sx q[0];
rz(2.0324043) q[0];
rz(-pi) q[1];
rz(-2.2081991) q[2];
sx q[2];
rz(-1.2167756) q[2];
sx q[2];
rz(-2.6824981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.69520187) q[1];
sx q[1];
rz(-0.099041136) q[1];
sx q[1];
rz(0.42548577) q[1];
x q[2];
rz(-2.7553431) q[3];
sx q[3];
rz(-1.2674517) q[3];
sx q[3];
rz(-1.7084427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7865929) q[2];
sx q[2];
rz(-0.91675106) q[2];
sx q[2];
rz(1.5691441) q[2];
rz(2.3908424) q[3];
sx q[3];
rz(-2.7995977) q[3];
sx q[3];
rz(-0.87740889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56562051) q[0];
sx q[0];
rz(-1.2979869) q[0];
sx q[0];
rz(-0.57814231) q[0];
rz(-1.3427973) q[1];
sx q[1];
rz(-0.31675757) q[1];
sx q[1];
rz(-2.996179) q[1];
rz(0.12231355) q[2];
sx q[2];
rz(-1.5707471) q[2];
sx q[2];
rz(-1.6483501) q[2];
rz(-0.20082898) q[3];
sx q[3];
rz(-0.26293892) q[3];
sx q[3];
rz(0.50504167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
