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
rz(-0.11197055) q[0];
sx q[0];
rz(-2.114871) q[0];
sx q[0];
rz(1.8486899) q[0];
rz(-1.0678043) q[1];
sx q[1];
rz(-0.81967241) q[1];
sx q[1];
rz(-1.9424865) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1333251) q[0];
sx q[0];
rz(-0.12524167) q[0];
sx q[0];
rz(-0.52405806) q[0];
rz(-0.22817518) q[2];
sx q[2];
rz(-0.68636218) q[2];
sx q[2];
rz(0.31191269) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0173229) q[1];
sx q[1];
rz(-0.60397479) q[1];
sx q[1];
rz(1.9159957) q[1];
rz(-0.61618038) q[3];
sx q[3];
rz(-0.95762223) q[3];
sx q[3];
rz(2.4626436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5045515) q[2];
sx q[2];
rz(-1.7402288) q[2];
sx q[2];
rz(-1.441628) q[2];
rz(-1.0124892) q[3];
sx q[3];
rz(-1.9952521) q[3];
sx q[3];
rz(-2.6700524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3016894) q[0];
sx q[0];
rz(-2.58044) q[0];
sx q[0];
rz(1.0294234) q[0];
rz(1.3599716) q[1];
sx q[1];
rz(-1.9554892) q[1];
sx q[1];
rz(-2.634868) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.862344) q[0];
sx q[0];
rz(-1.5515258) q[0];
sx q[0];
rz(-1.4665682) q[0];
rz(2.7756491) q[2];
sx q[2];
rz(-2.2323713) q[2];
sx q[2];
rz(-2.5876074) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22257993) q[1];
sx q[1];
rz(-1.6363385) q[1];
sx q[1];
rz(-0.99187851) q[1];
rz(-0.30015035) q[3];
sx q[3];
rz(-0.86816278) q[3];
sx q[3];
rz(2.5182078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3158675) q[2];
sx q[2];
rz(-2.5434912) q[2];
sx q[2];
rz(-0.21710795) q[2];
rz(-1.0202967) q[3];
sx q[3];
rz(-1.812457) q[3];
sx q[3];
rz(-2.1609658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26647907) q[0];
sx q[0];
rz(-1.7784235) q[0];
sx q[0];
rz(-2.3734221) q[0];
rz(2.8504596) q[1];
sx q[1];
rz(-1.7765287) q[1];
sx q[1];
rz(-1.1870144) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465837) q[0];
sx q[0];
rz(-1.4872941) q[0];
sx q[0];
rz(-2.9299987) q[0];
x q[1];
rz(-2.111258) q[2];
sx q[2];
rz(-0.085496453) q[2];
sx q[2];
rz(0.70921113) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1817081) q[1];
sx q[1];
rz(-1.925029) q[1];
sx q[1];
rz(-1.6003057) q[1];
rz(-0.94908917) q[3];
sx q[3];
rz(-1.177985) q[3];
sx q[3];
rz(1.2248685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27176303) q[2];
sx q[2];
rz(-0.84438476) q[2];
sx q[2];
rz(-2.1738906) q[2];
rz(-0.23396954) q[3];
sx q[3];
rz(-2.440019) q[3];
sx q[3];
rz(0.35681891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02803239) q[0];
sx q[0];
rz(-1.5910633) q[0];
sx q[0];
rz(-2.0618942) q[0];
rz(-0.30403852) q[1];
sx q[1];
rz(-1.1540776) q[1];
sx q[1];
rz(-1.8255723) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2000421) q[0];
sx q[0];
rz(-2.594875) q[0];
sx q[0];
rz(0.23970516) q[0];
x q[1];
rz(-2.166376) q[2];
sx q[2];
rz(-1.6395431) q[2];
sx q[2];
rz(-2.749328) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58319891) q[1];
sx q[1];
rz(-0.89413645) q[1];
sx q[1];
rz(2.3469902) q[1];
x q[2];
rz(-0.06240978) q[3];
sx q[3];
rz(-1.9963328) q[3];
sx q[3];
rz(1.1264914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6198081) q[2];
sx q[2];
rz(-1.9079756) q[2];
sx q[2];
rz(-3.0002777) q[2];
rz(-0.58051336) q[3];
sx q[3];
rz(-1.4760735) q[3];
sx q[3];
rz(-2.9132402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2572131) q[0];
sx q[0];
rz(-2.3452106) q[0];
sx q[0];
rz(1.6656531) q[0];
rz(-2.2085704) q[1];
sx q[1];
rz(-0.7111744) q[1];
sx q[1];
rz(1.3915541) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2403657) q[0];
sx q[0];
rz(-2.3963442) q[0];
sx q[0];
rz(1.3979493) q[0];
x q[1];
rz(1.7250762) q[2];
sx q[2];
rz(-1.9898459) q[2];
sx q[2];
rz(1.2487457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7793624) q[1];
sx q[1];
rz(-1.0500776) q[1];
sx q[1];
rz(-0.034265072) q[1];
rz(-pi) q[2];
rz(2.0972033) q[3];
sx q[3];
rz(-0.59848753) q[3];
sx q[3];
rz(-0.55803669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.89620227) q[2];
sx q[2];
rz(-2.7707477) q[2];
sx q[2];
rz(0.24946269) q[2];
rz(1.6438515) q[3];
sx q[3];
rz(-1.4062873) q[3];
sx q[3];
rz(2.7264285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45667085) q[0];
sx q[0];
rz(-2.6455854) q[0];
sx q[0];
rz(1.7556835) q[0];
rz(-1.2617525) q[1];
sx q[1];
rz(-1.933814) q[1];
sx q[1];
rz(0.52880803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.945707) q[0];
sx q[0];
rz(-1.5778827) q[0];
sx q[0];
rz(1.12557) q[0];
x q[1];
rz(0.63839913) q[2];
sx q[2];
rz(-2.0631472) q[2];
sx q[2];
rz(-1.3251208) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3439548) q[1];
sx q[1];
rz(-2.4993163) q[1];
sx q[1];
rz(1.4587059) q[1];
x q[2];
rz(0.15915079) q[3];
sx q[3];
rz(-2.4503631) q[3];
sx q[3];
rz(-0.42014141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.95278198) q[2];
sx q[2];
rz(-0.47241259) q[2];
sx q[2];
rz(1.7702276) q[2];
rz(2.93907) q[3];
sx q[3];
rz(-1.6731508) q[3];
sx q[3];
rz(-1.9523581) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8093569) q[0];
sx q[0];
rz(-1.7928596) q[0];
sx q[0];
rz(-0.15383823) q[0];
rz(0.73506749) q[1];
sx q[1];
rz(-1.091833) q[1];
sx q[1];
rz(0.14911266) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47220889) q[0];
sx q[0];
rz(-2.4824004) q[0];
sx q[0];
rz(1.7071595) q[0];
rz(-pi) q[1];
rz(-1.2066226) q[2];
sx q[2];
rz(-1.7330164) q[2];
sx q[2];
rz(-0.71031865) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0267549) q[1];
sx q[1];
rz(-0.85875612) q[1];
sx q[1];
rz(0.65178443) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.074769) q[3];
sx q[3];
rz(-1.4912581) q[3];
sx q[3];
rz(-0.37802896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4601124) q[2];
sx q[2];
rz(-1.6296547) q[2];
sx q[2];
rz(-0.45664772) q[2];
rz(1.418142) q[3];
sx q[3];
rz(-0.8816312) q[3];
sx q[3];
rz(2.452623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5099455) q[0];
sx q[0];
rz(-2.8896285) q[0];
sx q[0];
rz(3.0889567) q[0];
rz(1.8317728) q[1];
sx q[1];
rz(-0.24886623) q[1];
sx q[1];
rz(1.9356669) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9034957) q[0];
sx q[0];
rz(-2.4795929) q[0];
sx q[0];
rz(-2.2308626) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56059273) q[2];
sx q[2];
rz(-1.7686604) q[2];
sx q[2];
rz(2.3677111) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.29621848) q[1];
sx q[1];
rz(-2.4476909) q[1];
sx q[1];
rz(-2.8945663) q[1];
x q[2];
rz(0.94843978) q[3];
sx q[3];
rz(-1.3788169) q[3];
sx q[3];
rz(2.6485659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84411088) q[2];
sx q[2];
rz(-1.7606807) q[2];
sx q[2];
rz(-2.5816176) q[2];
rz(-2.0848134) q[3];
sx q[3];
rz(-2.4323075) q[3];
sx q[3];
rz(-2.3287676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63447222) q[0];
sx q[0];
rz(-1.7422603) q[0];
sx q[0];
rz(1.5647474) q[0];
rz(-1.5693846) q[1];
sx q[1];
rz(-1.1610616) q[1];
sx q[1];
rz(-2.3326468) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3472373) q[0];
sx q[0];
rz(-2.6449727) q[0];
sx q[0];
rz(3.0625694) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3303512) q[2];
sx q[2];
rz(-1.277207) q[2];
sx q[2];
rz(-0.21101235) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8802277) q[1];
sx q[1];
rz(-0.61101171) q[1];
sx q[1];
rz(-1.5935807) q[1];
rz(-pi) q[2];
rz(-1.21502) q[3];
sx q[3];
rz(-0.18863931) q[3];
sx q[3];
rz(1.961018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.481015) q[2];
sx q[2];
rz(-0.95260859) q[2];
sx q[2];
rz(0.048952254) q[2];
rz(-1.281721) q[3];
sx q[3];
rz(-0.47544026) q[3];
sx q[3];
rz(1.6767282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82185164) q[0];
sx q[0];
rz(-0.26092437) q[0];
sx q[0];
rz(-1.9191746) q[0];
rz(-1.6659196) q[1];
sx q[1];
rz(-1.2011352) q[1];
sx q[1];
rz(2.6840721) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37402651) q[0];
sx q[0];
rz(-1.4806642) q[0];
sx q[0];
rz(-0.73390168) q[0];
rz(-pi) q[1];
rz(-0.40167146) q[2];
sx q[2];
rz(-2.0435512) q[2];
sx q[2];
rz(0.061701802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2603915) q[1];
sx q[1];
rz(-2.7335751) q[1];
sx q[1];
rz(0.53149207) q[1];
rz(-2.4549834) q[3];
sx q[3];
rz(-1.2401738) q[3];
sx q[3];
rz(2.7329684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3931302) q[2];
sx q[2];
rz(-2.8905383) q[2];
sx q[2];
rz(1.040323) q[2];
rz(-2.6535502) q[3];
sx q[3];
rz(-1.7010242) q[3];
sx q[3];
rz(-0.53781167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77448612) q[0];
sx q[0];
rz(-1.0716866) q[0];
sx q[0];
rz(-2.5197784) q[0];
rz(1.8472916) q[1];
sx q[1];
rz(-0.58302561) q[1];
sx q[1];
rz(2.2517712) q[1];
rz(1.1414083) q[2];
sx q[2];
rz(-2.5038237) q[2];
sx q[2];
rz(0.43546168) q[2];
rz(0.46066416) q[3];
sx q[3];
rz(-0.74429269) q[3];
sx q[3];
rz(-0.85891354) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
