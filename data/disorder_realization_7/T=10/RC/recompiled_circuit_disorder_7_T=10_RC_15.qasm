OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(-0.20110826) q[0];
rz(-1.3970628) q[1];
sx q[1];
rz(-1.8048598) q[1];
sx q[1];
rz(-0.62682682) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39591046) q[0];
sx q[0];
rz(-2.9183309) q[0];
sx q[0];
rz(-1.0049055) q[0];
x q[1];
rz(0.5958545) q[2];
sx q[2];
rz(-2.7239954) q[2];
sx q[2];
rz(0.28996224) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1498897) q[1];
sx q[1];
rz(-2.457379) q[1];
sx q[1];
rz(-1.0690772) q[1];
rz(-pi) q[2];
rz(-1.2960035) q[3];
sx q[3];
rz(-1.7405675) q[3];
sx q[3];
rz(2.3613289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(-1.8480776) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-2.2938426) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936546) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(-0.98518103) q[0];
rz(0.39918104) q[1];
sx q[1];
rz(-1.8789411) q[1];
sx q[1];
rz(-3.0342297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3283008) q[0];
sx q[0];
rz(-1.0257226) q[0];
sx q[0];
rz(0.48304793) q[0];
x q[1];
rz(2.7254843) q[2];
sx q[2];
rz(-2.5118718) q[2];
sx q[2];
rz(2.1925418) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0641891) q[1];
sx q[1];
rz(-1.940155) q[1];
sx q[1];
rz(-2.2662152) q[1];
rz(-pi) q[2];
rz(1.4025877) q[3];
sx q[3];
rz(-1.9544) q[3];
sx q[3];
rz(2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(0.04034986) q[0];
rz(0.57998747) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(1.0823762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4989935) q[0];
sx q[0];
rz(-1.2046308) q[0];
sx q[0];
rz(-2.5901592) q[0];
rz(1.9627377) q[2];
sx q[2];
rz(-0.38092962) q[2];
sx q[2];
rz(2.4286963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.3324932) q[1];
sx q[1];
rz(0.2281245) q[1];
x q[2];
rz(2.9186967) q[3];
sx q[3];
rz(-1.2548903) q[3];
sx q[3];
rz(-0.99881682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0064156) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(-0.55389261) q[2];
rz(-1.6484377) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(2.5890787) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(-0.035382263) q[0];
rz(-2.1454051) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(-2.6534973) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29823829) q[0];
sx q[0];
rz(-0.47180628) q[0];
sx q[0];
rz(1.7422373) q[0];
rz(-pi) q[1];
rz(-0.98056356) q[2];
sx q[2];
rz(-1.0733986) q[2];
sx q[2];
rz(0.24888466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0779873) q[1];
sx q[1];
rz(-1.6997153) q[1];
sx q[1];
rz(1.414826) q[1];
rz(1.0933236) q[3];
sx q[3];
rz(-2.100088) q[3];
sx q[3];
rz(-2.570591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(-2.6341237) q[2];
rz(-1.9594225) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(1.2549843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59947157) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(-0.3702634) q[0];
rz(0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(-0.19651861) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.927658) q[0];
sx q[0];
rz(-1.2688583) q[0];
sx q[0];
rz(1.2769804) q[0];
rz(-pi) q[1];
rz(1.1097124) q[2];
sx q[2];
rz(-1.5421252) q[2];
sx q[2];
rz(-1.920514) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12381183) q[1];
sx q[1];
rz(-1.31243) q[1];
sx q[1];
rz(1.7523132) q[1];
x q[2];
rz(2.8653141) q[3];
sx q[3];
rz(-1.4071138) q[3];
sx q[3];
rz(2.7460263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(-2.2244942) q[2];
rz(-2.3069978) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.31345263) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(2.939558) q[0];
rz(2.0687885) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.7477759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48702792) q[0];
sx q[0];
rz(-1.2396887) q[0];
sx q[0];
rz(-1.6513255) q[0];
x q[1];
rz(-0.11153531) q[2];
sx q[2];
rz(-1.1984014) q[2];
sx q[2];
rz(2.6967449) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8762293) q[1];
sx q[1];
rz(-2.6637042) q[1];
sx q[1];
rz(-2.1761314) q[1];
rz(-pi) q[2];
rz(1.1145669) q[3];
sx q[3];
rz(-2.2531415) q[3];
sx q[3];
rz(-1.8519459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(-0.8423155) q[2];
rz(2.0541644) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6120646) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(2.1202309) q[0];
rz(1.1313324) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-0.21025118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049858) q[0];
sx q[0];
rz(-0.14963089) q[0];
sx q[0];
rz(-1.11087) q[0];
rz(0.37961752) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(-1.1773674) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.10340313) q[1];
sx q[1];
rz(-1.2545171) q[1];
sx q[1];
rz(1.3693387) q[1];
rz(-pi) q[2];
rz(-0.85150163) q[3];
sx q[3];
rz(-2.0264894) q[3];
sx q[3];
rz(-1.9853026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(-2.3794877) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(0.41346082) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.64637) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(-0.58724171) q[0];
rz(-0.44627055) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(0.97672021) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8551089) q[0];
sx q[0];
rz(-0.47874641) q[0];
sx q[0];
rz(-2.8141862) q[0];
rz(-3.061053) q[2];
sx q[2];
rz(-1.409568) q[2];
sx q[2];
rz(-1.5494407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50423065) q[1];
sx q[1];
rz(-2.0081875) q[1];
sx q[1];
rz(-0.38423844) q[1];
rz(-pi) q[2];
rz(2.5891853) q[3];
sx q[3];
rz(-2.2033785) q[3];
sx q[3];
rz(-1.6747024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7405159) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-0.55244279) q[2];
rz(-0.46594122) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66822806) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(0.92064944) q[0];
rz(0.051041516) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(-2.2424973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2313401) q[0];
sx q[0];
rz(-1.1465342) q[0];
sx q[0];
rz(0.92696855) q[0];
rz(-pi) q[1];
rz(-2.1295206) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(0.65288359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2194849) q[1];
sx q[1];
rz(-1.9003632) q[1];
sx q[1];
rz(0.12164128) q[1];
rz(-pi) q[2];
rz(-0.00021342834) q[3];
sx q[3];
rz(-1.2006239) q[3];
sx q[3];
rz(2.7406524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(2.7446279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2755651) q[0];
sx q[0];
rz(-2.4336928) q[0];
sx q[0];
rz(-0.78535725) q[0];
x q[1];
rz(0.88231477) q[2];
sx q[2];
rz(-2.6803662) q[2];
sx q[2];
rz(3.0416833) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46170235) q[1];
sx q[1];
rz(-1.9063408) q[1];
sx q[1];
rz(3.0479792) q[1];
rz(2.7361717) q[3];
sx q[3];
rz(-2.0339031) q[3];
sx q[3];
rz(-2.6904358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94840702) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(0.096654264) q[2];
rz(1.7276673) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5861355) q[0];
sx q[0];
rz(-1.1145034) q[0];
sx q[0];
rz(0.282423) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(-1.6554228) q[2];
sx q[2];
rz(-1.463269) q[2];
sx q[2];
rz(-0.29381897) q[2];
rz(0.11937033) q[3];
sx q[3];
rz(-1.3138249) q[3];
sx q[3];
rz(-1.8520595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];