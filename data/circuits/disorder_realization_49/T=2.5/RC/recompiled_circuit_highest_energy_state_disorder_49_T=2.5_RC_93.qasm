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
rz(0.024700392) q[0];
sx q[0];
rz(-1.2588809) q[0];
sx q[0];
rz(1.2727241) q[0];
rz(1.3920353) q[1];
sx q[1];
rz(-1.4104383) q[1];
sx q[1];
rz(2.2478204) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0964684) q[0];
sx q[0];
rz(-0.91334263) q[0];
sx q[0];
rz(-0.7492926) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8270764) q[2];
sx q[2];
rz(-1.2813142) q[2];
sx q[2];
rz(2.6319401) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8302248) q[1];
sx q[1];
rz(-1.1347949) q[1];
sx q[1];
rz(2.5495851) q[1];
rz(-0.36811604) q[3];
sx q[3];
rz(-1.5541346) q[3];
sx q[3];
rz(0.9986432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.36665234) q[2];
sx q[2];
rz(-1.4619991) q[2];
sx q[2];
rz(-2.6095663) q[2];
rz(0.73412791) q[3];
sx q[3];
rz(-0.2747772) q[3];
sx q[3];
rz(-2.0348569) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4183913) q[0];
sx q[0];
rz(-2.1558303) q[0];
sx q[0];
rz(-3.0611839) q[0];
rz(0.53781992) q[1];
sx q[1];
rz(-2.0729013) q[1];
sx q[1];
rz(-2.417876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88392576) q[0];
sx q[0];
rz(-0.91086713) q[0];
sx q[0];
rz(1.0038478) q[0];
x q[1];
rz(0.56494464) q[2];
sx q[2];
rz(-2.1568642) q[2];
sx q[2];
rz(-0.63462574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1538633) q[1];
sx q[1];
rz(-1.5692668) q[1];
sx q[1];
rz(-2.1592924) q[1];
rz(-pi) q[2];
rz(1.0243914) q[3];
sx q[3];
rz(-0.83163578) q[3];
sx q[3];
rz(-1.4354768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2878652) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(-2.2354324) q[2];
rz(-2.7791038) q[3];
sx q[3];
rz(-1.5972842) q[3];
sx q[3];
rz(0.046886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2694117) q[0];
sx q[0];
rz(-1.325664) q[0];
sx q[0];
rz(2.0500702) q[0];
rz(0.12256924) q[1];
sx q[1];
rz(-1.7054319) q[1];
sx q[1];
rz(1.3465808) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5751981) q[0];
sx q[0];
rz(-1.2148396) q[0];
sx q[0];
rz(3.1100789) q[0];
x q[1];
rz(-2.7287219) q[2];
sx q[2];
rz(-0.70253583) q[2];
sx q[2];
rz(-1.9091878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1344188) q[1];
sx q[1];
rz(-1.5190304) q[1];
sx q[1];
rz(1.365713) q[1];
rz(1.6835378) q[3];
sx q[3];
rz(-1.0442956) q[3];
sx q[3];
rz(-2.2248588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4284105) q[2];
sx q[2];
rz(-2.0117663) q[2];
sx q[2];
rz(0.23207363) q[2];
rz(-1.170916) q[3];
sx q[3];
rz(-0.62058312) q[3];
sx q[3];
rz(0.28461972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.93037) q[0];
sx q[0];
rz(-0.79293293) q[0];
sx q[0];
rz(-1.5027745) q[0];
rz(-0.59208313) q[1];
sx q[1];
rz(-1.9860257) q[1];
sx q[1];
rz(-0.88964644) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98418035) q[0];
sx q[0];
rz(-2.0286467) q[0];
sx q[0];
rz(0.97869398) q[0];
rz(-pi) q[1];
rz(-2.0621262) q[2];
sx q[2];
rz(-1.7867807) q[2];
sx q[2];
rz(-0.94397533) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0977323) q[1];
sx q[1];
rz(-1.3907038) q[1];
sx q[1];
rz(1.8921683) q[1];
rz(-pi) q[2];
rz(-0.82075228) q[3];
sx q[3];
rz(-2.6057557) q[3];
sx q[3];
rz(1.395063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64638102) q[2];
sx q[2];
rz(-1.4633598) q[2];
sx q[2];
rz(-2.4656673) q[2];
rz(1.4365139) q[3];
sx q[3];
rz(-0.33994514) q[3];
sx q[3];
rz(0.16726141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.2109569) q[0];
sx q[0];
rz(-2.8093331) q[0];
sx q[0];
rz(-0.19530547) q[0];
rz(3.0209814) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(2.9511071) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99283907) q[0];
sx q[0];
rz(-1.6354113) q[0];
sx q[0];
rz(3.0535327) q[0];
x q[1];
rz(-1.6570857) q[2];
sx q[2];
rz(-0.89601529) q[2];
sx q[2];
rz(-2.9099885) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0234649) q[1];
sx q[1];
rz(-1.6227229) q[1];
sx q[1];
rz(3.0408188) q[1];
rz(-pi) q[2];
rz(-2.5764774) q[3];
sx q[3];
rz(-0.69639528) q[3];
sx q[3];
rz(-0.90000641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5087937) q[2];
sx q[2];
rz(-1.6570647) q[2];
sx q[2];
rz(-1.18139) q[2];
rz(1.3440291) q[3];
sx q[3];
rz(-2.400178) q[3];
sx q[3];
rz(0.77176362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8089499) q[0];
sx q[0];
rz(-1.3430261) q[0];
sx q[0];
rz(-2.0020265) q[0];
rz(-0.66967669) q[1];
sx q[1];
rz(-0.91595903) q[1];
sx q[1];
rz(-1.12961) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07189508) q[0];
sx q[0];
rz(-2.7301333) q[0];
sx q[0];
rz(0.81019391) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8390719) q[2];
sx q[2];
rz(-0.74069689) q[2];
sx q[2];
rz(0.64756913) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5796389) q[1];
sx q[1];
rz(-2.7241251) q[1];
sx q[1];
rz(-2.6446656) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9608742) q[3];
sx q[3];
rz(-1.4152495) q[3];
sx q[3];
rz(0.29744086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7149522) q[2];
sx q[2];
rz(-1.7317829) q[2];
sx q[2];
rz(0.40531522) q[2];
rz(0.82695588) q[3];
sx q[3];
rz(-2.5155641) q[3];
sx q[3];
rz(-2.2820182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3655132) q[0];
sx q[0];
rz(-0.022553355) q[0];
sx q[0];
rz(-2.205701) q[0];
rz(0.86839688) q[1];
sx q[1];
rz(-0.52803841) q[1];
sx q[1];
rz(1.7220928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4286242) q[0];
sx q[0];
rz(-1.0626864) q[0];
sx q[0];
rz(-1.9789805) q[0];
x q[1];
rz(-2.9090857) q[2];
sx q[2];
rz(-2.0662796) q[2];
sx q[2];
rz(-1.2459618) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3277927) q[1];
sx q[1];
rz(-2.2918502) q[1];
sx q[1];
rz(2.2071597) q[1];
rz(-pi) q[2];
rz(-0.0017599864) q[3];
sx q[3];
rz(-1.5711745) q[3];
sx q[3];
rz(-0.019824713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3203656) q[2];
sx q[2];
rz(-1.3996404) q[2];
sx q[2];
rz(1.0199176) q[2];
rz(0.50926456) q[3];
sx q[3];
rz(-1.7002707) q[3];
sx q[3];
rz(3.063108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7996063) q[0];
sx q[0];
rz(-1.5155563) q[0];
sx q[0];
rz(1.1588143) q[0];
rz(-0.47053567) q[1];
sx q[1];
rz(-1.7262986) q[1];
sx q[1];
rz(-1.6758957) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138595) q[0];
sx q[0];
rz(-3.0854221) q[0];
sx q[0];
rz(0.43561952) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7297425) q[2];
sx q[2];
rz(-0.73397103) q[2];
sx q[2];
rz(-3.0628772) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1294687) q[1];
sx q[1];
rz(-2.2290342) q[1];
sx q[1];
rz(1.6679428) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37181446) q[3];
sx q[3];
rz(-0.68228693) q[3];
sx q[3];
rz(0.90250795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9965685) q[2];
sx q[2];
rz(-1.327876) q[2];
sx q[2];
rz(0.81929755) q[2];
rz(-0.79646349) q[3];
sx q[3];
rz(-2.7345246) q[3];
sx q[3];
rz(-1.4888633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9374989) q[0];
sx q[0];
rz(-3.0452073) q[0];
sx q[0];
rz(-2.3199484) q[0];
rz(-0.67283982) q[1];
sx q[1];
rz(-2.5655589) q[1];
sx q[1];
rz(-2.0988665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4194787) q[0];
sx q[0];
rz(-1.4983777) q[0];
sx q[0];
rz(2.9361322) q[0];
rz(-pi) q[1];
rz(0.94093948) q[2];
sx q[2];
rz(-0.90412882) q[2];
sx q[2];
rz(-0.43064865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.753841) q[1];
sx q[1];
rz(-2.4341704) q[1];
sx q[1];
rz(2.161987) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8886077) q[3];
sx q[3];
rz(-1.4232985) q[3];
sx q[3];
rz(1.9852888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3599856) q[2];
sx q[2];
rz(-1.2722509) q[2];
sx q[2];
rz(2.2992415) q[2];
rz(2.602747) q[3];
sx q[3];
rz(-1.4875965) q[3];
sx q[3];
rz(-0.52411383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46776029) q[0];
sx q[0];
rz(-2.9957275) q[0];
sx q[0];
rz(1.2354596) q[0];
rz(0.94888672) q[1];
sx q[1];
rz(-0.83273879) q[1];
sx q[1];
rz(-1.4804776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2713) q[0];
sx q[0];
rz(-2.2402431) q[0];
sx q[0];
rz(2.7842194) q[0];
rz(2.3034978) q[2];
sx q[2];
rz(-0.75586719) q[2];
sx q[2];
rz(2.2264293) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9292641) q[1];
sx q[1];
rz(-0.18106279) q[1];
sx q[1];
rz(2.7140815) q[1];
rz(-0.62678316) q[3];
sx q[3];
rz(-2.7434182) q[3];
sx q[3];
rz(2.6836723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.71198717) q[2];
sx q[2];
rz(-1.4236071) q[2];
sx q[2];
rz(-2.831366) q[2];
rz(-2.6817536) q[3];
sx q[3];
rz(-2.583677) q[3];
sx q[3];
rz(-1.3742113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9643758) q[0];
sx q[0];
rz(-2.0073267) q[0];
sx q[0];
rz(1.0649756) q[0];
rz(-0.89827697) q[1];
sx q[1];
rz(-0.54294642) q[1];
sx q[1];
rz(-0.44566659) q[1];
rz(1.0679792) q[2];
sx q[2];
rz(-0.72952727) q[2];
sx q[2];
rz(-2.9222957) q[2];
rz(0.38705714) q[3];
sx q[3];
rz(-2.195812) q[3];
sx q[3];
rz(0.14433911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
