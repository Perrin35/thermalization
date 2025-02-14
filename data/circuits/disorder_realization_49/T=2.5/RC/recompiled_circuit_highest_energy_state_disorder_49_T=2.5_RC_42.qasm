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
rz(-3.1168923) q[0];
sx q[0];
rz(-1.8827117) q[0];
sx q[0];
rz(-1.2727241) q[0];
rz(-1.7495573) q[1];
sx q[1];
rz(-1.7311544) q[1];
sx q[1];
rz(-2.2478204) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0860505) q[0];
sx q[0];
rz(-2.1889733) q[0];
sx q[0];
rz(-2.2937141) q[0];
rz(-pi) q[1];
rz(-0.70510257) q[2];
sx q[2];
rz(-2.7573708) q[2];
sx q[2];
rz(-2.9085858) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8402579) q[1];
sx q[1];
rz(-2.4221759) q[1];
sx q[1];
rz(0.69566984) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5886542) q[3];
sx q[3];
rz(-1.2027338) q[3];
sx q[3];
rz(-2.5758655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36665234) q[2];
sx q[2];
rz(-1.6795936) q[2];
sx q[2];
rz(-0.53202638) q[2];
rz(2.4074647) q[3];
sx q[3];
rz(-0.2747772) q[3];
sx q[3];
rz(-1.1067357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4183913) q[0];
sx q[0];
rz(-2.1558303) q[0];
sx q[0];
rz(3.0611839) q[0];
rz(-2.6037727) q[1];
sx q[1];
rz(-1.0686914) q[1];
sx q[1];
rz(2.417876) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2576669) q[0];
sx q[0];
rz(-0.91086713) q[0];
sx q[0];
rz(1.0038478) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2369464) q[2];
sx q[2];
rz(-1.1084741) q[2];
sx q[2];
rz(1.8682299) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5608175) q[1];
sx q[1];
rz(-2.5530949) q[1];
sx q[1];
rz(-1.568041) q[1];
x q[2];
rz(-0.81775093) q[3];
sx q[3];
rz(-1.9649385) q[3];
sx q[3];
rz(0.5241636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2878652) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(0.90616027) q[2];
rz(0.3624889) q[3];
sx q[3];
rz(-1.5972842) q[3];
sx q[3];
rz(0.046886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.872181) q[0];
sx q[0];
rz(-1.325664) q[0];
sx q[0];
rz(-2.0500702) q[0];
rz(3.0190234) q[1];
sx q[1];
rz(-1.4361607) q[1];
sx q[1];
rz(-1.7950119) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99341644) q[0];
sx q[0];
rz(-1.5412586) q[0];
sx q[0];
rz(1.2146773) q[0];
rz(-1.8982688) q[2];
sx q[2];
rz(-0.93743126) q[2];
sx q[2];
rz(-2.4302389) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57438499) q[1];
sx q[1];
rz(-1.775601) q[1];
sx q[1];
rz(-3.0887207) q[1];
rz(-pi) q[2];
rz(1.4580549) q[3];
sx q[3];
rz(-1.0442956) q[3];
sx q[3];
rz(2.2248588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.71318212) q[2];
sx q[2];
rz(-1.1298263) q[2];
sx q[2];
rz(-2.909519) q[2];
rz(-1.9706767) q[3];
sx q[3];
rz(-0.62058312) q[3];
sx q[3];
rz(-0.28461972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21122268) q[0];
sx q[0];
rz(-2.3486597) q[0];
sx q[0];
rz(1.5027745) q[0];
rz(0.59208313) q[1];
sx q[1];
rz(-1.1555669) q[1];
sx q[1];
rz(-0.88964644) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29763313) q[0];
sx q[0];
rz(-2.0951162) q[0];
sx q[0];
rz(-0.53589937) q[0];
rz(-pi) q[1];
rz(-2.0060894) q[2];
sx q[2];
rz(-2.6084628) q[2];
sx q[2];
rz(1.0077623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1210416) q[1];
sx q[1];
rz(-2.7747323) q[1];
sx q[1];
rz(1.0479142) q[1];
rz(-pi) q[2];
rz(-0.82075228) q[3];
sx q[3];
rz(-2.6057557) q[3];
sx q[3];
rz(1.395063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4952116) q[2];
sx q[2];
rz(-1.4633598) q[2];
sx q[2];
rz(-2.4656673) q[2];
rz(1.4365139) q[3];
sx q[3];
rz(-2.8016475) q[3];
sx q[3];
rz(-0.16726141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9306358) q[0];
sx q[0];
rz(-0.3322596) q[0];
sx q[0];
rz(2.9462872) q[0];
rz(-0.12061128) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(-0.19048555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57225655) q[0];
sx q[0];
rz(-1.4829206) q[0];
sx q[0];
rz(1.5059307) q[0];
x q[1];
rz(-2.4649925) q[2];
sx q[2];
rz(-1.6381421) q[2];
sx q[2];
rz(-1.8563893) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.021564158) q[1];
sx q[1];
rz(-3.0282674) q[1];
sx q[1];
rz(-2.6647407) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56511527) q[3];
sx q[3];
rz(-0.69639528) q[3];
sx q[3];
rz(0.90000641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5087937) q[2];
sx q[2];
rz(-1.4845279) q[2];
sx q[2];
rz(-1.18139) q[2];
rz(1.7975636) q[3];
sx q[3];
rz(-0.7414147) q[3];
sx q[3];
rz(0.77176362) q[3];
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
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8089499) q[0];
sx q[0];
rz(-1.7985666) q[0];
sx q[0];
rz(-1.1395662) q[0];
rz(2.471916) q[1];
sx q[1];
rz(-0.91595903) q[1];
sx q[1];
rz(2.0119827) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2655339) q[0];
sx q[0];
rz(-1.8647412) q[0];
sx q[0];
rz(2.8493899) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8390719) q[2];
sx q[2];
rz(-0.74069689) q[2];
sx q[2];
rz(-0.64756913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1150951) q[1];
sx q[1];
rz(-1.9352176) q[1];
sx q[1];
rz(-1.3624191) q[1];
rz(-pi) q[2];
rz(-2.4242006) q[3];
sx q[3];
rz(-2.903707) q[3];
sx q[3];
rz(0.57008509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7149522) q[2];
sx q[2];
rz(-1.7317829) q[2];
sx q[2];
rz(2.7362774) q[2];
rz(2.3146368) q[3];
sx q[3];
rz(-2.5155641) q[3];
sx q[3];
rz(2.2820182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3655132) q[0];
sx q[0];
rz(-3.1190393) q[0];
sx q[0];
rz(-2.205701) q[0];
rz(2.2731958) q[1];
sx q[1];
rz(-0.52803841) q[1];
sx q[1];
rz(1.4194999) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1552733) q[0];
sx q[0];
rz(-2.5012448) q[0];
sx q[0];
rz(0.61926432) q[0];
rz(-2.9090857) q[2];
sx q[2];
rz(-1.075313) q[2];
sx q[2];
rz(1.2459618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3277927) q[1];
sx q[1];
rz(-2.2918502) q[1];
sx q[1];
rz(0.93443296) q[1];
x q[2];
rz(-1.5704182) q[3];
sx q[3];
rz(-1.5725563) q[3];
sx q[3];
rz(-1.5906217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.821227) q[2];
sx q[2];
rz(-1.7419523) q[2];
sx q[2];
rz(-1.0199176) q[2];
rz(-2.6323281) q[3];
sx q[3];
rz(-1.441322) q[3];
sx q[3];
rz(0.078484623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7996063) q[0];
sx q[0];
rz(-1.6260363) q[0];
sx q[0];
rz(-1.9827783) q[0];
rz(0.47053567) q[1];
sx q[1];
rz(-1.7262986) q[1];
sx q[1];
rz(1.6758957) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3780791) q[0];
sx q[0];
rz(-1.5944885) q[0];
sx q[0];
rz(-0.050934249) q[0];
rz(-1.2242555) q[2];
sx q[2];
rz(-2.2316791) q[2];
sx q[2];
rz(-0.61049547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.012124) q[1];
sx q[1];
rz(-2.2290342) q[1];
sx q[1];
rz(1.4736498) q[1];
rz(-0.64792525) q[3];
sx q[3];
rz(-1.8019391) q[3];
sx q[3];
rz(2.1794139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9965685) q[2];
sx q[2];
rz(-1.8137167) q[2];
sx q[2];
rz(2.3222951) q[2];
rz(-0.79646349) q[3];
sx q[3];
rz(-0.4070681) q[3];
sx q[3];
rz(1.4888633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2040937) q[0];
sx q[0];
rz(-0.096385328) q[0];
sx q[0];
rz(-0.82164422) q[0];
rz(-0.67283982) q[1];
sx q[1];
rz(-2.5655589) q[1];
sx q[1];
rz(1.0427262) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0053528) q[0];
sx q[0];
rz(-1.3658821) q[0];
sx q[0];
rz(1.4968275) q[0];
rz(-pi) q[1];
rz(2.2006532) q[2];
sx q[2];
rz(-2.2374638) q[2];
sx q[2];
rz(2.710944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.38775165) q[1];
sx q[1];
rz(-0.70742224) q[1];
sx q[1];
rz(-2.161987) q[1];
rz(-pi) q[2];
rz(-0.53570536) q[3];
sx q[3];
rz(-2.8495478) q[3];
sx q[3];
rz(-3.0391708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78160703) q[2];
sx q[2];
rz(-1.8693417) q[2];
sx q[2];
rz(-2.2992415) q[2];
rz(-2.602747) q[3];
sx q[3];
rz(-1.4875965) q[3];
sx q[3];
rz(-2.6174788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6738324) q[0];
sx q[0];
rz(-0.14586511) q[0];
sx q[0];
rz(1.2354596) q[0];
rz(-2.1927059) q[1];
sx q[1];
rz(-0.83273879) q[1];
sx q[1];
rz(-1.4804776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87029261) q[0];
sx q[0];
rz(-2.2402431) q[0];
sx q[0];
rz(-0.35737329) q[0];
rz(-pi) q[1];
rz(-0.5625426) q[2];
sx q[2];
rz(-1.0357366) q[2];
sx q[2];
rz(-1.8059274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9292641) q[1];
sx q[1];
rz(-0.18106279) q[1];
sx q[1];
rz(-0.4275112) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8126902) q[3];
sx q[3];
rz(-1.8902361) q[3];
sx q[3];
rz(-1.1238568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.71198717) q[2];
sx q[2];
rz(-1.4236071) q[2];
sx q[2];
rz(-0.31022662) q[2];
rz(0.45983908) q[3];
sx q[3];
rz(-2.583677) q[3];
sx q[3];
rz(-1.3742113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1772169) q[0];
sx q[0];
rz(-1.1342659) q[0];
sx q[0];
rz(-2.076617) q[0];
rz(-2.2433157) q[1];
sx q[1];
rz(-2.5986462) q[1];
sx q[1];
rz(2.6959261) q[1];
rz(1.0679792) q[2];
sx q[2];
rz(-0.72952727) q[2];
sx q[2];
rz(-2.9222957) q[2];
rz(-2.2326917) q[3];
sx q[3];
rz(-1.881897) q[3];
sx q[3];
rz(1.9492634) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
