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
rz(0.9917292) q[0];
sx q[0];
rz(-0.60757929) q[0];
sx q[0];
rz(-2.1460331) q[0];
rz(0.52752703) q[1];
sx q[1];
rz(-2.1105284) q[1];
sx q[1];
rz(0.48479015) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3178892) q[0];
sx q[0];
rz(-1.4394338) q[0];
sx q[0];
rz(2.0844368) q[0];
x q[1];
rz(0.77930696) q[2];
sx q[2];
rz(-1.8911084) q[2];
sx q[2];
rz(0.28398289) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97580662) q[1];
sx q[1];
rz(-0.96250223) q[1];
sx q[1];
rz(-3.0700141) q[1];
x q[2];
rz(0.23956128) q[3];
sx q[3];
rz(-0.29920855) q[3];
sx q[3];
rz(-0.20209683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63554865) q[2];
sx q[2];
rz(-2.1940993) q[2];
sx q[2];
rz(0.49723899) q[2];
rz(2.5908568) q[3];
sx q[3];
rz(-2.4516055) q[3];
sx q[3];
rz(-2.0042888) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9841442) q[0];
sx q[0];
rz(-1.7406311) q[0];
sx q[0];
rz(0.76341158) q[0];
rz(0.58452559) q[1];
sx q[1];
rz(-1.3817363) q[1];
sx q[1];
rz(-2.1302628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69788591) q[0];
sx q[0];
rz(-2.9666565) q[0];
sx q[0];
rz(-0.849443) q[0];
x q[1];
rz(2.9052034) q[2];
sx q[2];
rz(-1.4589785) q[2];
sx q[2];
rz(-0.37416247) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3015131) q[1];
sx q[1];
rz(-0.40427819) q[1];
sx q[1];
rz(1.0227281) q[1];
rz(-2.743558) q[3];
sx q[3];
rz(-1.0938489) q[3];
sx q[3];
rz(1.2881035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6611019) q[2];
sx q[2];
rz(-1.9384408) q[2];
sx q[2];
rz(1.0002381) q[2];
rz(0.58146042) q[3];
sx q[3];
rz(-0.73271078) q[3];
sx q[3];
rz(-1.2015517) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30705273) q[0];
sx q[0];
rz(-2.5426799) q[0];
sx q[0];
rz(-0.55759984) q[0];
rz(-0.3071951) q[1];
sx q[1];
rz(-2.5418044) q[1];
sx q[1];
rz(-0.44816005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9014943) q[0];
sx q[0];
rz(-1.0567825) q[0];
sx q[0];
rz(-1.4185403) q[0];
x q[1];
rz(2.9286341) q[2];
sx q[2];
rz(-2.0816275) q[2];
sx q[2];
rz(1.2868483) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8413755) q[1];
sx q[1];
rz(-0.63134495) q[1];
sx q[1];
rz(-1.5247099) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43436111) q[3];
sx q[3];
rz(-1.5944527) q[3];
sx q[3];
rz(1.0756186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8223411) q[2];
sx q[2];
rz(-2.4288869) q[2];
sx q[2];
rz(2.5437497) q[2];
rz(0.44761014) q[3];
sx q[3];
rz(-1.2667344) q[3];
sx q[3];
rz(-3.083057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41281259) q[0];
sx q[0];
rz(-0.030981177) q[0];
sx q[0];
rz(-2.431751) q[0];
rz(-2.1624883) q[1];
sx q[1];
rz(-0.85936463) q[1];
sx q[1];
rz(1.8202579) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6431491) q[0];
sx q[0];
rz(-1.6263026) q[0];
sx q[0];
rz(-0.70141354) q[0];
rz(-pi) q[1];
rz(-2.9843281) q[2];
sx q[2];
rz(-0.7866106) q[2];
sx q[2];
rz(0.32341126) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85985124) q[1];
sx q[1];
rz(-0.70003372) q[1];
sx q[1];
rz(1.1477877) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0783441) q[3];
sx q[3];
rz(-2.5406917) q[3];
sx q[3];
rz(2.7067824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9909624) q[2];
sx q[2];
rz(-1.4272828) q[2];
sx q[2];
rz(-2.4347351) q[2];
rz(-1.3458716) q[3];
sx q[3];
rz(-1.462505) q[3];
sx q[3];
rz(0.70485392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3020346) q[0];
sx q[0];
rz(-0.78080451) q[0];
sx q[0];
rz(-1.6395521) q[0];
rz(-1.543965) q[1];
sx q[1];
rz(-0.85999703) q[1];
sx q[1];
rz(1.4942716) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2236007) q[0];
sx q[0];
rz(-0.57044464) q[0];
sx q[0];
rz(1.5738652) q[0];
rz(-pi) q[1];
rz(1.5635033) q[2];
sx q[2];
rz(-2.1661758) q[2];
sx q[2];
rz(-3.1356168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0493023) q[1];
sx q[1];
rz(-0.43152025) q[1];
sx q[1];
rz(-2.3445361) q[1];
rz(3.1190514) q[3];
sx q[3];
rz(-2.1858741) q[3];
sx q[3];
rz(-0.32451567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8813701) q[2];
sx q[2];
rz(-1.2726731) q[2];
sx q[2];
rz(-0.55848813) q[2];
rz(2.6327366) q[3];
sx q[3];
rz(-1.5285834) q[3];
sx q[3];
rz(1.2671027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59638554) q[0];
sx q[0];
rz(-1.243243) q[0];
sx q[0];
rz(0.16045706) q[0];
rz(-2.4682553) q[1];
sx q[1];
rz(-1.2168177) q[1];
sx q[1];
rz(2.014726) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1864439) q[0];
sx q[0];
rz(-1.8732244) q[0];
sx q[0];
rz(-1.9496098) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8931863) q[2];
sx q[2];
rz(-0.83288542) q[2];
sx q[2];
rz(2.5423342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2963449) q[1];
sx q[1];
rz(-0.61734491) q[1];
sx q[1];
rz(1.2165353) q[1];
x q[2];
rz(2.1127719) q[3];
sx q[3];
rz(-1.3888956) q[3];
sx q[3];
rz(-1.0396837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0141853) q[2];
sx q[2];
rz(-1.5421474) q[2];
sx q[2];
rz(0.7737774) q[2];
rz(0.989178) q[3];
sx q[3];
rz(-0.4072322) q[3];
sx q[3];
rz(2.0686843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9147375) q[0];
sx q[0];
rz(-1.4829153) q[0];
sx q[0];
rz(-2.3684655) q[0];
rz(-1.585656) q[1];
sx q[1];
rz(-1.8525886) q[1];
sx q[1];
rz(-1.9645436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0601412) q[0];
sx q[0];
rz(-2.3075275) q[0];
sx q[0];
rz(3.0990776) q[0];
rz(-pi) q[1];
rz(-2.8297205) q[2];
sx q[2];
rz(-2.1215028) q[2];
sx q[2];
rz(-2.4518509) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6625632) q[1];
sx q[1];
rz(-1.8445043) q[1];
sx q[1];
rz(0.91027963) q[1];
rz(-pi) q[2];
rz(-0.66399948) q[3];
sx q[3];
rz(-2.917508) q[3];
sx q[3];
rz(-0.9730207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.78544632) q[2];
sx q[2];
rz(-0.36474228) q[2];
sx q[2];
rz(0.44343239) q[2];
rz(2.7555079) q[3];
sx q[3];
rz(-1.4177136) q[3];
sx q[3];
rz(2.940787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76483738) q[0];
sx q[0];
rz(-1.5489464) q[0];
sx q[0];
rz(0.03373294) q[0];
rz(2.4434166) q[1];
sx q[1];
rz(-0.59711421) q[1];
sx q[1];
rz(0.66506344) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18298223) q[0];
sx q[0];
rz(-2.0183051) q[0];
sx q[0];
rz(-0.29194269) q[0];
rz(-pi) q[1];
rz(-2.6680634) q[2];
sx q[2];
rz(-0.90190926) q[2];
sx q[2];
rz(-2.7471126) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4951952) q[1];
sx q[1];
rz(-2.2943999) q[1];
sx q[1];
rz(-0.88314914) q[1];
rz(-1.5499209) q[3];
sx q[3];
rz(-1.2513147) q[3];
sx q[3];
rz(0.11939458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72208059) q[2];
sx q[2];
rz(-0.91700143) q[2];
sx q[2];
rz(-1.3631932) q[2];
rz(0.52115399) q[3];
sx q[3];
rz(-1.3241973) q[3];
sx q[3];
rz(0.49191973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9687965) q[0];
sx q[0];
rz(-1.3461312) q[0];
sx q[0];
rz(2.3222493) q[0];
rz(-0.68483812) q[1];
sx q[1];
rz(-1.8543154) q[1];
sx q[1];
rz(-3.0192764) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4222833) q[0];
sx q[0];
rz(-2.3580812) q[0];
sx q[0];
rz(2.1859403) q[0];
rz(-1.0367461) q[2];
sx q[2];
rz(-2.4754562) q[2];
sx q[2];
rz(0.18830794) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6284991) q[1];
sx q[1];
rz(-2.4446396) q[1];
sx q[1];
rz(-0.51817399) q[1];
rz(-pi) q[2];
rz(-2.0360394) q[3];
sx q[3];
rz(-1.0330345) q[3];
sx q[3];
rz(-0.28013438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43716064) q[2];
sx q[2];
rz(-1.8081343) q[2];
sx q[2];
rz(2.2740347) q[2];
rz(-1.7433172) q[3];
sx q[3];
rz(-2.7531392) q[3];
sx q[3];
rz(-2.5625663) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.687998) q[0];
sx q[0];
rz(-1.2247676) q[0];
sx q[0];
rz(2.0523742) q[0];
rz(-3.094063) q[1];
sx q[1];
rz(-1.0462953) q[1];
sx q[1];
rz(0.68738031) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6733765) q[0];
sx q[0];
rz(-0.49096732) q[0];
sx q[0];
rz(-1.6223905) q[0];
x q[1];
rz(2.0014761) q[2];
sx q[2];
rz(-1.8269369) q[2];
sx q[2];
rz(1.3821951) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8721622) q[1];
sx q[1];
rz(-2.3331642) q[1];
sx q[1];
rz(0.6489268) q[1];
rz(3.0484588) q[3];
sx q[3];
rz(-0.93619213) q[3];
sx q[3];
rz(-2.840691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.99531949) q[2];
sx q[2];
rz(-0.61890382) q[2];
sx q[2];
rz(-1.4825561) q[2];
rz(0.65685529) q[3];
sx q[3];
rz(-1.1494136) q[3];
sx q[3];
rz(-2.759554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.252608) q[0];
sx q[0];
rz(-1.1532619) q[0];
sx q[0];
rz(1.1070195) q[0];
rz(0.57869115) q[1];
sx q[1];
rz(-2.3042669) q[1];
sx q[1];
rz(-0.28490983) q[1];
rz(1.4877697) q[2];
sx q[2];
rz(-1.9258693) q[2];
sx q[2];
rz(-3.0167962) q[2];
rz(-2.6246351) q[3];
sx q[3];
rz(-1.6421972) q[3];
sx q[3];
rz(0.35785892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
