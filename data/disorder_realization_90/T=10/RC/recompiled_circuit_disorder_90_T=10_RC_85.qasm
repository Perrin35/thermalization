OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(-0.087021526) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4092769) q[0];
sx q[0];
rz(-1.9415932) q[0];
sx q[0];
rz(2.2229574) q[0];
rz(0.1331698) q[2];
sx q[2];
rz(-0.77169092) q[2];
sx q[2];
rz(2.0193677) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7056071) q[1];
sx q[1];
rz(-1.4325303) q[1];
sx q[1];
rz(-2.0558753) q[1];
x q[2];
rz(-0.41214715) q[3];
sx q[3];
rz(-1.5282359) q[3];
sx q[3];
rz(1.3679078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13880754) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(-2.7677317) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(2.9132304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8841298) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(1.9467547) q[0];
rz(-0.082611235) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(3.1412178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21391348) q[0];
sx q[0];
rz(-1.1535026) q[0];
sx q[0];
rz(2.7620865) q[0];
x q[1];
rz(2.2194901) q[2];
sx q[2];
rz(-2.8173335) q[2];
sx q[2];
rz(-0.97933724) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.61987585) q[1];
sx q[1];
rz(-2.4881425) q[1];
sx q[1];
rz(0.39342777) q[1];
rz(1.625929) q[3];
sx q[3];
rz(-1.1109567) q[3];
sx q[3];
rz(-2.9210747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3327545) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(-2.9648798) q[2];
rz(-0.79408944) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90213838) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(-2.7368271) q[0];
rz(-1.8602712) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.8331029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3151911) q[0];
sx q[0];
rz(-2.0668525) q[0];
sx q[0];
rz(1.0415003) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8070418) q[2];
sx q[2];
rz(-2.0958732) q[2];
sx q[2];
rz(0.46307785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1723459) q[1];
sx q[1];
rz(-1.8241276) q[1];
sx q[1];
rz(-2.2971056) q[1];
rz(0.19034068) q[3];
sx q[3];
rz(-2.6958709) q[3];
sx q[3];
rz(2.1325071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9222766) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(-3.1211839) q[2];
rz(1.071788) q[3];
sx q[3];
rz(-1.1276378) q[3];
sx q[3];
rz(2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9999009) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(0.91598696) q[0];
rz(2.6782716) q[1];
sx q[1];
rz(-1.0659734) q[1];
sx q[1];
rz(-1.0571009) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0650478) q[0];
sx q[0];
rz(-1.1335982) q[0];
sx q[0];
rz(-1.5426427) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3154503) q[2];
sx q[2];
rz(-0.75781265) q[2];
sx q[2];
rz(-3.0405424) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6009439) q[1];
sx q[1];
rz(-1.800866) q[1];
sx q[1];
rz(-0.01407108) q[1];
x q[2];
rz(1.7980868) q[3];
sx q[3];
rz(-0.82633457) q[3];
sx q[3];
rz(-2.8518761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.89381924) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(2.7988953) q[2];
rz(1.6977067) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(-2.4263583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7856359) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(0.086439565) q[0];
rz(-1.3899639) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(-2.6729029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8353032) q[0];
sx q[0];
rz(-1.9625184) q[0];
sx q[0];
rz(2.3098371) q[0];
rz(-1.7519978) q[2];
sx q[2];
rz(-1.2672408) q[2];
sx q[2];
rz(-0.62523491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7680109) q[1];
sx q[1];
rz(-1.4410767) q[1];
sx q[1];
rz(-2.9380161) q[1];
rz(-pi) q[2];
rz(3.0607871) q[3];
sx q[3];
rz(-2.3256133) q[3];
sx q[3];
rz(-2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9995352) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(0.86432499) q[2];
rz(0.21720973) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(-0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42751673) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(0.19700225) q[0];
rz(-1.3621832) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(2.8053455) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4649481) q[0];
sx q[0];
rz(-1.1332382) q[0];
sx q[0];
rz(0.69533555) q[0];
rz(-pi) q[1];
rz(-2.5480812) q[2];
sx q[2];
rz(-2.0715908) q[2];
sx q[2];
rz(1.3154495) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5583401) q[1];
sx q[1];
rz(-1.5234158) q[1];
sx q[1];
rz(-0.58981311) q[1];
rz(0.69643173) q[3];
sx q[3];
rz(-1.4381572) q[3];
sx q[3];
rz(1.9385251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(2.2952648) q[2];
rz(1.2396631) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3939312) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(0.25892648) q[0];
rz(1.3461643) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-2.0475725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0712229) q[0];
sx q[0];
rz(-2.5734512) q[0];
sx q[0];
rz(2.2413261) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1292186) q[2];
sx q[2];
rz(-1.7933328) q[2];
sx q[2];
rz(2.6518133) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5309696) q[1];
sx q[1];
rz(-1.5429284) q[1];
sx q[1];
rz(1.4762957) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.592698) q[3];
sx q[3];
rz(-1.6688445) q[3];
sx q[3];
rz(2.736562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1743494) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(-2.8209177) q[2];
rz(0.43618068) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(-0.44803739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2748579) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(1.0048237) q[0];
rz(0.59016219) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(0.80387962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76669508) q[0];
sx q[0];
rz(-0.98537966) q[0];
sx q[0];
rz(-1.5234408) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7558394) q[2];
sx q[2];
rz(-2.1171326) q[2];
sx q[2];
rz(-1.9667369) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1357321) q[1];
sx q[1];
rz(-1.5605968) q[1];
sx q[1];
rz(-1.6594396) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6611093) q[3];
sx q[3];
rz(-0.90925018) q[3];
sx q[3];
rz(2.6722398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33264318) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-2.1311029) q[2];
rz(-0.60339749) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(-0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.5378961) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(0.4075152) q[0];
rz(-0.28911668) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(-2.3908652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.663561) q[0];
sx q[0];
rz(-0.41304092) q[0];
sx q[0];
rz(2.3703299) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27955555) q[2];
sx q[2];
rz(-1.4923555) q[2];
sx q[2];
rz(2.560175) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1322861) q[1];
sx q[1];
rz(-0.3711776) q[1];
sx q[1];
rz(0.98307857) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0575033) q[3];
sx q[3];
rz(-0.30232271) q[3];
sx q[3];
rz(1.7526527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0728545) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(1.1661952) q[2];
rz(-1.3646305) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5584548) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(-1.7512084) q[0];
rz(2.8109) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(1.6814544) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18171039) q[0];
sx q[0];
rz(-2.070817) q[0];
sx q[0];
rz(2.659003) q[0];
rz(-pi) q[1];
rz(-1.6129458) q[2];
sx q[2];
rz(-0.50439207) q[2];
sx q[2];
rz(-1.1500037) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9863661) q[1];
sx q[1];
rz(-0.44253293) q[1];
sx q[1];
rz(-0.44154422) q[1];
rz(-pi) q[2];
rz(-0.71698935) q[3];
sx q[3];
rz(-2.2319712) q[3];
sx q[3];
rz(-2.9644074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34974393) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(1.6798518) q[2];
rz(-2.0215624) q[3];
sx q[3];
rz(-2.5199065) q[3];
sx q[3];
rz(2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.5158952) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(-1.760578) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(2.3537221) q[2];
sx q[2];
rz(-1.315016) q[2];
sx q[2];
rz(-1.6223326) q[2];
rz(2.667726) q[3];
sx q[3];
rz(-0.73054536) q[3];
sx q[3];
rz(1.2252145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
