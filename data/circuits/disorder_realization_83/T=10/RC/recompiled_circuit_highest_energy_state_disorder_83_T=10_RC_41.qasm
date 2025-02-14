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
rz(-0.42521617) q[0];
sx q[0];
rz(4.4758237) q[0];
sx q[0];
rz(9.0444179) q[0];
rz(-1.3588139) q[1];
sx q[1];
rz(-2.9549197) q[1];
sx q[1];
rz(2.2169854) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43792576) q[0];
sx q[0];
rz(-0.80090085) q[0];
sx q[0];
rz(-0.4842224) q[0];
rz(-pi) q[1];
rz(-2.086928) q[2];
sx q[2];
rz(-1.7164111) q[2];
sx q[2];
rz(2.9235385) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8421665) q[1];
sx q[1];
rz(-1.6101082) q[1];
sx q[1];
rz(2.2996603) q[1];
rz(3.0713586) q[3];
sx q[3];
rz(-3.0133504) q[3];
sx q[3];
rz(-1.2508891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0164464) q[2];
sx q[2];
rz(-0.95637286) q[2];
sx q[2];
rz(-2.1535786) q[2];
rz(-0.2068578) q[3];
sx q[3];
rz(-1.7346953) q[3];
sx q[3];
rz(0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.67966953) q[0];
sx q[0];
rz(-1.672687) q[0];
sx q[0];
rz(2.8894506) q[0];
rz(3.120046) q[1];
sx q[1];
rz(-0.69308678) q[1];
sx q[1];
rz(-2.2861939) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48777521) q[0];
sx q[0];
rz(-0.027055351) q[0];
sx q[0];
rz(-1.6271126) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5163745) q[2];
sx q[2];
rz(-2.4897235) q[2];
sx q[2];
rz(-0.21833459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.508042) q[1];
sx q[1];
rz(-1.7910069) q[1];
sx q[1];
rz(-1.0255662) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68287373) q[3];
sx q[3];
rz(-1.3078017) q[3];
sx q[3];
rz(0.20786404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97027332) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(-0.26101905) q[2];
rz(3.1018992) q[3];
sx q[3];
rz(-1.7429765) q[3];
sx q[3];
rz(-0.54350129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4259341) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(-2.4470827) q[0];
rz(2.014324) q[1];
sx q[1];
rz(-2.290461) q[1];
sx q[1];
rz(0.99004254) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3304402) q[0];
sx q[0];
rz(-1.2227204) q[0];
sx q[0];
rz(1.483782) q[0];
x q[1];
rz(0.19845138) q[2];
sx q[2];
rz(-0.99764148) q[2];
sx q[2];
rz(1.0400553) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5824008) q[1];
sx q[1];
rz(-1.2123931) q[1];
sx q[1];
rz(-0.47689516) q[1];
x q[2];
rz(-0.53931358) q[3];
sx q[3];
rz(-2.1997872) q[3];
sx q[3];
rz(-1.1719538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8433044) q[2];
sx q[2];
rz(-2.1683606) q[2];
sx q[2];
rz(2.6514371) q[2];
rz(-2.7847024) q[3];
sx q[3];
rz(-2.7481952) q[3];
sx q[3];
rz(1.2662158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056034293) q[0];
sx q[0];
rz(-1.1857251) q[0];
sx q[0];
rz(-2.2991142) q[0];
rz(0.8017686) q[1];
sx q[1];
rz(-2.8623878) q[1];
sx q[1];
rz(0.78559771) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2548837) q[0];
sx q[0];
rz(-0.87866106) q[0];
sx q[0];
rz(-1.9002302) q[0];
rz(-pi) q[1];
rz(-2.9688492) q[2];
sx q[2];
rz(-2.4170791) q[2];
sx q[2];
rz(0.24562626) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.468049) q[1];
sx q[1];
rz(-2.1077413) q[1];
sx q[1];
rz(-1.931744) q[1];
rz(-1.7895643) q[3];
sx q[3];
rz(-0.76319088) q[3];
sx q[3];
rz(-1.4783875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0164612) q[2];
sx q[2];
rz(-2.3931914) q[2];
sx q[2];
rz(1.008519) q[2];
rz(-2.290001) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(-2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1763024) q[0];
sx q[0];
rz(-2.1602614) q[0];
sx q[0];
rz(1.0705795) q[0];
rz(-0.77752441) q[1];
sx q[1];
rz(-0.91064015) q[1];
sx q[1];
rz(1.0106962) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1309388) q[0];
sx q[0];
rz(-1.21216) q[0];
sx q[0];
rz(-0.50827311) q[0];
rz(-pi) q[1];
rz(-1.7456876) q[2];
sx q[2];
rz(-1.6351467) q[2];
sx q[2];
rz(-1.627587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3123828) q[1];
sx q[1];
rz(-2.0264033) q[1];
sx q[1];
rz(-0.73789755) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56212496) q[3];
sx q[3];
rz(-1.6549806) q[3];
sx q[3];
rz(0.69417324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8396478) q[2];
sx q[2];
rz(-0.10493111) q[2];
sx q[2];
rz(1.5555752) q[2];
rz(2.1258866) q[3];
sx q[3];
rz(-2.000122) q[3];
sx q[3];
rz(1.4122081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.79065901) q[0];
sx q[0];
rz(-2.9582294) q[0];
sx q[0];
rz(0.80107981) q[0];
rz(0.13310295) q[1];
sx q[1];
rz(-1.8201273) q[1];
sx q[1];
rz(2.7395693) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4426081) q[0];
sx q[0];
rz(-2.3056718) q[0];
sx q[0];
rz(1.1135654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5779675) q[2];
sx q[2];
rz(-0.29080393) q[2];
sx q[2];
rz(-0.70722843) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.0084486246) q[1];
sx q[1];
rz(-2.5767972) q[1];
sx q[1];
rz(1.0393591) q[1];
rz(-pi) q[2];
rz(2.9803863) q[3];
sx q[3];
rz(-2.4771002) q[3];
sx q[3];
rz(-2.4624069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63783995) q[2];
sx q[2];
rz(-1.736234) q[2];
sx q[2];
rz(-2.3657738) q[2];
rz(-1.4555629) q[3];
sx q[3];
rz(-1.173923) q[3];
sx q[3];
rz(2.1337401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044947226) q[0];
sx q[0];
rz(-2.2438887) q[0];
sx q[0];
rz(-2.4626379) q[0];
rz(1.2848162) q[1];
sx q[1];
rz(-2.4285451) q[1];
sx q[1];
rz(1.7624034) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0984707) q[0];
sx q[0];
rz(-0.80821191) q[0];
sx q[0];
rz(2.9007382) q[0];
rz(-pi) q[1];
rz(-1.2871615) q[2];
sx q[2];
rz(-0.29046187) q[2];
sx q[2];
rz(-2.1112372) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0022565) q[1];
sx q[1];
rz(-2.4508307) q[1];
sx q[1];
rz(0.73378566) q[1];
x q[2];
rz(1.3016765) q[3];
sx q[3];
rz(-2.9342071) q[3];
sx q[3];
rz(3.0548422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8046367) q[2];
sx q[2];
rz(-0.75504428) q[2];
sx q[2];
rz(-1.9026559) q[2];
rz(1.3321446) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(2.8968887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181344) q[0];
sx q[0];
rz(-1.1661538) q[0];
sx q[0];
rz(0.13846692) q[0];
rz(2.9578517) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(2.8750681) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375232) q[0];
sx q[0];
rz(-2.7752987) q[0];
sx q[0];
rz(-2.0384602) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1600003) q[2];
sx q[2];
rz(-0.60705429) q[2];
sx q[2];
rz(2.0701054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3764972) q[1];
sx q[1];
rz(-1.9899211) q[1];
sx q[1];
rz(0.061582743) q[1];
rz(0.0098465482) q[3];
sx q[3];
rz(-2.4146663) q[3];
sx q[3];
rz(1.7444514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0515685) q[2];
sx q[2];
rz(-1.8847621) q[2];
sx q[2];
rz(-2.9070692) q[2];
rz(-1.4759493) q[3];
sx q[3];
rz(-1.539295) q[3];
sx q[3];
rz(-2.3852111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.804857) q[0];
sx q[0];
rz(-2.2754301) q[0];
sx q[0];
rz(-2.3418703) q[0];
rz(-2.5526478) q[1];
sx q[1];
rz(-2.2942693) q[1];
sx q[1];
rz(0.2074997) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070415592) q[0];
sx q[0];
rz(-1.1780329) q[0];
sx q[0];
rz(1.7123187) q[0];
rz(-pi) q[1];
rz(-1.9249179) q[2];
sx q[2];
rz(-2.3975044) q[2];
sx q[2];
rz(0.085372774) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5821857) q[1];
sx q[1];
rz(-1.1731802) q[1];
sx q[1];
rz(-0.84487265) q[1];
x q[2];
rz(0.028771632) q[3];
sx q[3];
rz(-1.1011916) q[3];
sx q[3];
rz(2.4503277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3346682) q[2];
sx q[2];
rz(-0.88403264) q[2];
sx q[2];
rz(-2.7743288) q[2];
rz(1.4372829) q[3];
sx q[3];
rz(-1.1742914) q[3];
sx q[3];
rz(1.9678763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87840286) q[0];
sx q[0];
rz(-1.0028361) q[0];
sx q[0];
rz(2.3314085) q[0];
rz(-2.0267678) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(2.5883163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1076521) q[0];
sx q[0];
rz(-2.2080511) q[0];
sx q[0];
rz(3.0725543) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79375879) q[2];
sx q[2];
rz(-2.0362034) q[2];
sx q[2];
rz(2.2528354) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1571036) q[1];
sx q[1];
rz(-1.8525288) q[1];
sx q[1];
rz(2.9473122) q[1];
rz(1.9235885) q[3];
sx q[3];
rz(-0.81656352) q[3];
sx q[3];
rz(-1.0366131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0193923) q[2];
sx q[2];
rz(-2.4330752) q[2];
sx q[2];
rz(-0.23966399) q[2];
rz(1.8137118) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(-0.46752587) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0620621) q[0];
sx q[0];
rz(-1.7876328) q[0];
sx q[0];
rz(-2.7476516) q[0];
rz(2.486034) q[1];
sx q[1];
rz(-1.9656904) q[1];
sx q[1];
rz(0.94737731) q[1];
rz(1.9026359) q[2];
sx q[2];
rz(-1.9658026) q[2];
sx q[2];
rz(2.6323242) q[2];
rz(1.3169133) q[3];
sx q[3];
rz(-1.0051654) q[3];
sx q[3];
rz(-1.2011423) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
