OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0527394) q[0];
sx q[0];
rz(-2.7272447) q[0];
sx q[0];
rz(-2.156884) q[0];
rz(-0.86496487) q[1];
sx q[1];
rz(-0.84671658) q[1];
sx q[1];
rz(-2.4196978) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979209) q[0];
sx q[0];
rz(-2.0204442) q[0];
sx q[0];
rz(-1.1051154) q[0];
rz(-pi) q[1];
rz(-2.7593132) q[2];
sx q[2];
rz(-0.77300249) q[2];
sx q[2];
rz(-2.892008) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.940168) q[1];
sx q[1];
rz(-0.076751953) q[1];
sx q[1];
rz(1.9606664) q[1];
rz(-pi) q[2];
rz(-0.15318449) q[3];
sx q[3];
rz(-1.9537874) q[3];
sx q[3];
rz(2.9944978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7559173) q[2];
sx q[2];
rz(-0.52272457) q[2];
sx q[2];
rz(0.44504607) q[2];
rz(1.2612777) q[3];
sx q[3];
rz(-1.6877561) q[3];
sx q[3];
rz(-1.6311579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234914) q[0];
sx q[0];
rz(-2.0913251) q[0];
sx q[0];
rz(2.4871248) q[0];
rz(-2.0596313) q[1];
sx q[1];
rz(-1.8257717) q[1];
sx q[1];
rz(-1.6771603) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0102756) q[0];
sx q[0];
rz(-0.47476381) q[0];
sx q[0];
rz(1.6807589) q[0];
x q[1];
rz(-3.0686321) q[2];
sx q[2];
rz(-0.96809298) q[2];
sx q[2];
rz(-2.631244) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1480746) q[1];
sx q[1];
rz(-1.52964) q[1];
sx q[1];
rz(2.4012474) q[1];
rz(-pi) q[2];
rz(0.62276472) q[3];
sx q[3];
rz(-0.55985427) q[3];
sx q[3];
rz(2.6948158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4945041) q[2];
sx q[2];
rz(-3.0454128) q[2];
sx q[2];
rz(0.75867009) q[2];
rz(1.8132973) q[3];
sx q[3];
rz(-1.0814861) q[3];
sx q[3];
rz(-1.7648511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7741622) q[0];
sx q[0];
rz(-1.6663015) q[0];
sx q[0];
rz(-3.1040763) q[0];
rz(1.0388177) q[1];
sx q[1];
rz(-0.33375868) q[1];
sx q[1];
rz(-2.2822101) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20221449) q[0];
sx q[0];
rz(-2.4945745) q[0];
sx q[0];
rz(0.31536343) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0065303) q[2];
sx q[2];
rz(-0.76530582) q[2];
sx q[2];
rz(-2.0489592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22166477) q[1];
sx q[1];
rz(-1.8916191) q[1];
sx q[1];
rz(1.1236091) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.061305) q[3];
sx q[3];
rz(-2.1119364) q[3];
sx q[3];
rz(3.0681572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3694156) q[2];
sx q[2];
rz(-0.93537664) q[2];
sx q[2];
rz(-2.3773362) q[2];
rz(2.1956826) q[3];
sx q[3];
rz(-1.9352244) q[3];
sx q[3];
rz(-2.2435772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94364828) q[0];
sx q[0];
rz(-2.2479489) q[0];
sx q[0];
rz(-1.3125032) q[0];
rz(-0.30458826) q[1];
sx q[1];
rz(-2.1423788) q[1];
sx q[1];
rz(-0.33590683) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0806341) q[0];
sx q[0];
rz(-1.4984382) q[0];
sx q[0];
rz(-3.1273187) q[0];
rz(-pi) q[1];
rz(0.96520378) q[2];
sx q[2];
rz(-1.9450257) q[2];
sx q[2];
rz(0.40520129) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0842753) q[1];
sx q[1];
rz(-2.0641987) q[1];
sx q[1];
rz(-2.0396903) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5139989) q[3];
sx q[3];
rz(-1.3496163) q[3];
sx q[3];
rz(-3.055453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37561068) q[2];
sx q[2];
rz(-0.42372647) q[2];
sx q[2];
rz(1.83164) q[2];
rz(-1.2951819) q[3];
sx q[3];
rz(-2.3056307) q[3];
sx q[3];
rz(-1.0999934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64422166) q[0];
sx q[0];
rz(-0.27232429) q[0];
sx q[0];
rz(-0.8154794) q[0];
rz(-2.4244335) q[1];
sx q[1];
rz(-1.6970789) q[1];
sx q[1];
rz(-1.1517634) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718147) q[0];
sx q[0];
rz(-1.0564532) q[0];
sx q[0];
rz(1.8215979) q[0];
rz(1.8134591) q[2];
sx q[2];
rz(-1.2753701) q[2];
sx q[2];
rz(-2.5385419) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5260701) q[1];
sx q[1];
rz(-2.263219) q[1];
sx q[1];
rz(-0.39151616) q[1];
x q[2];
rz(1.5342496) q[3];
sx q[3];
rz(-1.8724955) q[3];
sx q[3];
rz(0.73153618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.069245) q[2];
sx q[2];
rz(-0.77584156) q[2];
sx q[2];
rz(-1.5838712) q[2];
rz(0.97258687) q[3];
sx q[3];
rz(-1.2310622) q[3];
sx q[3];
rz(-1.7044273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69777456) q[0];
sx q[0];
rz(-2.8688718) q[0];
sx q[0];
rz(-0.66993237) q[0];
rz(2.2386235) q[1];
sx q[1];
rz(-0.65330708) q[1];
sx q[1];
rz(0.088937581) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085468) q[0];
sx q[0];
rz(-2.4886697) q[0];
sx q[0];
rz(1.9050951) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8208262) q[2];
sx q[2];
rz(-0.93968117) q[2];
sx q[2];
rz(-2.4970412) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6429421) q[1];
sx q[1];
rz(-1.6410259) q[1];
sx q[1];
rz(-1.8204017) q[1];
x q[2];
rz(-2.0518725) q[3];
sx q[3];
rz(-1.4566453) q[3];
sx q[3];
rz(2.6237918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6637806) q[2];
sx q[2];
rz(-2.0974443) q[2];
sx q[2];
rz(0.021154724) q[2];
rz(-0.92709213) q[3];
sx q[3];
rz(-1.5320211) q[3];
sx q[3];
rz(-2.482614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053719036) q[0];
sx q[0];
rz(-1.3608195) q[0];
sx q[0];
rz(2.0045643) q[0];
rz(-2.673705) q[1];
sx q[1];
rz(-1.7158022) q[1];
sx q[1];
rz(-0.55380026) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1374986) q[0];
sx q[0];
rz(-1.2826254) q[0];
sx q[0];
rz(2.1779446) q[0];
x q[1];
rz(-1.1365898) q[2];
sx q[2];
rz(-1.8946664) q[2];
sx q[2];
rz(-0.5262652) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3674161) q[1];
sx q[1];
rz(-1.8715845) q[1];
sx q[1];
rz(-2.1117626) q[1];
x q[2];
rz(1.475369) q[3];
sx q[3];
rz(-0.57668873) q[3];
sx q[3];
rz(-0.63729034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0849453) q[2];
sx q[2];
rz(-2.2979996) q[2];
sx q[2];
rz(-0.6960558) q[2];
rz(2.7737235) q[3];
sx q[3];
rz(-2.0335679) q[3];
sx q[3];
rz(2.8554816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87357658) q[0];
sx q[0];
rz(-1.3883256) q[0];
sx q[0];
rz(0.8152813) q[0];
rz(-2.6878327) q[1];
sx q[1];
rz(-2.2397857) q[1];
sx q[1];
rz(2.544983) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0946569) q[0];
sx q[0];
rz(-0.8917745) q[0];
sx q[0];
rz(-0.89767098) q[0];
rz(-pi) q[1];
rz(-1.6758133) q[2];
sx q[2];
rz(-2.3987282) q[2];
sx q[2];
rz(-0.32203963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4441159) q[1];
sx q[1];
rz(-1.8702862) q[1];
sx q[1];
rz(1.0018574) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1659307) q[3];
sx q[3];
rz(-1.941276) q[3];
sx q[3];
rz(-2.3080829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0348908) q[2];
sx q[2];
rz(-1.4406349) q[2];
sx q[2];
rz(1.2197257) q[2];
rz(-1.2446416) q[3];
sx q[3];
rz(-1.605426) q[3];
sx q[3];
rz(0.19967782) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.225746) q[0];
sx q[0];
rz(-0.93996489) q[0];
sx q[0];
rz(-1.655727) q[0];
rz(-2.0303717) q[1];
sx q[1];
rz(-1.2948371) q[1];
sx q[1];
rz(-1.3908386) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6538453) q[0];
sx q[0];
rz(-1.5731205) q[0];
sx q[0];
rz(0.33727686) q[0];
rz(-pi) q[1];
rz(0.31220372) q[2];
sx q[2];
rz(-2.7810568) q[2];
sx q[2];
rz(-2.3106239) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0842485) q[1];
sx q[1];
rz(-2.930713) q[1];
sx q[1];
rz(-1.4046304) q[1];
x q[2];
rz(1.9377356) q[3];
sx q[3];
rz(-1.8979567) q[3];
sx q[3];
rz(-0.1042455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86604649) q[2];
sx q[2];
rz(-1.6547497) q[2];
sx q[2];
rz(-0.1869959) q[2];
rz(-0.20728076) q[3];
sx q[3];
rz(-2.5076301) q[3];
sx q[3];
rz(2.0980339) q[3];
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
rz(-2.2279376) q[0];
sx q[0];
rz(-0.74988237) q[0];
sx q[0];
rz(3.1392745) q[0];
rz(1.2222611) q[1];
sx q[1];
rz(-1.5536676) q[1];
sx q[1];
rz(1.4171756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42444705) q[0];
sx q[0];
rz(-1.7758649) q[0];
sx q[0];
rz(-1.7030832) q[0];
rz(-pi) q[1];
rz(0.50062407) q[2];
sx q[2];
rz(-0.60030314) q[2];
sx q[2];
rz(0.17498091) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8479653) q[1];
sx q[1];
rz(-2.9305998) q[1];
sx q[1];
rz(-1.3354098) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3145261) q[3];
sx q[3];
rz(-0.60171222) q[3];
sx q[3];
rz(0.0046241143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2554539) q[2];
sx q[2];
rz(-1.6243434) q[2];
sx q[2];
rz(-2.9296866) q[2];
rz(0.40031561) q[3];
sx q[3];
rz(-0.90462697) q[3];
sx q[3];
rz(1.3306085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0935852) q[0];
sx q[0];
rz(-1.891991) q[0];
sx q[0];
rz(1.4461507) q[0];
rz(0.49225898) q[1];
sx q[1];
rz(-1.6620363) q[1];
sx q[1];
rz(-1.5580039) q[1];
rz(-1.7854431) q[2];
sx q[2];
rz(-1.3673906) q[2];
sx q[2];
rz(-2.2833952) q[2];
rz(0.51385469) q[3];
sx q[3];
rz(-2.0217635) q[3];
sx q[3];
rz(1.3213984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
