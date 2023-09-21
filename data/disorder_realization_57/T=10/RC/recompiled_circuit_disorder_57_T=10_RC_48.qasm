OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52743071) q[0];
sx q[0];
rz(3.951374) q[0];
sx q[0];
rz(9.9561719) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(-1.7778492) q[1];
sx q[1];
rz(-1.2030503) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42703585) q[0];
sx q[0];
rz(-1.3792975) q[0];
sx q[0];
rz(1.3546076) q[0];
rz(-pi) q[1];
rz(1.7877903) q[2];
sx q[2];
rz(-1.8197682) q[2];
sx q[2];
rz(0.72460246) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6984451) q[1];
sx q[1];
rz(-1.5138211) q[1];
sx q[1];
rz(-2.9064102) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12227998) q[3];
sx q[3];
rz(-2.8465726) q[3];
sx q[3];
rz(0.89000851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.27665859) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-0.19908389) q[2];
rz(-1.6137326) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(-0.73408192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1428225) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(-0.81940991) q[0];
rz(0.2858513) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(1.2664638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74840141) q[0];
sx q[0];
rz(-1.9494434) q[0];
sx q[0];
rz(-2.4875531) q[0];
rz(-1.9752713) q[2];
sx q[2];
rz(-1.3504488) q[2];
sx q[2];
rz(-1.5420367) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.61050615) q[1];
sx q[1];
rz(-0.55263457) q[1];
sx q[1];
rz(-0.83341877) q[1];
x q[2];
rz(2.500202) q[3];
sx q[3];
rz(-2.4603599) q[3];
sx q[3];
rz(0.89279803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9849898) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(-1.9796237) q[2];
rz(0.087163838) q[3];
sx q[3];
rz(-1.6379084) q[3];
sx q[3];
rz(3.0676837) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6067628) q[0];
sx q[0];
rz(-2.02089) q[0];
sx q[0];
rz(-0.30360046) q[0];
rz(-1.3820232) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(0.64750013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.203813) q[0];
sx q[0];
rz(-2.505629) q[0];
sx q[0];
rz(-3.0657835) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0258255) q[2];
sx q[2];
rz(-0.93169824) q[2];
sx q[2];
rz(-2.7884723) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4422569) q[1];
sx q[1];
rz(-1.2067716) q[1];
sx q[1];
rz(-2.060021) q[1];
x q[2];
rz(2.5472574) q[3];
sx q[3];
rz(-2.1777275) q[3];
sx q[3];
rz(-2.2567574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8422164) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(1.5807318) q[2];
rz(2.1598699) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9075539) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(2.2696944) q[0];
rz(2.7543228) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(0.29104582) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0204578) q[0];
sx q[0];
rz(-0.99636787) q[0];
sx q[0];
rz(2.8681884) q[0];
x q[1];
rz(1.8918858) q[2];
sx q[2];
rz(-1.2920657) q[2];
sx q[2];
rz(0.89821399) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17109045) q[1];
sx q[1];
rz(-0.88741747) q[1];
sx q[1];
rz(0.6716397) q[1];
rz(-pi) q[2];
rz(-1.8157186) q[3];
sx q[3];
rz(-2.5513253) q[3];
sx q[3];
rz(0.11412379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7724093) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(-2.5975361) q[2];
rz(2.6323075) q[3];
sx q[3];
rz(-2.0777168) q[3];
sx q[3];
rz(0.86597401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13957025) q[0];
sx q[0];
rz(-2.1152571) q[0];
sx q[0];
rz(-0.64506662) q[0];
rz(-2.5158665) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(-2.5114139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1030578) q[0];
sx q[0];
rz(-1.4487421) q[0];
sx q[0];
rz(-0.5794258) q[0];
rz(-1.1516045) q[2];
sx q[2];
rz(-2.345511) q[2];
sx q[2];
rz(-2.5615356) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.899257) q[1];
sx q[1];
rz(-1.3471654) q[1];
sx q[1];
rz(-1.0244589) q[1];
rz(0.24124055) q[3];
sx q[3];
rz(-2.5706228) q[3];
sx q[3];
rz(-0.2824479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75382918) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(-1.7774263) q[2];
rz(-1.2498614) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(-0.92145872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(1.0867778) q[0];
sx q[0];
rz(-0.58371109) q[0];
sx q[0];
rz(2.4247647) q[0];
rz(-2.7038799) q[1];
sx q[1];
rz(-0.68044674) q[1];
sx q[1];
rz(2.3278918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1610634) q[0];
sx q[0];
rz(-1.6335532) q[0];
sx q[0];
rz(1.784523) q[0];
x q[1];
rz(-2.5021624) q[2];
sx q[2];
rz(-0.8317906) q[2];
sx q[2];
rz(-2.8384428) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4390347) q[1];
sx q[1];
rz(-2.3137989) q[1];
sx q[1];
rz(1.1101767) q[1];
rz(-1.2171451) q[3];
sx q[3];
rz(-1.2815223) q[3];
sx q[3];
rz(1.7464964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2563236) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(1.9539333) q[2];
rz(2.2096283) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(0.37117547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1535783) q[0];
sx q[0];
rz(-2.1126641) q[0];
sx q[0];
rz(-0.6860835) q[0];
rz(-1.5230806) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(-2.4783321) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8535186) q[0];
sx q[0];
rz(-2.220904) q[0];
sx q[0];
rz(-3.1115565) q[0];
x q[1];
rz(-0.63753193) q[2];
sx q[2];
rz(-1.7469284) q[2];
sx q[2];
rz(-0.1917563) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.38695947) q[1];
sx q[1];
rz(-2.6135751) q[1];
sx q[1];
rz(1.9792884) q[1];
rz(-pi) q[2];
rz(2.456326) q[3];
sx q[3];
rz(-2.4589834) q[3];
sx q[3];
rz(2.9311789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66701165) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(0.015080301) q[2];
rz(-1.0117426) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(-2.570178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.052977234) q[0];
sx q[0];
rz(-0.73260728) q[0];
sx q[0];
rz(0.10738871) q[0];
rz(-2.836851) q[1];
sx q[1];
rz(-1.3184897) q[1];
sx q[1];
rz(-2.6838141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7055571) q[0];
sx q[0];
rz(-1.2495263) q[0];
sx q[0];
rz(0.9106439) q[0];
x q[1];
rz(2.5758178) q[2];
sx q[2];
rz(-2.3104295) q[2];
sx q[2];
rz(2.4772252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4630423) q[1];
sx q[1];
rz(-1.2150587) q[1];
sx q[1];
rz(2.8929391) q[1];
rz(-pi) q[2];
rz(-1.9321953) q[3];
sx q[3];
rz(-1.0482338) q[3];
sx q[3];
rz(-2.2613329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36499873) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.7101074) q[2];
rz(1.7163904) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(1.2861929) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15437056) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(-0.95348683) q[0];
rz(1.3257239) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(-0.58473933) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1138211) q[0];
sx q[0];
rz(-2.6391533) q[0];
sx q[0];
rz(0.2343982) q[0];
rz(-pi) q[1];
rz(-1.2051177) q[2];
sx q[2];
rz(-2.0846016) q[2];
sx q[2];
rz(-0.50525451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.19496275) q[1];
sx q[1];
rz(-1.1319036) q[1];
sx q[1];
rz(1.6473325) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9986266) q[3];
sx q[3];
rz(-1.3060095) q[3];
sx q[3];
rz(-0.20204443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8087625) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(-2.8653223) q[2];
rz(-0.55109465) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(2.7579257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0693414) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(-1.1570702) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(1.6190593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9690543) q[0];
sx q[0];
rz(-1.8803055) q[0];
sx q[0];
rz(0.68185101) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68799512) q[2];
sx q[2];
rz(-1.8383467) q[2];
sx q[2];
rz(-2.5431395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5787705) q[1];
sx q[1];
rz(-2.4914503) q[1];
sx q[1];
rz(2.8679718) q[1];
rz(-pi) q[2];
rz(-1.1679681) q[3];
sx q[3];
rz(-1.9625003) q[3];
sx q[3];
rz(-0.21564461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9860501) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(0.84214169) q[2];
rz(1.5367674) q[3];
sx q[3];
rz(-1.3804881) q[3];
sx q[3];
rz(3.1150637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90606541) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(2.4189667) q[1];
sx q[1];
rz(-2.2650748) q[1];
sx q[1];
rz(1.5320019) q[1];
rz(-1.4459544) q[2];
sx q[2];
rz(-2.9308133) q[2];
sx q[2];
rz(1.4598893) q[2];
rz(3.0620861) q[3];
sx q[3];
rz(-2.2600997) q[3];
sx q[3];
rz(-3.1117677) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
