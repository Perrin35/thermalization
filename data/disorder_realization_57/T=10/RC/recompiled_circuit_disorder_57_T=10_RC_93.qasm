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
rz(1.3637435) q[1];
sx q[1];
rz(10.627828) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7145568) q[0];
sx q[0];
rz(-1.3792975) q[0];
sx q[0];
rz(1.7869851) q[0];
rz(-0.7026303) q[2];
sx q[2];
rz(-2.812817) q[2];
sx q[2];
rz(1.4544912) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.10567084) q[1];
sx q[1];
rz(-0.24186132) q[1];
sx q[1];
rz(2.9015404) q[1];
rz(-pi) q[2];
rz(3.0193127) q[3];
sx q[3];
rz(-2.8465726) q[3];
sx q[3];
rz(2.2515841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8649341) q[2];
sx q[2];
rz(-2.0099137) q[2];
sx q[2];
rz(-2.9425088) q[2];
rz(-1.52786) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(-2.4075107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1428225) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(2.3221827) q[0];
rz(0.2858513) q[1];
sx q[1];
rz(-0.91393036) q[1];
sx q[1];
rz(-1.2664638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74840141) q[0];
sx q[0];
rz(-1.9494434) q[0];
sx q[0];
rz(-0.65403954) q[0];
x q[1];
rz(-2.9026047) q[2];
sx q[2];
rz(-1.9649437) q[2];
sx q[2];
rz(3.0195401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5310865) q[1];
sx q[1];
rz(-2.5889581) q[1];
sx q[1];
rz(2.3081739) q[1];
rz(-pi) q[2];
rz(2.0224138) q[3];
sx q[3];
rz(-2.0997117) q[3];
sx q[3];
rz(-3.0144453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1566029) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(-1.9796237) q[2];
rz(0.087163838) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(0.073908977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6067628) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(-0.30360046) q[0];
rz(1.3820232) q[1];
sx q[1];
rz(-1.8654114) q[1];
sx q[1];
rz(-2.4940925) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.203813) q[0];
sx q[0];
rz(-0.63596361) q[0];
sx q[0];
rz(3.0657835) q[0];
x q[1];
rz(-0.60909231) q[2];
sx q[2];
rz(-0.81431544) q[2];
sx q[2];
rz(-0.44037214) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6806879) q[1];
sx q[1];
rz(-0.60084963) q[1];
sx q[1];
rz(-2.2520573) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59433523) q[3];
sx q[3];
rz(-2.1777275) q[3];
sx q[3];
rz(-2.2567574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2993762) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(-1.5608609) q[2];
rz(-0.98172274) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9075539) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(2.2696944) q[0];
rz(0.38726989) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(-0.29104582) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5440118) q[0];
sx q[0];
rz(-0.62950069) q[0];
sx q[0];
rz(-1.9660216) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8918858) q[2];
sx q[2];
rz(-1.849527) q[2];
sx q[2];
rz(0.89821399) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9705022) q[1];
sx q[1];
rz(-0.88741747) q[1];
sx q[1];
rz(-0.6716397) q[1];
rz(-pi) q[2];
rz(1.8157186) q[3];
sx q[3];
rz(-0.59026736) q[3];
sx q[3];
rz(-3.0274689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7724093) q[2];
sx q[2];
rz(-0.79493752) q[2];
sx q[2];
rz(-0.54405653) q[2];
rz(-0.50928515) q[3];
sx q[3];
rz(-2.0777168) q[3];
sx q[3];
rz(-2.2756186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0020224) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(-0.64506662) q[0];
rz(-0.6257261) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(2.5114139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1030578) q[0];
sx q[0];
rz(-1.4487421) q[0];
sx q[0];
rz(-2.5621668) q[0];
rz(-pi) q[1];
rz(-0.39406392) q[2];
sx q[2];
rz(-0.85959083) q[2];
sx q[2];
rz(-1.1472536) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.899257) q[1];
sx q[1];
rz(-1.7944272) q[1];
sx q[1];
rz(1.0244589) q[1];
rz(-pi) q[2];
rz(0.24124055) q[3];
sx q[3];
rz(-0.57096982) q[3];
sx q[3];
rz(-2.8591448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3877635) q[2];
sx q[2];
rz(-1.8176273) q[2];
sx q[2];
rz(-1.3641664) q[2];
rz(1.2498614) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(0.92145872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0548148) q[0];
sx q[0];
rz(-0.58371109) q[0];
sx q[0];
rz(-2.4247647) q[0];
rz(2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(2.3278918) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87156224) q[0];
sx q[0];
rz(-0.22261482) q[0];
sx q[0];
rz(-1.2827669) q[0];
rz(0.6394303) q[2];
sx q[2];
rz(-2.3098021) q[2];
sx q[2];
rz(-0.30314988) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3354213) q[1];
sx q[1];
rz(-2.2911982) q[1];
sx q[1];
rz(-2.6909188) q[1];
rz(2.2807062) q[3];
sx q[3];
rz(-0.45300278) q[3];
sx q[3];
rz(-2.3080254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8852691) q[2];
sx q[2];
rz(-0.43748125) q[2];
sx q[2];
rz(1.9539333) q[2];
rz(2.2096283) q[3];
sx q[3];
rz(-2.0075802) q[3];
sx q[3];
rz(2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1535783) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(-2.4555092) q[0];
rz(-1.6185121) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(2.4783321) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2384757) q[0];
sx q[0];
rz(-0.65070063) q[0];
sx q[0];
rz(-1.531321) q[0];
x q[1];
rz(-2.5040607) q[2];
sx q[2];
rz(-1.7469284) q[2];
sx q[2];
rz(0.1917563) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.315553) q[1];
sx q[1];
rz(-1.7722881) q[1];
sx q[1];
rz(-1.0793346) q[1];
rz(-0.68526666) q[3];
sx q[3];
rz(-0.68260926) q[3];
sx q[3];
rz(-2.9311789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66701165) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(0.015080301) q[2];
rz(-2.1298501) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(-0.57141465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886154) q[0];
sx q[0];
rz(-0.73260728) q[0];
sx q[0];
rz(-3.0342039) q[0];
rz(0.30474162) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(2.6838141) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37516884) q[0];
sx q[0];
rz(-2.1918115) q[0];
sx q[0];
rz(0.39874886) q[0];
rz(-2.1019943) q[2];
sx q[2];
rz(-2.2441412) q[2];
sx q[2];
rz(-1.7216059) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6785504) q[1];
sx q[1];
rz(-1.926534) q[1];
sx q[1];
rz(-0.24865351) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55191603) q[3];
sx q[3];
rz(-1.8822) q[3];
sx q[3];
rz(0.50406721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7765939) q[2];
sx q[2];
rz(-1.4564617) q[2];
sx q[2];
rz(-1.4314852) q[2];
rz(-1.7163904) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(1.8553998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15437056) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(2.1881058) q[0];
rz(1.3257239) q[1];
sx q[1];
rz(-2.6486501) q[1];
sx q[1];
rz(-2.5568533) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7493233) q[0];
sx q[0];
rz(-1.4587147) q[0];
sx q[0];
rz(-2.6507676) q[0];
rz(-pi) q[1];
rz(-1.936475) q[2];
sx q[2];
rz(-1.056991) q[2];
sx q[2];
rz(-0.50525451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7331839) q[1];
sx q[1];
rz(-1.5015263) q[1];
sx q[1];
rz(0.44002156) q[1];
rz(-pi) q[2];
rz(-1.0869914) q[3];
sx q[3];
rz(-0.30011794) q[3];
sx q[3];
rz(-2.8407612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8087625) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(2.8653223) q[2];
rz(-0.55109465) q[3];
sx q[3];
rz(-1.3994183) q[3];
sx q[3];
rz(0.38366693) q[3];
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
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0693414) q[0];
sx q[0];
rz(-2.6273917) q[0];
sx q[0];
rz(0.65548354) q[0];
rz(-1.9845225) q[1];
sx q[1];
rz(-1.3815222) q[1];
sx q[1];
rz(-1.5225333) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039278395) q[0];
sx q[0];
rz(-0.73844665) q[0];
sx q[0];
rz(0.46955903) q[0];
x q[1];
rz(0.40752878) q[2];
sx q[2];
rz(-2.4113843) q[2];
sx q[2];
rz(-2.4804295) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.22776991) q[1];
sx q[1];
rz(-1.4064944) q[1];
sx q[1];
rz(2.5096202) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3823062) q[3];
sx q[3];
rz(-2.5873103) q[3];
sx q[3];
rz(0.62461531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15554252) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(0.84214169) q[2];
rz(1.5367674) q[3];
sx q[3];
rz(-1.3804881) q[3];
sx q[3];
rz(-0.026528927) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2355272) q[0];
sx q[0];
rz(-1.9530095) q[0];
sx q[0];
rz(-0.45146913) q[0];
rz(-2.4189667) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(-0.026635219) q[2];
sx q[2];
rz(-1.7799108) q[2];
sx q[2];
rz(-1.5540661) q[2];
rz(2.261653) q[3];
sx q[3];
rz(-1.6321245) q[3];
sx q[3];
rz(1.5499916) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];