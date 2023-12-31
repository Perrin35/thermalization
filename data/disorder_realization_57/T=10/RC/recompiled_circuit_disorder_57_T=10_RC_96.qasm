OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6141619) q[0];
sx q[0];
rz(-0.80978137) q[0];
sx q[0];
rz(2.6101987) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(-1.7778492) q[1];
sx q[1];
rz(-1.2030503) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1019842) q[0];
sx q[0];
rz(-1.3586205) q[0];
sx q[0];
rz(2.9456445) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7877903) q[2];
sx q[2];
rz(-1.8197682) q[2];
sx q[2];
rz(2.4169902) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0002999) q[1];
sx q[1];
rz(-1.3360026) q[1];
sx q[1];
rz(1.512212) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0193127) q[3];
sx q[3];
rz(-0.29502007) q[3];
sx q[3];
rz(-0.89000851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27665859) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(2.9425088) q[2];
rz(1.6137326) q[3];
sx q[3];
rz(-0.49694967) q[3];
sx q[3];
rz(-0.73408192) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99877015) q[0];
sx q[0];
rz(-3.0144189) q[0];
sx q[0];
rz(2.3221827) q[0];
rz(2.8557414) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(-1.2664638) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768515) q[0];
sx q[0];
rz(-0.74155945) q[0];
sx q[0];
rz(0.57912796) q[0];
x q[1];
rz(1.0533633) q[2];
sx q[2];
rz(-2.6839163) q[2];
sx q[2];
rz(0.44331726) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30217583) q[1];
sx q[1];
rz(-1.9315047) q[1];
sx q[1];
rz(1.9990666) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1191788) q[3];
sx q[3];
rz(-2.0997117) q[3];
sx q[3];
rz(0.12714735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1566029) q[2];
sx q[2];
rz(-1.9282324) q[2];
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
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53482985) q[0];
sx q[0];
rz(-2.02089) q[0];
sx q[0];
rz(-0.30360046) q[0];
rz(-1.3820232) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(0.64750013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8436463) q[0];
sx q[0];
rz(-2.2046411) q[0];
sx q[0];
rz(-1.6266536) q[0];
rz(2.5325003) q[2];
sx q[2];
rz(-2.3272772) q[2];
sx q[2];
rz(-2.7012205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69933575) q[1];
sx q[1];
rz(-1.2067716) q[1];
sx q[1];
rz(-2.060021) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87326761) q[3];
sx q[3];
rz(-2.0487361) q[3];
sx q[3];
rz(-2.0877116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2993762) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(-1.5608609) q[2];
rz(0.98172274) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(2.4981807) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2340387) q[0];
sx q[0];
rz(-0.83775318) q[0];
sx q[0];
rz(0.87189829) q[0];
rz(0.38726989) q[1];
sx q[1];
rz(-0.64278066) q[1];
sx q[1];
rz(0.29104582) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29845333) q[0];
sx q[0];
rz(-1.342134) q[0];
sx q[0];
rz(0.97897162) q[0];
rz(-pi) q[1];
rz(-0.29291885) q[2];
sx q[2];
rz(-1.2625164) q[2];
sx q[2];
rz(2.3777547) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9705022) q[1];
sx q[1];
rz(-2.2541752) q[1];
sx q[1];
rz(-2.469953) q[1];
x q[2];
rz(-1.325874) q[3];
sx q[3];
rz(-0.59026736) q[3];
sx q[3];
rz(0.11412379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7724093) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(2.5975361) q[2];
rz(2.6323075) q[3];
sx q[3];
rz(-2.0777168) q[3];
sx q[3];
rz(0.86597401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0020224) q[0];
sx q[0];
rz(-2.1152571) q[0];
sx q[0];
rz(0.64506662) q[0];
rz(-0.6257261) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(-0.63017875) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7932927) q[0];
sx q[0];
rz(-0.59069809) q[0];
sx q[0];
rz(2.9212055) q[0];
rz(-0.39406392) q[2];
sx q[2];
rz(-2.2820018) q[2];
sx q[2];
rz(1.1472536) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4634134) q[1];
sx q[1];
rz(-0.58600512) q[1];
sx q[1];
rz(1.1581808) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7230677) q[3];
sx q[3];
rz(-2.1232743) q[3];
sx q[3];
rz(-0.0020364062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3877635) q[2];
sx q[2];
rz(-1.8176273) q[2];
sx q[2];
rz(1.7774263) q[2];
rz(1.2498614) q[3];
sx q[3];
rz(-1.2156237) q[3];
sx q[3];
rz(-0.92145872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0867778) q[0];
sx q[0];
rz(-0.58371109) q[0];
sx q[0];
rz(0.71682799) q[0];
rz(-0.43771276) q[1];
sx q[1];
rz(-0.68044674) q[1];
sx q[1];
rz(-2.3278918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87156224) q[0];
sx q[0];
rz(-0.22261482) q[0];
sx q[0];
rz(-1.8588258) q[0];
rz(-2.4196163) q[2];
sx q[2];
rz(-2.0275896) q[2];
sx q[2];
rz(-1.4095838) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70255792) q[1];
sx q[1];
rz(-0.82779373) q[1];
sx q[1];
rz(-2.031416) q[1];
rz(0.86088647) q[3];
sx q[3];
rz(-0.45300278) q[3];
sx q[3];
rz(-0.83356726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8852691) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(1.9539333) q[2];
rz(-2.2096283) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1535783) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(0.6860835) q[0];
rz(1.5230806) q[1];
sx q[1];
rz(-0.97421092) q[1];
sx q[1];
rz(0.66326052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8535186) q[0];
sx q[0];
rz(-0.92068866) q[0];
sx q[0];
rz(0.030036146) q[0];
rz(0.29055039) q[2];
sx q[2];
rz(-2.4834589) q[2];
sx q[2];
rz(-1.611329) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82603964) q[1];
sx q[1];
rz(-1.3693046) q[1];
sx q[1];
rz(1.0793346) q[1];
rz(1.0955986) q[3];
sx q[3];
rz(-1.0605269) q[3];
sx q[3];
rz(2.1197532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.474581) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(-1.0117426) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(-2.570178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(2.836851) q[1];
sx q[1];
rz(-1.3184897) q[1];
sx q[1];
rz(2.6838141) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4360355) q[0];
sx q[0];
rz(-1.8920664) q[0];
sx q[0];
rz(0.9106439) q[0];
rz(-pi) q[1];
rz(2.5758178) q[2];
sx q[2];
rz(-0.83116311) q[2];
sx q[2];
rz(0.6643675) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0486054) q[1];
sx q[1];
rz(-0.43097207) q[1];
sx q[1];
rz(0.98577568) q[1];
rz(-pi) q[2];
rz(-1.9321953) q[3];
sx q[3];
rz(-2.0933588) q[3];
sx q[3];
rz(-0.88025974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7765939) q[2];
sx q[2];
rz(-1.4564617) q[2];
sx q[2];
rz(1.4314852) q[2];
rz(-1.7163904) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(-1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9872221) q[0];
sx q[0];
rz(-0.4796589) q[0];
sx q[0];
rz(2.1881058) q[0];
rz(-1.8158688) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(-0.58473933) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23823243) q[0];
sx q[0];
rz(-2.0582709) q[0];
sx q[0];
rz(-1.6977298) q[0];
rz(-pi) q[1];
x q[1];
rz(2.598001) q[2];
sx q[2];
rz(-1.2541176) q[2];
sx q[2];
rz(1.8900332) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1251724) q[1];
sx q[1];
rz(-0.44508815) q[1];
sx q[1];
rz(-0.16146407) q[1];
rz(2.9986266) q[3];
sx q[3];
rz(-1.8355832) q[3];
sx q[3];
rz(0.20204443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3328302) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(2.8653223) q[2];
rz(0.55109465) q[3];
sx q[3];
rz(-1.3994183) q[3];
sx q[3];
rz(2.7579257) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0693414) q[0];
sx q[0];
rz(-2.6273917) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(1.1570702) q[1];
sx q[1];
rz(-1.3815222) q[1];
sx q[1];
rz(-1.5225333) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1023143) q[0];
sx q[0];
rz(-2.403146) q[0];
sx q[0];
rz(-2.6720336) q[0];
rz(-pi) q[1];
rz(-1.9117781) q[2];
sx q[2];
rz(-2.2298862) q[2];
sx q[2];
rz(-1.9552719) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5628221) q[1];
sx q[1];
rz(-2.4914503) q[1];
sx q[1];
rz(2.8679718) q[1];
rz(-pi) q[2];
rz(-2.3823062) q[3];
sx q[3];
rz(-2.5873103) q[3];
sx q[3];
rz(2.5169773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9860501) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(0.84214169) q[2];
rz(1.6048253) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(-0.026528927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90606541) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(-2.4189667) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(3.1149574) q[2];
sx q[2];
rz(-1.7799108) q[2];
sx q[2];
rz(-1.5540661) q[2];
rz(-3.0620861) q[3];
sx q[3];
rz(-0.88149298) q[3];
sx q[3];
rz(0.029824921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
