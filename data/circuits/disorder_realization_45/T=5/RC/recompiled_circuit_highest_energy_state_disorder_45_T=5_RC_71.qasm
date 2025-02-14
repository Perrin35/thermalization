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
rz(1.6765321) q[0];
sx q[0];
rz(-2.9002011) q[0];
sx q[0];
rz(-0.13225947) q[0];
rz(-0.013068696) q[1];
sx q[1];
rz(-0.71813923) q[1];
sx q[1];
rz(3.1225966) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.976124) q[0];
sx q[0];
rz(-0.69688334) q[0];
sx q[0];
rz(-1.2319618) q[0];
rz(-1.7573518) q[2];
sx q[2];
rz(-1.6636408) q[2];
sx q[2];
rz(-1.873119) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4035114) q[1];
sx q[1];
rz(-0.36952239) q[1];
sx q[1];
rz(0.50003894) q[1];
rz(-pi) q[2];
rz(-1.8034987) q[3];
sx q[3];
rz(-1.4178172) q[3];
sx q[3];
rz(2.1683954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0211109) q[2];
sx q[2];
rz(-2.3349473) q[2];
sx q[2];
rz(-0.79180229) q[2];
rz(-1.0857438) q[3];
sx q[3];
rz(-1.9224242) q[3];
sx q[3];
rz(2.048548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87078142) q[0];
sx q[0];
rz(-2.9830611) q[0];
sx q[0];
rz(-2.812401) q[0];
rz(-1.5022494) q[1];
sx q[1];
rz(-0.90198016) q[1];
sx q[1];
rz(-2.7194068) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49445186) q[0];
sx q[0];
rz(-1.9943155) q[0];
sx q[0];
rz(-2.5472872) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3866502) q[2];
sx q[2];
rz(-2.1111541) q[2];
sx q[2];
rz(-0.042405142) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.364155) q[1];
sx q[1];
rz(-1.5244487) q[1];
sx q[1];
rz(0.32407659) q[1];
rz(0.77236891) q[3];
sx q[3];
rz(-1.4968781) q[3];
sx q[3];
rz(-2.071601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72055703) q[2];
sx q[2];
rz(-1.1853508) q[2];
sx q[2];
rz(-1.1400918) q[2];
rz(-3.1371878) q[3];
sx q[3];
rz(-1.5811698) q[3];
sx q[3];
rz(-3.0018023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7797101) q[0];
sx q[0];
rz(-3.0411868) q[0];
sx q[0];
rz(-2.7959339) q[0];
rz(-2.0480305) q[1];
sx q[1];
rz(-0.82229096) q[1];
sx q[1];
rz(-1.0543157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7525268) q[0];
sx q[0];
rz(-1.9574165) q[0];
sx q[0];
rz(-1.2091314) q[0];
rz(-pi) q[1];
rz(-1.6796175) q[2];
sx q[2];
rz(-0.76078712) q[2];
sx q[2];
rz(-1.6016122) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9413175) q[1];
sx q[1];
rz(-1.6345342) q[1];
sx q[1];
rz(0.45301183) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1183447) q[3];
sx q[3];
rz(-2.23684) q[3];
sx q[3];
rz(-2.9262528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67636079) q[2];
sx q[2];
rz(-2.5939442) q[2];
sx q[2];
rz(2.5993627) q[2];
rz(-2.9690361) q[3];
sx q[3];
rz(-1.4901284) q[3];
sx q[3];
rz(1.1057314) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8050103) q[0];
sx q[0];
rz(-1.7886826) q[0];
sx q[0];
rz(-1.2611058) q[0];
rz(-1.1359967) q[1];
sx q[1];
rz(-2.0295862) q[1];
sx q[1];
rz(-0.13066185) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4379556) q[0];
sx q[0];
rz(-1.2017631) q[0];
sx q[0];
rz(-0.28965182) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7950141) q[2];
sx q[2];
rz(-1.0750689) q[2];
sx q[2];
rz(0.6065281) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9798292) q[1];
sx q[1];
rz(-1.6024622) q[1];
sx q[1];
rz(-2.0791441) q[1];
rz(-pi) q[2];
rz(2.6819508) q[3];
sx q[3];
rz(-2.423624) q[3];
sx q[3];
rz(1.2470064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7589492) q[2];
sx q[2];
rz(-0.073315695) q[2];
sx q[2];
rz(-2.9200413) q[2];
rz(-2.4394636) q[3];
sx q[3];
rz(-1.5626855) q[3];
sx q[3];
rz(2.4041972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57109433) q[0];
sx q[0];
rz(-1.8966738) q[0];
sx q[0];
rz(-3.0076497) q[0];
rz(-0.36930034) q[1];
sx q[1];
rz(-1.1993273) q[1];
sx q[1];
rz(1.7514924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5969191) q[0];
sx q[0];
rz(-1.3256761) q[0];
sx q[0];
rz(2.0952203) q[0];
x q[1];
rz(0.26525396) q[2];
sx q[2];
rz(-1.2758288) q[2];
sx q[2];
rz(2.2504928) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0005324) q[1];
sx q[1];
rz(-0.62126011) q[1];
sx q[1];
rz(-0.48807524) q[1];
x q[2];
rz(-2.4485314) q[3];
sx q[3];
rz(-2.6821838) q[3];
sx q[3];
rz(2.1341679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8451346) q[2];
sx q[2];
rz(-0.90201169) q[2];
sx q[2];
rz(-1.5560537) q[2];
rz(2.8499917) q[3];
sx q[3];
rz(-0.4762989) q[3];
sx q[3];
rz(-2.8628023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45048243) q[0];
sx q[0];
rz(-2.1463558) q[0];
sx q[0];
rz(-0.52571785) q[0];
rz(-0.37824962) q[1];
sx q[1];
rz(-1.6098166) q[1];
sx q[1];
rz(0.27413109) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4286684) q[0];
sx q[0];
rz(-0.98085058) q[0];
sx q[0];
rz(1.601786) q[0];
rz(-pi) q[1];
rz(-1.1172764) q[2];
sx q[2];
rz(-1.1764865) q[2];
sx q[2];
rz(-2.0298983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8822599) q[1];
sx q[1];
rz(-1.8253528) q[1];
sx q[1];
rz(-1.1544636) q[1];
x q[2];
rz(1.317362) q[3];
sx q[3];
rz(-1.4259286) q[3];
sx q[3];
rz(2.1180017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5452925) q[2];
sx q[2];
rz(-1.9175074) q[2];
sx q[2];
rz(0.25823414) q[2];
rz(1.5958512) q[3];
sx q[3];
rz(-1.3629379) q[3];
sx q[3];
rz(0.405092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79893583) q[0];
sx q[0];
rz(-0.5181784) q[0];
sx q[0];
rz(-0.10547353) q[0];
rz(2.8666829) q[1];
sx q[1];
rz(-1.7173488) q[1];
sx q[1];
rz(-2.8555433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28562322) q[0];
sx q[0];
rz(-0.21846314) q[0];
sx q[0];
rz(-0.56665786) q[0];
rz(-pi) q[1];
rz(-1.8931383) q[2];
sx q[2];
rz(-2.3791056) q[2];
sx q[2];
rz(2.40138) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0512684) q[1];
sx q[1];
rz(-2.7503138) q[1];
sx q[1];
rz(-1.7404187) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55283847) q[3];
sx q[3];
rz(-1.9561319) q[3];
sx q[3];
rz(-0.60147731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1470571) q[2];
sx q[2];
rz(-0.68238443) q[2];
sx q[2];
rz(-0.4500173) q[2];
rz(2.5543645) q[3];
sx q[3];
rz(-1.8000032) q[3];
sx q[3];
rz(-1.2500866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29379544) q[0];
sx q[0];
rz(-1.4601409) q[0];
sx q[0];
rz(1.1867123) q[0];
rz(2.8222491) q[1];
sx q[1];
rz(-1.0217383) q[1];
sx q[1];
rz(2.879338) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866383) q[0];
sx q[0];
rz(-1.0269181) q[0];
sx q[0];
rz(-1.8094814) q[0];
rz(-pi) q[1];
rz(-1.8732583) q[2];
sx q[2];
rz(-2.710963) q[2];
sx q[2];
rz(-0.18948775) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1848534) q[1];
sx q[1];
rz(-1.0295086) q[1];
sx q[1];
rz(0.095603099) q[1];
rz(-pi) q[2];
rz(-2.6899509) q[3];
sx q[3];
rz(-1.181895) q[3];
sx q[3];
rz(1.9967784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.78476) q[2];
sx q[2];
rz(-2.6748952) q[2];
sx q[2];
rz(-2.812815) q[2];
rz(1.5152991) q[3];
sx q[3];
rz(-1.7833775) q[3];
sx q[3];
rz(-2.7330107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3560155) q[0];
sx q[0];
rz(-0.51078904) q[0];
sx q[0];
rz(0.27895862) q[0];
rz(-0.51697671) q[1];
sx q[1];
rz(-0.37716436) q[1];
sx q[1];
rz(-1.4252211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036439786) q[0];
sx q[0];
rz(-0.62314674) q[0];
sx q[0];
rz(0.07499174) q[0];
rz(-pi) q[1];
rz(0.65966971) q[2];
sx q[2];
rz(-2.2222509) q[2];
sx q[2];
rz(-1.4209117) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8502073) q[1];
sx q[1];
rz(-1.3818023) q[1];
sx q[1];
rz(-0.50598617) q[1];
rz(-pi) q[2];
rz(-1.2012977) q[3];
sx q[3];
rz(-2.4505121) q[3];
sx q[3];
rz(1.7352417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4031389) q[2];
sx q[2];
rz(-2.8454915) q[2];
sx q[2];
rz(-0.1599172) q[2];
rz(-0.81965172) q[3];
sx q[3];
rz(-1.9194226) q[3];
sx q[3];
rz(-2.405449) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6920456) q[0];
sx q[0];
rz(-1.291438) q[0];
sx q[0];
rz(-0.29577574) q[0];
rz(-2.6462818) q[1];
sx q[1];
rz(-0.43423978) q[1];
sx q[1];
rz(-1.5089418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6497632) q[0];
sx q[0];
rz(-1.9252281) q[0];
sx q[0];
rz(-1.9948122) q[0];
rz(1.7118218) q[2];
sx q[2];
rz(-2.786028) q[2];
sx q[2];
rz(-0.18294683) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7612696) q[1];
sx q[1];
rz(-1.4810303) q[1];
sx q[1];
rz(-1.5029546) q[1];
x q[2];
rz(2.6466683) q[3];
sx q[3];
rz(-1.2240181) q[3];
sx q[3];
rz(2.1174255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32790023) q[2];
sx q[2];
rz(-1.6715965) q[2];
sx q[2];
rz(2.5076765) q[2];
rz(-0.085112326) q[3];
sx q[3];
rz(-1.2462933) q[3];
sx q[3];
rz(-2.3688721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2278628) q[0];
sx q[0];
rz(-1.5416523) q[0];
sx q[0];
rz(-0.1097485) q[0];
rz(1.2211424) q[1];
sx q[1];
rz(-2.2962062) q[1];
sx q[1];
rz(-2.5795945) q[1];
rz(-1.132125) q[2];
sx q[2];
rz(-1.7747468) q[2];
sx q[2];
rz(-2.940831) q[2];
rz(1.4681592) q[3];
sx q[3];
rz(-2.9087421) q[3];
sx q[3];
rz(-1.4431492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
