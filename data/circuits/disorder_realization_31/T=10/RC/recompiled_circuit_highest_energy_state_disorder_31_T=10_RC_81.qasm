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
rz(1.0072768) q[0];
sx q[0];
rz(-0.64581031) q[0];
sx q[0];
rz(0.3663775) q[0];
rz(0.74908787) q[1];
sx q[1];
rz(-1.5146511) q[1];
sx q[1];
rz(-0.14263022) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8170022) q[0];
sx q[0];
rz(-2.3329371) q[0];
sx q[0];
rz(2.3456744) q[0];
rz(0.5515302) q[2];
sx q[2];
rz(-0.26187927) q[2];
sx q[2];
rz(2.8079845) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4971338) q[1];
sx q[1];
rz(-0.49127775) q[1];
sx q[1];
rz(1.3158448) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5589988) q[3];
sx q[3];
rz(-2.5911463) q[3];
sx q[3];
rz(-2.1728476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7819405) q[2];
sx q[2];
rz(-2.6617229) q[2];
sx q[2];
rz(1.5894319) q[2];
rz(1.3945256) q[3];
sx q[3];
rz(-2.7701869) q[3];
sx q[3];
rz(-0.68215218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8438016) q[0];
sx q[0];
rz(-0.60972917) q[0];
sx q[0];
rz(0.30526701) q[0];
rz(1.3097395) q[1];
sx q[1];
rz(-0.42116183) q[1];
sx q[1];
rz(3.0895244) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16927707) q[0];
sx q[0];
rz(-2.9575363) q[0];
sx q[0];
rz(2.4993624) q[0];
x q[1];
rz(2.4737055) q[2];
sx q[2];
rz(-2.6781262) q[2];
sx q[2];
rz(-0.98834044) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0525103) q[1];
sx q[1];
rz(-1.4064612) q[1];
sx q[1];
rz(1.8402064) q[1];
x q[2];
rz(3.0092485) q[3];
sx q[3];
rz(-2.2346046) q[3];
sx q[3];
rz(-2.9470766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0569968) q[2];
sx q[2];
rz(-1.3330385) q[2];
sx q[2];
rz(-1.4196654) q[2];
rz(-2.2325884) q[3];
sx q[3];
rz(-0.82908583) q[3];
sx q[3];
rz(-1.7349294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1664355) q[0];
sx q[0];
rz(-1.1844013) q[0];
sx q[0];
rz(-1.3982406) q[0];
rz(2.1947529) q[1];
sx q[1];
rz(-2.0473174) q[1];
sx q[1];
rz(-1.2907226) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4104267) q[0];
sx q[0];
rz(-1.7058347) q[0];
sx q[0];
rz(0.076083994) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3554736) q[2];
sx q[2];
rz(-1.2183471) q[2];
sx q[2];
rz(0.64482433) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0086338) q[1];
sx q[1];
rz(-1.5040288) q[1];
sx q[1];
rz(0.029032727) q[1];
rz(-2.7644058) q[3];
sx q[3];
rz(-0.50807774) q[3];
sx q[3];
rz(2.7722904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1281841) q[2];
sx q[2];
rz(-2.7161697) q[2];
sx q[2];
rz(-0.95995963) q[2];
rz(-0.32402447) q[3];
sx q[3];
rz(-2.6257381) q[3];
sx q[3];
rz(-2.7014151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1337166) q[0];
sx q[0];
rz(-0.77819264) q[0];
sx q[0];
rz(2.9414951) q[0];
rz(-2.5307185) q[1];
sx q[1];
rz(-2.4433544) q[1];
sx q[1];
rz(-2.3470338) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.824423) q[0];
sx q[0];
rz(-0.18704913) q[0];
sx q[0];
rz(-0.6424128) q[0];
rz(-pi) q[1];
x q[1];
rz(1.346454) q[2];
sx q[2];
rz(-2.4560929) q[2];
sx q[2];
rz(-2.3239674) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.17758372) q[1];
sx q[1];
rz(-1.4191886) q[1];
sx q[1];
rz(1.8880009) q[1];
rz(-pi) q[2];
rz(0.67312397) q[3];
sx q[3];
rz(-0.84852058) q[3];
sx q[3];
rz(-0.22658843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8841298) q[2];
sx q[2];
rz(-2.5778759) q[2];
sx q[2];
rz(1.5719315) q[2];
rz(-1.7852768) q[3];
sx q[3];
rz(-1.1607728) q[3];
sx q[3];
rz(-0.13264382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86303478) q[0];
sx q[0];
rz(-2.886241) q[0];
sx q[0];
rz(-0.72364664) q[0];
rz(0.10069314) q[1];
sx q[1];
rz(-1.0483402) q[1];
sx q[1];
rz(0.06631276) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4935987) q[0];
sx q[0];
rz(-2.1382522) q[0];
sx q[0];
rz(-0.074851224) q[0];
x q[1];
rz(2.5939717) q[2];
sx q[2];
rz(-0.98973083) q[2];
sx q[2];
rz(-0.18085322) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4206754) q[1];
sx q[1];
rz(-2.2128792) q[1];
sx q[1];
rz(-2.9319619) q[1];
rz(-pi) q[2];
rz(-0.67541404) q[3];
sx q[3];
rz(-0.49478045) q[3];
sx q[3];
rz(0.29008415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.940332) q[2];
sx q[2];
rz(-1.9278434) q[2];
sx q[2];
rz(-2.9230389) q[2];
rz(-0.39441937) q[3];
sx q[3];
rz(-0.68587488) q[3];
sx q[3];
rz(0.39943892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047693096) q[0];
sx q[0];
rz(-2.221929) q[0];
sx q[0];
rz(3.1374875) q[0];
rz(-0.63402367) q[1];
sx q[1];
rz(-1.5019633) q[1];
sx q[1];
rz(2.1852469) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8391181) q[0];
sx q[0];
rz(-1.7540521) q[0];
sx q[0];
rz(1.5292581) q[0];
rz(-pi) q[1];
rz(2.4615131) q[2];
sx q[2];
rz(-2.5181327) q[2];
sx q[2];
rz(2.4768444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9117212) q[1];
sx q[1];
rz(-1.8485862) q[1];
sx q[1];
rz(-1.5650657) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6878674) q[3];
sx q[3];
rz(-1.9170951) q[3];
sx q[3];
rz(0.1629783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6213106) q[2];
sx q[2];
rz(-1.7063528) q[2];
sx q[2];
rz(-2.3069646) q[2];
rz(-2.7708715) q[3];
sx q[3];
rz(-0.70086896) q[3];
sx q[3];
rz(2.4247775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-1.9653559) q[0];
sx q[0];
rz(-3.1270202) q[0];
sx q[0];
rz(0.45005774) q[0];
rz(-0.4484446) q[1];
sx q[1];
rz(-2.3095755) q[1];
sx q[1];
rz(-0.73743302) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56313228) q[0];
sx q[0];
rz(-1.1473343) q[0];
sx q[0];
rz(2.9388301) q[0];
rz(3.1011821) q[2];
sx q[2];
rz(-1.7120265) q[2];
sx q[2];
rz(0.34858957) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0998675) q[1];
sx q[1];
rz(-1.2507572) q[1];
sx q[1];
rz(-0.36497023) q[1];
x q[2];
rz(2.8705863) q[3];
sx q[3];
rz(-2.3428863) q[3];
sx q[3];
rz(0.17548156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2927148) q[2];
sx q[2];
rz(-0.60056168) q[2];
sx q[2];
rz(2.4265477) q[2];
rz(2.6333366) q[3];
sx q[3];
rz(-1.4628937) q[3];
sx q[3];
rz(-1.3983294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7364863) q[0];
sx q[0];
rz(-2.5953601) q[0];
sx q[0];
rz(0.35838321) q[0];
rz(2.7690923) q[1];
sx q[1];
rz(-1.7820216) q[1];
sx q[1];
rz(0.43982664) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43482237) q[0];
sx q[0];
rz(-2.3447038) q[0];
sx q[0];
rz(-3.0838941) q[0];
rz(-pi) q[1];
rz(-0.8586712) q[2];
sx q[2];
rz(-2.5617848) q[2];
sx q[2];
rz(-0.64263868) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3258626) q[1];
sx q[1];
rz(-1.3943895) q[1];
sx q[1];
rz(2.8730064) q[1];
rz(-0.16086264) q[3];
sx q[3];
rz(-0.32718143) q[3];
sx q[3];
rz(0.88674816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0899352) q[2];
sx q[2];
rz(-3.0445485) q[2];
sx q[2];
rz(-2.4329321) q[2];
rz(2.8900201) q[3];
sx q[3];
rz(-2.1422062) q[3];
sx q[3];
rz(0.034041762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2947023) q[0];
sx q[0];
rz(-1.0346233) q[0];
sx q[0];
rz(2.5501472) q[0];
rz(3.1239608) q[1];
sx q[1];
rz(-2.089274) q[1];
sx q[1];
rz(-3.0406521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936241) q[0];
sx q[0];
rz(-1.7059479) q[0];
sx q[0];
rz(-1.253932) q[0];
rz(-2.5048961) q[2];
sx q[2];
rz(-1.6023738) q[2];
sx q[2];
rz(-0.65177887) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.832022) q[1];
sx q[1];
rz(-2.8831318) q[1];
sx q[1];
rz(2.0548106) q[1];
rz(-pi) q[2];
rz(-1.798257) q[3];
sx q[3];
rz(-0.69256562) q[3];
sx q[3];
rz(2.8253959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.65323222) q[2];
sx q[2];
rz(-2.4182352) q[2];
sx q[2];
rz(-1.6548659) q[2];
rz(-0.15445736) q[3];
sx q[3];
rz(-1.1052701) q[3];
sx q[3];
rz(2.7753593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3514997) q[0];
sx q[0];
rz(-0.14227754) q[0];
sx q[0];
rz(-1.0737786) q[0];
rz(1.0723266) q[1];
sx q[1];
rz(-2.8876979) q[1];
sx q[1];
rz(3.0825739) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85262839) q[0];
sx q[0];
rz(-1.5375332) q[0];
sx q[0];
rz(0.19927523) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8779936) q[2];
sx q[2];
rz(-2.0661743) q[2];
sx q[2];
rz(1.4857298) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4764875) q[1];
sx q[1];
rz(-0.8554994) q[1];
sx q[1];
rz(2.4060898) q[1];
rz(-pi) q[2];
rz(-2.288184) q[3];
sx q[3];
rz(-2.0730258) q[3];
sx q[3];
rz(2.2706007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8129639) q[2];
sx q[2];
rz(-1.2920516) q[2];
sx q[2];
rz(2.9126419) q[2];
rz(-1.1936584) q[3];
sx q[3];
rz(-0.3242068) q[3];
sx q[3];
rz(1.6875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-3.1345632) q[0];
sx q[0];
rz(-2.3352191) q[0];
sx q[0];
rz(2.4051608) q[0];
rz(-1.3855343) q[1];
sx q[1];
rz(-1.3722739) q[1];
sx q[1];
rz(-1.3442232) q[1];
rz(1.6177542) q[2];
sx q[2];
rz(-2.6018819) q[2];
sx q[2];
rz(0.56503208) q[2];
rz(0.61290716) q[3];
sx q[3];
rz(-1.404543) q[3];
sx q[3];
rz(2.3096655) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
