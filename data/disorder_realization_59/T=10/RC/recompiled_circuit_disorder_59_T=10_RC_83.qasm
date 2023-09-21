OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(-0.78200114) q[0];
sx q[0];
rz(1.8703823) q[0];
rz(-2.864569) q[1];
sx q[1];
rz(-2.6695873) q[1];
sx q[1];
rz(-0.0013874887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9936375) q[0];
sx q[0];
rz(-1.7339098) q[0];
sx q[0];
rz(1.1975343) q[0];
rz(0.087287993) q[2];
sx q[2];
rz(-2.6929571) q[2];
sx q[2];
rz(2.0729614) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.063349799) q[1];
sx q[1];
rz(-1.500505) q[1];
sx q[1];
rz(0.49777828) q[1];
rz(-pi) q[2];
rz(2.1625159) q[3];
sx q[3];
rz(-0.43833971) q[3];
sx q[3];
rz(-2.0548267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15443054) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(-1.0162214) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(-2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7063023) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(-0.97066561) q[0];
rz(2.1043815) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(-2.326139) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97074189) q[0];
sx q[0];
rz(-0.65119699) q[0];
sx q[0];
rz(-3.0424776) q[0];
rz(-2.8380727) q[2];
sx q[2];
rz(-0.75280658) q[2];
sx q[2];
rz(0.34740651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.193589) q[1];
sx q[1];
rz(-1.5352328) q[1];
sx q[1];
rz(0.30102945) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29173298) q[3];
sx q[3];
rz(-2.4749304) q[3];
sx q[3];
rz(3.047903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4619535) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(-1.1535545) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3011424) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(-0.87483037) q[0];
rz(-1.3300928) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(0.99951807) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777269) q[0];
sx q[0];
rz(-0.3361055) q[0];
sx q[0];
rz(-2.0396114) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3109762) q[2];
sx q[2];
rz(-2.3306371) q[2];
sx q[2];
rz(-0.22731552) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0083025) q[1];
sx q[1];
rz(-1.3160719) q[1];
sx q[1];
rz(2.6363274) q[1];
rz(-2.0796892) q[3];
sx q[3];
rz(-1.5012) q[3];
sx q[3];
rz(-0.059046179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(-2.5615454) q[2];
rz(-0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.9918611) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(2.2312009) q[0];
rz(-2.6903649) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(-0.27483637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4764458) q[0];
sx q[0];
rz(-3.0229212) q[0];
sx q[0];
rz(0.84400405) q[0];
x q[1];
rz(0.0015953548) q[2];
sx q[2];
rz(-2.0312107) q[2];
sx q[2];
rz(2.7574725) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53510016) q[1];
sx q[1];
rz(-2.1537158) q[1];
sx q[1];
rz(-2.5938354) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12675385) q[3];
sx q[3];
rz(-2.2634014) q[3];
sx q[3];
rz(-0.17514378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3466907) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(-1.9968962) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(0.16170734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-0.65790025) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(2.143798) q[0];
rz(-2.9580341) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(1.6246187) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.330634) q[0];
sx q[0];
rz(-1.2067544) q[0];
sx q[0];
rz(2.2383658) q[0];
rz(-0.60913182) q[2];
sx q[2];
rz(-2.5172148) q[2];
sx q[2];
rz(0.49027157) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9285674) q[1];
sx q[1];
rz(-1.9145609) q[1];
sx q[1];
rz(-2.7108971) q[1];
rz(0.59668221) q[3];
sx q[3];
rz(-2.0697069) q[3];
sx q[3];
rz(2.1614206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8020442) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(-2.8175763) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-0.75606212) q[3];
sx q[3];
rz(1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.3762387) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(0.34067672) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37939385) q[0];
sx q[0];
rz(-1.6025935) q[0];
sx q[0];
rz(1.5048774) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2773877) q[2];
sx q[2];
rz(-0.83653203) q[2];
sx q[2];
rz(-2.4794527) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1250455) q[1];
sx q[1];
rz(-0.53495896) q[1];
sx q[1];
rz(-0.46334456) q[1];
rz(-0.42231456) q[3];
sx q[3];
rz(-1.1186244) q[3];
sx q[3];
rz(2.5715668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1910151) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(0.77159709) q[2];
rz(-0.54780444) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(-0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(0.46494928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9879887) q[0];
sx q[0];
rz(-1.1055595) q[0];
sx q[0];
rz(-2.9265755) q[0];
rz(0.49352383) q[2];
sx q[2];
rz(-2.0888622) q[2];
sx q[2];
rz(-1.7871737) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9629434) q[1];
sx q[1];
rz(-1.6766251) q[1];
sx q[1];
rz(-0.21957285) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5246478) q[3];
sx q[3];
rz(-1.5966291) q[3];
sx q[3];
rz(1.3290562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1376301) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(0.31759343) q[2];
rz(-2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5193609) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(0.28433329) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5057482) q[0];
sx q[0];
rz(-2.3730179) q[0];
sx q[0];
rz(1.4710674) q[0];
rz(1.3252844) q[2];
sx q[2];
rz(-1.3405372) q[2];
sx q[2];
rz(0.85207176) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17532149) q[1];
sx q[1];
rz(-0.92370719) q[1];
sx q[1];
rz(-0.21827571) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4140698) q[3];
sx q[3];
rz(-1.6046451) q[3];
sx q[3];
rz(-2.1050997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4006965) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(-0.25137869) q[2];
rz(2.5583983) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(1.4096227) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(-0.92357606) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(2.267568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4955935) q[0];
sx q[0];
rz(-0.80276239) q[0];
sx q[0];
rz(-1.5828703) q[0];
x q[1];
rz(1.1633515) q[2];
sx q[2];
rz(-1.2375087) q[2];
sx q[2];
rz(2.8826706) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6898432) q[1];
sx q[1];
rz(-0.50536957) q[1];
sx q[1];
rz(1.4228574) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60871082) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(-1.1876719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.727227) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(1.184458) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(-1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0062362) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(0.7243048) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(-1.9627409) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2154685) q[0];
sx q[0];
rz(-1.3660396) q[0];
sx q[0];
rz(1.8490851) q[0];
rz(-pi) q[1];
rz(0.63209052) q[2];
sx q[2];
rz(-0.1569911) q[2];
sx q[2];
rz(2.0096411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.19503838) q[1];
sx q[1];
rz(-2.185501) q[1];
sx q[1];
rz(0.22919319) q[1];
x q[2];
rz(-2.7614818) q[3];
sx q[3];
rz(-0.5553402) q[3];
sx q[3];
rz(1.4243319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.05802352) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(-0.87456885) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.6806867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(1.6451251) q[2];
sx q[2];
rz(-2.7004514) q[2];
sx q[2];
rz(-2.2218291) q[2];
rz(-1.7189797) q[3];
sx q[3];
rz(-1.2978745) q[3];
sx q[3];
rz(0.10980448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
