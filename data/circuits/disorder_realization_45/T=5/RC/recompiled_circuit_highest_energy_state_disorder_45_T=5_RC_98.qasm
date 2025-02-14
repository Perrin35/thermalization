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
rz(2.4234534) q[1];
sx q[1];
rz(9.443774) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16546862) q[0];
sx q[0];
rz(-0.69688334) q[0];
sx q[0];
rz(-1.9096309) q[0];
rz(-pi) q[1];
rz(-2.0360581) q[2];
sx q[2];
rz(-0.20814116) q[2];
sx q[2];
rz(-0.15424745) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.87352301) q[1];
sx q[1];
rz(-1.8933081) q[1];
sx q[1];
rz(1.7544062) q[1];
x q[2];
rz(2.1601342) q[3];
sx q[3];
rz(-0.27772003) q[3];
sx q[3];
rz(3.115417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12048177) q[2];
sx q[2];
rz(-0.80664539) q[2];
sx q[2];
rz(-0.79180229) q[2];
rz(-2.0558489) q[3];
sx q[3];
rz(-1.2191685) q[3];
sx q[3];
rz(-1.0930446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2708112) q[0];
sx q[0];
rz(-0.15853156) q[0];
sx q[0];
rz(-2.812401) q[0];
rz(1.5022494) q[1];
sx q[1];
rz(-2.2396125) q[1];
sx q[1];
rz(0.42218581) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6471408) q[0];
sx q[0];
rz(-1.9943155) q[0];
sx q[0];
rz(-2.5472872) q[0];
x q[1];
rz(0.29624002) q[2];
sx q[2];
rz(-0.56791091) q[2];
sx q[2];
rz(-2.8366249) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77743769) q[1];
sx q[1];
rz(-1.6171439) q[1];
sx q[1];
rz(-0.32407659) q[1];
rz(-pi) q[2];
rz(1.6738191) q[3];
sx q[3];
rz(-0.80108445) q[3];
sx q[3];
rz(0.57263206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4210356) q[2];
sx q[2];
rz(-1.1853508) q[2];
sx q[2];
rz(-2.0015008) q[2];
rz(-0.0044048443) q[3];
sx q[3];
rz(-1.5811698) q[3];
sx q[3];
rz(-0.13979039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36188257) q[0];
sx q[0];
rz(-0.10040586) q[0];
sx q[0];
rz(-0.34565872) q[0];
rz(-2.0480305) q[1];
sx q[1];
rz(-0.82229096) q[1];
sx q[1];
rz(2.0872769) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525268) q[0];
sx q[0];
rz(-1.9574165) q[0];
sx q[0];
rz(1.2091314) q[0];
rz(0.81297154) q[2];
sx q[2];
rz(-1.6457498) q[2];
sx q[2];
rz(3.093442) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40151628) q[1];
sx q[1];
rz(-1.1187727) q[1];
sx q[1];
rz(1.4999309) q[1];
rz(-3.1183447) q[3];
sx q[3];
rz(-2.23684) q[3];
sx q[3];
rz(-2.9262528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4652319) q[2];
sx q[2];
rz(-2.5939442) q[2];
sx q[2];
rz(2.5993627) q[2];
rz(-2.9690361) q[3];
sx q[3];
rz(-1.6514643) q[3];
sx q[3];
rz(-1.1057314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8050103) q[0];
sx q[0];
rz(-1.35291) q[0];
sx q[0];
rz(1.8804869) q[0];
rz(-1.1359967) q[1];
sx q[1];
rz(-1.1120064) q[1];
sx q[1];
rz(0.13066185) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9016584) q[0];
sx q[0];
rz(-1.3011509) q[0];
sx q[0];
rz(-1.95437) q[0];
x q[1];
rz(-2.7515) q[2];
sx q[2];
rz(-2.60139) q[2];
sx q[2];
rz(0.15946968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16176341) q[1];
sx q[1];
rz(-1.5391304) q[1];
sx q[1];
rz(-2.0791441) q[1];
rz(-pi) q[2];
rz(-1.9404802) q[3];
sx q[3];
rz(-2.2013328) q[3];
sx q[3];
rz(-0.66555221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7589492) q[2];
sx q[2];
rz(-0.073315695) q[2];
sx q[2];
rz(2.9200413) q[2];
rz(2.4394636) q[3];
sx q[3];
rz(-1.5626855) q[3];
sx q[3];
rz(-2.4041972) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57109433) q[0];
sx q[0];
rz(-1.2449188) q[0];
sx q[0];
rz(3.0076497) q[0];
rz(0.36930034) q[1];
sx q[1];
rz(-1.9422653) q[1];
sx q[1];
rz(-1.3901002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62890118) q[0];
sx q[0];
rz(-2.5675964) q[0];
sx q[0];
rz(2.0341134) q[0];
x q[1];
rz(-0.85890074) q[2];
sx q[2];
rz(-2.7475069) q[2];
sx q[2];
rz(-1.6426298) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.83727969) q[1];
sx q[1];
rz(-1.8472478) q[1];
sx q[1];
rz(-2.5778092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36365328) q[3];
sx q[3];
rz(-1.2835652) q[3];
sx q[3];
rz(-1.2032697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.2964581) q[2];
sx q[2];
rz(-0.90201169) q[2];
sx q[2];
rz(1.5560537) q[2];
rz(-2.8499917) q[3];
sx q[3];
rz(-0.4762989) q[3];
sx q[3];
rz(-0.27879032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6911102) q[0];
sx q[0];
rz(-0.99523681) q[0];
sx q[0];
rz(-0.52571785) q[0];
rz(2.763343) q[1];
sx q[1];
rz(-1.6098166) q[1];
sx q[1];
rz(-2.8674616) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0167086) q[0];
sx q[0];
rz(-1.5965466) q[0];
sx q[0];
rz(-2.5514249) q[0];
x q[1];
rz(2.3304105) q[2];
sx q[2];
rz(-2.5497782) q[2];
sx q[2];
rz(2.9331911) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3124481) q[1];
sx q[1];
rz(-0.48407468) q[1];
sx q[1];
rz(-2.1425254) q[1];
rz(-pi) q[2];
rz(-0.14957803) q[3];
sx q[3];
rz(-1.8215186) q[3];
sx q[3];
rz(0.58457812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5452925) q[2];
sx q[2];
rz(-1.2240852) q[2];
sx q[2];
rz(-0.25823414) q[2];
rz(-1.5457414) q[3];
sx q[3];
rz(-1.3629379) q[3];
sx q[3];
rz(-2.7365007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3426568) q[0];
sx q[0];
rz(-2.6234143) q[0];
sx q[0];
rz(3.0361191) q[0];
rz(2.8666829) q[1];
sx q[1];
rz(-1.4242438) q[1];
sx q[1];
rz(-0.28604937) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8559694) q[0];
sx q[0];
rz(-0.21846314) q[0];
sx q[0];
rz(0.56665786) q[0];
rz(-pi) q[1];
rz(-0.83467612) q[2];
sx q[2];
rz(-1.7913941) q[2];
sx q[2];
rz(2.5479864) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0512684) q[1];
sx q[1];
rz(-0.39127884) q[1];
sx q[1];
rz(1.7404187) q[1];
rz(-pi) q[2];
rz(-0.55283847) q[3];
sx q[3];
rz(-1.1854608) q[3];
sx q[3];
rz(2.5401153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99453551) q[2];
sx q[2];
rz(-0.68238443) q[2];
sx q[2];
rz(-0.4500173) q[2];
rz(-0.58722812) q[3];
sx q[3];
rz(-1.8000032) q[3];
sx q[3];
rz(-1.2500866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8477972) q[0];
sx q[0];
rz(-1.4601409) q[0];
sx q[0];
rz(-1.9548804) q[0];
rz(2.8222491) q[1];
sx q[1];
rz(-1.0217383) q[1];
sx q[1];
rz(-0.26225463) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866383) q[0];
sx q[0];
rz(-2.1146746) q[0];
sx q[0];
rz(1.8094814) q[0];
x q[1];
rz(-1.9840711) q[2];
sx q[2];
rz(-1.6954633) q[2];
sx q[2];
rz(1.1050129) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1848534) q[1];
sx q[1];
rz(-1.0295086) q[1];
sx q[1];
rz(-0.095603099) q[1];
x q[2];
rz(1.1434302) q[3];
sx q[3];
rz(-1.1550723) q[3];
sx q[3];
rz(0.24408578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.78476) q[2];
sx q[2];
rz(-2.6748952) q[2];
sx q[2];
rz(-0.32877767) q[2];
rz(-1.6262936) q[3];
sx q[3];
rz(-1.7833775) q[3];
sx q[3];
rz(-2.7330107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78557712) q[0];
sx q[0];
rz(-2.6308036) q[0];
sx q[0];
rz(-2.862634) q[0];
rz(0.51697671) q[1];
sx q[1];
rz(-0.37716436) q[1];
sx q[1];
rz(-1.7163716) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0128929) q[0];
sx q[0];
rz(-2.1919247) q[0];
sx q[0];
rz(1.5170044) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3383612) q[2];
sx q[2];
rz(-2.0798426) q[2];
sx q[2];
rz(2.5521297) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9659075) q[1];
sx q[1];
rz(-1.07465) q[1];
sx q[1];
rz(-1.7860852) q[1];
rz(-2.2278085) q[3];
sx q[3];
rz(-1.803064) q[3];
sx q[3];
rz(-0.454458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4031389) q[2];
sx q[2];
rz(-0.29610115) q[2];
sx q[2];
rz(-2.9816755) q[2];
rz(-2.3219409) q[3];
sx q[3];
rz(-1.22217) q[3];
sx q[3];
rz(-2.405449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6920456) q[0];
sx q[0];
rz(-1.291438) q[0];
sx q[0];
rz(-0.29577574) q[0];
rz(2.6462818) q[1];
sx q[1];
rz(-0.43423978) q[1];
sx q[1];
rz(-1.6326509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49182941) q[0];
sx q[0];
rz(-1.9252281) q[0];
sx q[0];
rz(1.9948122) q[0];
rz(-0.052148722) q[2];
sx q[2];
rz(-1.9226769) q[2];
sx q[2];
rz(-0.032648409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1965643) q[1];
sx q[1];
rz(-1.5032282) q[1];
sx q[1];
rz(-3.0516207) q[1];
rz(-pi) q[2];
rz(-1.9604574) q[3];
sx q[3];
rz(-1.1077322) q[3];
sx q[3];
rz(2.4135426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32790023) q[2];
sx q[2];
rz(-1.6715965) q[2];
sx q[2];
rz(2.5076765) q[2];
rz(0.085112326) q[3];
sx q[3];
rz(-1.2462933) q[3];
sx q[3];
rz(2.3688721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2278628) q[0];
sx q[0];
rz(-1.5999404) q[0];
sx q[0];
rz(3.0318442) q[0];
rz(-1.2211424) q[1];
sx q[1];
rz(-0.84538645) q[1];
sx q[1];
rz(0.56199817) q[1];
rz(2.0094677) q[2];
sx q[2];
rz(-1.7747468) q[2];
sx q[2];
rz(-2.940831) q[2];
rz(-1.3391277) q[3];
sx q[3];
rz(-1.5471519) q[3];
sx q[3];
rz(0.027761264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
