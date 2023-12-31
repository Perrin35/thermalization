OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(-0.68343502) q[0];
sx q[0];
rz(-2.662822) q[0];
rz(-3.1105644) q[1];
sx q[1];
rz(-1.9801158) q[1];
sx q[1];
rz(2.5002313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35533479) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(-1.3134365) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37748572) q[2];
sx q[2];
rz(-1.1698206) q[2];
sx q[2];
rz(0.12667835) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8909059) q[1];
sx q[1];
rz(-2.2911934) q[1];
sx q[1];
rz(-2.5198063) q[1];
rz(-pi) q[2];
rz(-1.9786505) q[3];
sx q[3];
rz(-2.1601094) q[3];
sx q[3];
rz(-2.4401963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.550094) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(-2.6634898) q[2];
rz(1.6889307) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96145445) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(-1.0189198) q[0];
rz(-1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(2.5085124) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7019254) q[0];
sx q[0];
rz(-2.3283726) q[0];
sx q[0];
rz(1.4600091) q[0];
rz(0.98302977) q[2];
sx q[2];
rz(-0.9126185) q[2];
sx q[2];
rz(2.7764729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9004904) q[1];
sx q[1];
rz(-1.2304658) q[1];
sx q[1];
rz(1.9963005) q[1];
x q[2];
rz(0.88926104) q[3];
sx q[3];
rz(-1.325843) q[3];
sx q[3];
rz(2.3219061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7818266) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(-0.11745545) q[2];
rz(-0.30101267) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(-1.3628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2837219) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(3.0531847) q[0];
rz(-1.1075426) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(-2.450768) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82949626) q[0];
sx q[0];
rz(-1.7236241) q[0];
sx q[0];
rz(1.6468847) q[0];
rz(-pi) q[1];
rz(2.3766741) q[2];
sx q[2];
rz(-1.7419445) q[2];
sx q[2];
rz(-2.314832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55331992) q[1];
sx q[1];
rz(-1.5968482) q[1];
sx q[1];
rz(-1.3887029) q[1];
x q[2];
rz(-1.9224123) q[3];
sx q[3];
rz(-1.2120486) q[3];
sx q[3];
rz(-0.84695942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2174125) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(0.10647354) q[2];
rz(-1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(-1.4829372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0687662) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(-2.6191214) q[0];
rz(0.32896313) q[1];
sx q[1];
rz(-1.4986228) q[1];
sx q[1];
rz(2.5879588) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41933665) q[0];
sx q[0];
rz(-2.0291078) q[0];
sx q[0];
rz(-2.3769747) q[0];
rz(-pi) q[1];
rz(-0.81981084) q[2];
sx q[2];
rz(-2.1852583) q[2];
sx q[2];
rz(-2.7053506) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8071027) q[1];
sx q[1];
rz(-2.346171) q[1];
sx q[1];
rz(2.9544178) q[1];
rz(-pi) q[2];
rz(-1.5418566) q[3];
sx q[3];
rz(-2.3674298) q[3];
sx q[3];
rz(2.418747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5840977) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(-0.49475691) q[2];
rz(-0.90302145) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(1.8516699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6506127) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(-1.0850798) q[0];
rz(2.1814573) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(2.9575612) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9404011) q[0];
sx q[0];
rz(-1.2815676) q[0];
sx q[0];
rz(1.679923) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8475624) q[2];
sx q[2];
rz(-1.9230611) q[2];
sx q[2];
rz(-2.2500452) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6726471) q[1];
sx q[1];
rz(-2.1429859) q[1];
sx q[1];
rz(-1.0020301) q[1];
rz(0.72599756) q[3];
sx q[3];
rz(-1.4193924) q[3];
sx q[3];
rz(-0.21522537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2400143) q[2];
sx q[2];
rz(-0.34873909) q[2];
sx q[2];
rz(-2.0007755) q[2];
rz(-0.42282894) q[3];
sx q[3];
rz(-1.4774277) q[3];
sx q[3];
rz(-2.0146577) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35247701) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(0.88678962) q[0];
rz(-0.24041644) q[1];
sx q[1];
rz(-1.0188894) q[1];
sx q[1];
rz(-0.14850798) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1214971) q[0];
sx q[0];
rz(-1.7317061) q[0];
sx q[0];
rz(-0.73385977) q[0];
rz(-pi) q[1];
rz(0.44273419) q[2];
sx q[2];
rz(-0.25207439) q[2];
sx q[2];
rz(2.4661951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6993466) q[1];
sx q[1];
rz(-1.3117426) q[1];
sx q[1];
rz(-0.11286347) q[1];
rz(1.1735736) q[3];
sx q[3];
rz(-1.1363315) q[3];
sx q[3];
rz(-1.43515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.285816) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(-2.6532069) q[2];
rz(2.6521902) q[3];
sx q[3];
rz(-2.2198052) q[3];
sx q[3];
rz(1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69328904) q[0];
sx q[0];
rz(-1.8087837) q[0];
sx q[0];
rz(2.3235902) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(1.12524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97911191) q[0];
sx q[0];
rz(-1.9444124) q[0];
sx q[0];
rz(-0.11066779) q[0];
rz(-pi) q[1];
rz(-1.869404) q[2];
sx q[2];
rz(-1.4223137) q[2];
sx q[2];
rz(0.4543002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4459671) q[1];
sx q[1];
rz(-1.8672767) q[1];
sx q[1];
rz(-2.8545024) q[1];
x q[2];
rz(2.1358228) q[3];
sx q[3];
rz(-2.1566448) q[3];
sx q[3];
rz(0.98423959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0685048) q[2];
sx q[2];
rz(-2.1606074) q[2];
sx q[2];
rz(-2.4970064) q[2];
rz(1.6623496) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(-1.7361599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1709764) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(1.9203141) q[1];
sx q[1];
rz(-1.5131283) q[1];
sx q[1];
rz(-2.1309526) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0518274) q[0];
sx q[0];
rz(-0.83501378) q[0];
sx q[0];
rz(0.21659539) q[0];
rz(-pi) q[1];
rz(2.0660731) q[2];
sx q[2];
rz(-1.3066548) q[2];
sx q[2];
rz(-1.8758945) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.752754) q[1];
sx q[1];
rz(-0.32143444) q[1];
sx q[1];
rz(-0.80498098) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40057064) q[3];
sx q[3];
rz(-1.4860324) q[3];
sx q[3];
rz(-0.29464196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6179787) q[2];
sx q[2];
rz(-0.31948677) q[2];
sx q[2];
rz(2.4475205) q[2];
rz(2.5726035) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(-0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086031832) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(-0.49474299) q[0];
rz(2.5231979) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(3.0659952) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48329096) q[0];
sx q[0];
rz(-1.0169944) q[0];
sx q[0];
rz(-0.36853405) q[0];
rz(-pi) q[1];
rz(-1.2443301) q[2];
sx q[2];
rz(-0.18507659) q[2];
sx q[2];
rz(-1.9290123) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52121431) q[1];
sx q[1];
rz(-1.673939) q[1];
sx q[1];
rz(-0.68876536) q[1];
rz(-2.5209386) q[3];
sx q[3];
rz(-1.7194347) q[3];
sx q[3];
rz(1.8295446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8081234) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(-1.3405651) q[2];
rz(0.30424413) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(-2.5568967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8463523) q[0];
sx q[0];
rz(-1.8363991) q[0];
sx q[0];
rz(0.17679581) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(0.44100824) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6498043) q[0];
sx q[0];
rz(-1.7735964) q[0];
sx q[0];
rz(-2.9359398) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44956019) q[2];
sx q[2];
rz(-2.9712147) q[2];
sx q[2];
rz(0.61194387) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9778206) q[1];
sx q[1];
rz(-1.2308916) q[1];
sx q[1];
rz(2.0274859) q[1];
rz(-pi) q[2];
rz(1.2714083) q[3];
sx q[3];
rz(-2.1736439) q[3];
sx q[3];
rz(2.833948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56197721) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(-0.30187541) q[2];
rz(0.89312303) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4176035) q[0];
sx q[0];
rz(-1.8483193) q[0];
sx q[0];
rz(1.666477) q[0];
rz(0.026731116) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(-1.0529636) q[2];
sx q[2];
rz(-0.79635194) q[2];
sx q[2];
rz(1.2780381) q[2];
rz(0.55587739) q[3];
sx q[3];
rz(-2.0049958) q[3];
sx q[3];
rz(-0.47331664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
