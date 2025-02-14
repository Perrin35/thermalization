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
rz(-1.4650605) q[0];
sx q[0];
rz(-0.2413916) q[0];
sx q[0];
rz(-3.0093332) q[0];
rz(-0.013068696) q[1];
sx q[1];
rz(-0.71813923) q[1];
sx q[1];
rz(-0.018996039) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.976124) q[0];
sx q[0];
rz(-0.69688334) q[0];
sx q[0];
rz(1.2319618) q[0];
rz(3.0471184) q[2];
sx q[2];
rz(-1.3850537) q[2];
sx q[2];
rz(2.8217725) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5031095) q[1];
sx q[1];
rz(-1.7448402) q[1];
sx q[1];
rz(-0.32764224) q[1];
x q[2];
rz(0.98145841) q[3];
sx q[3];
rz(-2.8638726) q[3];
sx q[3];
rz(-0.026175682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.12048177) q[2];
sx q[2];
rz(-0.80664539) q[2];
sx q[2];
rz(-2.3497904) q[2];
rz(-2.0558489) q[3];
sx q[3];
rz(-1.9224242) q[3];
sx q[3];
rz(-2.048548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2708112) q[0];
sx q[0];
rz(-2.9830611) q[0];
sx q[0];
rz(0.32919163) q[0];
rz(-1.5022494) q[1];
sx q[1];
rz(-0.90198016) q[1];
sx q[1];
rz(0.42218581) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80543431) q[0];
sx q[0];
rz(-2.1065188) q[0];
sx q[0];
rz(-2.0690919) q[0];
rz(-1.7549424) q[2];
sx q[2];
rz(-1.0304385) q[2];
sx q[2];
rz(-0.042405142) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3637963) q[1];
sx q[1];
rz(-1.8945122) q[1];
sx q[1];
rz(1.5219076) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77236891) q[3];
sx q[3];
rz(-1.4968781) q[3];
sx q[3];
rz(-1.0699917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72055703) q[2];
sx q[2];
rz(-1.9562419) q[2];
sx q[2];
rz(1.1400918) q[2];
rz(-0.0044048443) q[3];
sx q[3];
rz(-1.5604228) q[3];
sx q[3];
rz(0.13979039) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7797101) q[0];
sx q[0];
rz(-3.0411868) q[0];
sx q[0];
rz(0.34565872) q[0];
rz(2.0480305) q[1];
sx q[1];
rz(-0.82229096) q[1];
sx q[1];
rz(1.0543157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3234184) q[0];
sx q[0];
rz(-1.2369122) q[0];
sx q[0];
rz(2.7310577) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3286211) q[2];
sx q[2];
rz(-1.4958428) q[2];
sx q[2];
rz(-3.093442) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40151628) q[1];
sx q[1];
rz(-1.1187727) q[1];
sx q[1];
rz(1.6416618) q[1];
rz(-2.2369713) q[3];
sx q[3];
rz(-1.5525177) q[3];
sx q[3];
rz(1.3410904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4652319) q[2];
sx q[2];
rz(-0.54764843) q[2];
sx q[2];
rz(-2.5993627) q[2];
rz(-0.17255653) q[3];
sx q[3];
rz(-1.4901284) q[3];
sx q[3];
rz(2.0358613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3365823) q[0];
sx q[0];
rz(-1.7886826) q[0];
sx q[0];
rz(-1.2611058) q[0];
rz(-2.005596) q[1];
sx q[1];
rz(-2.0295862) q[1];
sx q[1];
rz(0.13066185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74742252) q[0];
sx q[0];
rz(-0.46499377) q[0];
sx q[0];
rz(0.93469145) q[0];
rz(-pi) q[1];
rz(-1.7950141) q[2];
sx q[2];
rz(-2.0665238) q[2];
sx q[2];
rz(-2.5350646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7893204) q[1];
sx q[1];
rz(-2.6323458) q[1];
sx q[1];
rz(1.5058084) q[1];
rz(2.6819508) q[3];
sx q[3];
rz(-2.423624) q[3];
sx q[3];
rz(-1.8945862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38264349) q[2];
sx q[2];
rz(-3.068277) q[2];
sx q[2];
rz(2.9200413) q[2];
rz(0.70212901) q[3];
sx q[3];
rz(-1.5626855) q[3];
sx q[3];
rz(2.4041972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5704983) q[0];
sx q[0];
rz(-1.2449188) q[0];
sx q[0];
rz(3.0076497) q[0];
rz(2.7722923) q[1];
sx q[1];
rz(-1.9422653) q[1];
sx q[1];
rz(1.3901002) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5126915) q[0];
sx q[0];
rz(-2.5675964) q[0];
sx q[0];
rz(1.1074793) q[0];
rz(2.8763387) q[2];
sx q[2];
rz(-1.2758288) q[2];
sx q[2];
rz(-2.2504928) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0005324) q[1];
sx q[1];
rz(-0.62126011) q[1];
sx q[1];
rz(0.48807524) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69306121) q[3];
sx q[3];
rz(-2.6821838) q[3];
sx q[3];
rz(1.0074248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2964581) q[2];
sx q[2];
rz(-0.90201169) q[2];
sx q[2];
rz(1.585539) q[2];
rz(0.29160094) q[3];
sx q[3];
rz(-2.6652938) q[3];
sx q[3];
rz(0.27879032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6911102) q[0];
sx q[0];
rz(-2.1463558) q[0];
sx q[0];
rz(-0.52571785) q[0];
rz(-0.37824962) q[1];
sx q[1];
rz(-1.6098166) q[1];
sx q[1];
rz(0.27413109) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7129242) q[0];
sx q[0];
rz(-2.1607421) q[0];
sx q[0];
rz(-1.601786) q[0];
rz(-1.1172764) q[2];
sx q[2];
rz(-1.9651061) q[2];
sx q[2];
rz(2.0298983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8291446) q[1];
sx q[1];
rz(-0.48407468) q[1];
sx q[1];
rz(2.1425254) q[1];
rz(-1.8242307) q[3];
sx q[3];
rz(-1.4259286) q[3];
sx q[3];
rz(-1.0235909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5963001) q[2];
sx q[2];
rz(-1.9175074) q[2];
sx q[2];
rz(0.25823414) q[2];
rz(-1.5457414) q[3];
sx q[3];
rz(-1.7786547) q[3];
sx q[3];
rz(2.7365007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79893583) q[0];
sx q[0];
rz(-0.5181784) q[0];
sx q[0];
rz(-3.0361191) q[0];
rz(-2.8666829) q[1];
sx q[1];
rz(-1.7173488) q[1];
sx q[1];
rz(2.8555433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8409922) q[0];
sx q[0];
rz(-1.4541885) q[0];
sx q[0];
rz(-0.18516115) q[0];
rz(-1.2484543) q[2];
sx q[2];
rz(-2.3791056) q[2];
sx q[2];
rz(0.74021268) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6374928) q[1];
sx q[1];
rz(-1.5063725) q[1];
sx q[1];
rz(-1.1845876) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4838832) q[3];
sx q[3];
rz(-0.66221395) q[3];
sx q[3];
rz(-2.7194104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1470571) q[2];
sx q[2];
rz(-2.4592082) q[2];
sx q[2];
rz(-0.4500173) q[2];
rz(2.5543645) q[3];
sx q[3];
rz(-1.3415895) q[3];
sx q[3];
rz(-1.8915061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8477972) q[0];
sx q[0];
rz(-1.4601409) q[0];
sx q[0];
rz(-1.9548804) q[0];
rz(-0.31934357) q[1];
sx q[1];
rz(-2.1198544) q[1];
sx q[1];
rz(0.26225463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5004999) q[0];
sx q[0];
rz(-1.3670792) q[0];
sx q[0];
rz(2.5848956) q[0];
rz(1.1575216) q[2];
sx q[2];
rz(-1.6954633) q[2];
sx q[2];
rz(1.1050129) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7062807) q[1];
sx q[1];
rz(-1.6526994) q[1];
sx q[1];
rz(1.0274853) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1434302) q[3];
sx q[3];
rz(-1.9865204) q[3];
sx q[3];
rz(-2.8975069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.78476) q[2];
sx q[2];
rz(-0.46669745) q[2];
sx q[2];
rz(-2.812815) q[2];
rz(-1.6262936) q[3];
sx q[3];
rz(-1.3582151) q[3];
sx q[3];
rz(-0.40858194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78557712) q[0];
sx q[0];
rz(-0.51078904) q[0];
sx q[0];
rz(-0.27895862) q[0];
rz(2.6246159) q[1];
sx q[1];
rz(-2.7644283) q[1];
sx q[1];
rz(-1.7163716) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681717) q[0];
sx q[0];
rz(-1.6145339) q[0];
sx q[0];
rz(2.5197791) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4819229) q[2];
sx q[2];
rz(-2.2222509) q[2];
sx q[2];
rz(1.4209117) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17568517) q[1];
sx q[1];
rz(-2.0669427) q[1];
sx q[1];
rz(-1.3555075) q[1];
rz(-2.2278085) q[3];
sx q[3];
rz(-1.3385286) q[3];
sx q[3];
rz(-2.6871347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4031389) q[2];
sx q[2];
rz(-0.29610115) q[2];
sx q[2];
rz(2.9816755) q[2];
rz(2.3219409) q[3];
sx q[3];
rz(-1.9194226) q[3];
sx q[3];
rz(0.73614365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6920456) q[0];
sx q[0];
rz(-1.291438) q[0];
sx q[0];
rz(-2.8458169) q[0];
rz(0.4953109) q[1];
sx q[1];
rz(-2.7073529) q[1];
sx q[1];
rz(-1.6326509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9072278) q[0];
sx q[0];
rz(-1.9669269) q[0];
sx q[0];
rz(0.38568003) q[0];
x q[1];
rz(-3.0894439) q[2];
sx q[2];
rz(-1.2189157) q[2];
sx q[2];
rz(3.1089442) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8733522) q[1];
sx q[1];
rz(-0.11246364) q[1];
sx q[1];
rz(-0.64551533) q[1];
rz(-2.4911777) q[3];
sx q[3];
rz(-2.5456508) q[3];
sx q[3];
rz(3.126248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32790023) q[2];
sx q[2];
rz(-1.4699961) q[2];
sx q[2];
rz(-0.63391614) q[2];
rz(-0.085112326) q[3];
sx q[3];
rz(-1.2462933) q[3];
sx q[3];
rz(0.77272052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91372981) q[0];
sx q[0];
rz(-1.5999404) q[0];
sx q[0];
rz(3.0318442) q[0];
rz(1.9204503) q[1];
sx q[1];
rz(-0.84538645) q[1];
sx q[1];
rz(0.56199817) q[1];
rz(-1.1176422) q[2];
sx q[2];
rz(-2.6606371) q[2];
sx q[2];
rz(1.3640192) q[2];
rz(-1.6734335) q[3];
sx q[3];
rz(-2.9087421) q[3];
sx q[3];
rz(-1.4431492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
