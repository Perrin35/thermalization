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
rz(1.849527) q[0];
sx q[0];
rz(5.4592291) q[0];
sx q[0];
rz(8.4973314) q[0];
rz(-1.0499586) q[1];
sx q[1];
rz(-2.1702622) q[1];
sx q[1];
rz(-2.3271022) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2789291) q[0];
sx q[0];
rz(-1.188108) q[0];
sx q[0];
rz(-2.0576059) q[0];
rz(-pi) q[1];
rz(0.94169803) q[2];
sx q[2];
rz(-0.88028958) q[2];
sx q[2];
rz(2.5494573) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3315288) q[1];
sx q[1];
rz(-2.0671856) q[1];
sx q[1];
rz(2.1138825) q[1];
rz(-pi) q[2];
rz(2.2839969) q[3];
sx q[3];
rz(-0.88123816) q[3];
sx q[3];
rz(-1.621472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77259511) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(-2.7102846) q[2];
rz(1.4946233) q[3];
sx q[3];
rz(-1.5458115) q[3];
sx q[3];
rz(-2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876708) q[0];
sx q[0];
rz(-2.9449154) q[0];
sx q[0];
rz(-1.3171296) q[0];
rz(-2.092284) q[1];
sx q[1];
rz(-0.53593719) q[1];
sx q[1];
rz(-0.61872331) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71402446) q[0];
sx q[0];
rz(-2.51665) q[0];
sx q[0];
rz(-2.8929536) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0486062) q[2];
sx q[2];
rz(-2.3903762) q[2];
sx q[2];
rz(-0.71556811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6438457) q[1];
sx q[1];
rz(-0.37393716) q[1];
sx q[1];
rz(-1.9356329) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57024184) q[3];
sx q[3];
rz(-2.5240353) q[3];
sx q[3];
rz(-2.9887426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72545663) q[2];
sx q[2];
rz(-1.340103) q[2];
sx q[2];
rz(-0.49187342) q[2];
rz(-1.5791996) q[3];
sx q[3];
rz(-1.4926566) q[3];
sx q[3];
rz(-0.57945848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014515011) q[0];
sx q[0];
rz(-1.0579695) q[0];
sx q[0];
rz(-2.8614817) q[0];
rz(-3.0666871) q[1];
sx q[1];
rz(-0.43498755) q[1];
sx q[1];
rz(1.3704971) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89557392) q[0];
sx q[0];
rz(-1.7798806) q[0];
sx q[0];
rz(0.81935291) q[0];
rz(-2.2417775) q[2];
sx q[2];
rz(-2.0986989) q[2];
sx q[2];
rz(-1.4769276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0468415) q[1];
sx q[1];
rz(-2.1194723) q[1];
sx q[1];
rz(1.2612753) q[1];
x q[2];
rz(1.5519484) q[3];
sx q[3];
rz(-1.7030088) q[3];
sx q[3];
rz(2.978505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.28477731) q[2];
sx q[2];
rz(-1.4692401) q[2];
sx q[2];
rz(-2.3397297) q[2];
rz(0.91184584) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(2.1626933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.3367679) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(2.5208852) q[0];
rz(0.2785109) q[1];
sx q[1];
rz(-1.4418437) q[1];
sx q[1];
rz(-1.0034358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0517563) q[0];
sx q[0];
rz(-0.47430719) q[0];
sx q[0];
rz(-0.82292212) q[0];
x q[1];
rz(3.0714967) q[2];
sx q[2];
rz(-1.150944) q[2];
sx q[2];
rz(1.9934838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.089503376) q[1];
sx q[1];
rz(-1.6253396) q[1];
sx q[1];
rz(0.2691582) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8068321) q[3];
sx q[3];
rz(-0.98525648) q[3];
sx q[3];
rz(-1.745831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8371381) q[2];
sx q[2];
rz(-2.0292323) q[2];
sx q[2];
rz(-3.0755074) q[2];
rz(-1.1262013) q[3];
sx q[3];
rz(-2.3812713) q[3];
sx q[3];
rz(-2.5662305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5758301) q[0];
sx q[0];
rz(-0.081871651) q[0];
sx q[0];
rz(2.3542812) q[0];
rz(0.85834223) q[1];
sx q[1];
rz(-0.97802496) q[1];
sx q[1];
rz(-1.0816921) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3560548) q[0];
sx q[0];
rz(-0.83844664) q[0];
sx q[0];
rz(-2.4466956) q[0];
rz(2.2901847) q[2];
sx q[2];
rz(-0.59307304) q[2];
sx q[2];
rz(2.3874758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.2277775) q[1];
sx q[1];
rz(-2.8568513) q[1];
sx q[1];
rz(0.022241994) q[1];
rz(-2.5852935) q[3];
sx q[3];
rz(-2.2101758) q[3];
sx q[3];
rz(-2.104708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.16296884) q[2];
sx q[2];
rz(-1.9917515) q[2];
sx q[2];
rz(-2.2990189) q[2];
rz(-2.8652371) q[3];
sx q[3];
rz(-0.98139757) q[3];
sx q[3];
rz(1.3493376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.855298) q[0];
sx q[0];
rz(-1.3110302) q[0];
sx q[0];
rz(1.1899765) q[0];
rz(-1.9189934) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(0.58293265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8609498) q[0];
sx q[0];
rz(-1.1034729) q[0];
sx q[0];
rz(-1.1790465) q[0];
rz(-pi) q[1];
rz(2.1443679) q[2];
sx q[2];
rz(-1.5825446) q[2];
sx q[2];
rz(0.16596101) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5546023) q[1];
sx q[1];
rz(-1.490056) q[1];
sx q[1];
rz(2.9029773) q[1];
rz(-2.0058026) q[3];
sx q[3];
rz(-2.1564751) q[3];
sx q[3];
rz(-3.1068561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0588093) q[2];
sx q[2];
rz(-1.8913816) q[2];
sx q[2];
rz(2.1018551) q[2];
rz(-0.58935634) q[3];
sx q[3];
rz(-1.6683234) q[3];
sx q[3];
rz(1.0065494) q[3];
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
rz(0.37650087) q[0];
sx q[0];
rz(-0.14987513) q[0];
sx q[0];
rz(-1.804922) q[0];
rz(-0.22557766) q[1];
sx q[1];
rz(-1.0712737) q[1];
sx q[1];
rz(0.76990661) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1223567) q[0];
sx q[0];
rz(-1.4154288) q[0];
sx q[0];
rz(0.47433406) q[0];
x q[1];
rz(-0.90183423) q[2];
sx q[2];
rz(-2.2775501) q[2];
sx q[2];
rz(1.2400986) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6698031) q[1];
sx q[1];
rz(-2.9297929) q[1];
sx q[1];
rz(2.7752911) q[1];
rz(2.5363471) q[3];
sx q[3];
rz(-1.5063933) q[3];
sx q[3];
rz(1.6689672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.3930438) q[2];
sx q[2];
rz(-2.1262271) q[2];
rz(0.20396248) q[3];
sx q[3];
rz(-1.5488307) q[3];
sx q[3];
rz(1.7502194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2400951) q[0];
sx q[0];
rz(-0.63145852) q[0];
sx q[0];
rz(-0.39328662) q[0];
rz(2.6914864) q[1];
sx q[1];
rz(-1.422912) q[1];
sx q[1];
rz(0.030287655) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05370985) q[0];
sx q[0];
rz(-1.5397414) q[0];
sx q[0];
rz(2.8939918) q[0];
rz(-pi) q[1];
rz(-2.3322163) q[2];
sx q[2];
rz(-1.7580108) q[2];
sx q[2];
rz(-0.16448122) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.92851621) q[1];
sx q[1];
rz(-1.5476481) q[1];
sx q[1];
rz(-0.54575463) q[1];
x q[2];
rz(-2.6859517) q[3];
sx q[3];
rz(-0.29248387) q[3];
sx q[3];
rz(2.4258421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0344737) q[2];
sx q[2];
rz(-1.5566885) q[2];
sx q[2];
rz(-1.1220453) q[2];
rz(1.0551039) q[3];
sx q[3];
rz(-0.3873581) q[3];
sx q[3];
rz(-1.8900227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9686862) q[0];
sx q[0];
rz(-0.98772573) q[0];
sx q[0];
rz(-2.3671142) q[0];
rz(-0.6340181) q[1];
sx q[1];
rz(-2.1114383) q[1];
sx q[1];
rz(2.074362) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7726583) q[0];
sx q[0];
rz(-2.308929) q[0];
sx q[0];
rz(1.5940109) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.048914476) q[2];
sx q[2];
rz(-1.8777784) q[2];
sx q[2];
rz(-0.129667) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53505234) q[1];
sx q[1];
rz(-2.7379704) q[1];
sx q[1];
rz(0.89151793) q[1];
x q[2];
rz(2.0499633) q[3];
sx q[3];
rz(-1.7288343) q[3];
sx q[3];
rz(-0.8534067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5344703) q[2];
sx q[2];
rz(-3.0202713) q[2];
sx q[2];
rz(-1.8758476) q[2];
rz(-3.1014465) q[3];
sx q[3];
rz(-2.7669192) q[3];
sx q[3];
rz(3.1025187) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6604615) q[0];
sx q[0];
rz(-0.88812319) q[0];
sx q[0];
rz(1.8817425) q[0];
rz(1.4866359) q[1];
sx q[1];
rz(-1.0340374) q[1];
sx q[1];
rz(-3.1390417) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0028062971) q[0];
sx q[0];
rz(-1.696179) q[0];
sx q[0];
rz(-0.12055293) q[0];
rz(-pi) q[1];
rz(-0.0095944842) q[2];
sx q[2];
rz(-1.6348607) q[2];
sx q[2];
rz(-1.4256776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3348744) q[1];
sx q[1];
rz(-0.93466991) q[1];
sx q[1];
rz(-1.0735372) q[1];
rz(-2.0052913) q[3];
sx q[3];
rz(-1.5948203) q[3];
sx q[3];
rz(-2.0498118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.5551787) q[2];
sx q[2];
rz(-1.5191583) q[2];
sx q[2];
rz(0.67443887) q[2];
rz(-2.0241578) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(0.31291541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-3.0565223) q[0];
sx q[0];
rz(-2.0105579) q[0];
sx q[0];
rz(2.6149268) q[0];
rz(-2.7936735) q[1];
sx q[1];
rz(-1.9031453) q[1];
sx q[1];
rz(-1.3487945) q[1];
rz(-1.3580218) q[2];
sx q[2];
rz(-1.909303) q[2];
sx q[2];
rz(2.252966) q[2];
rz(0.017853768) q[3];
sx q[3];
rz(-1.717106) q[3];
sx q[3];
rz(1.794869) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
