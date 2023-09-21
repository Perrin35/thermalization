OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61385566) q[0];
sx q[0];
rz(4.6392415) q[0];
sx q[0];
rz(8.5949329) q[0];
rz(-2.3614376) q[1];
sx q[1];
rz(-1.0649788) q[1];
sx q[1];
rz(0.87632626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60183817) q[0];
sx q[0];
rz(-1.1798501) q[0];
sx q[0];
rz(1.8200309) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.949602) q[2];
sx q[2];
rz(-0.99388323) q[2];
sx q[2];
rz(2.752395) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4727002) q[1];
sx q[1];
rz(-1.2408222) q[1];
sx q[1];
rz(3.1256691) q[1];
rz(-pi) q[2];
x q[2];
rz(1.417744) q[3];
sx q[3];
rz(-1.9733841) q[3];
sx q[3];
rz(-0.11868417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(-0.7286287) q[2];
rz(-0.5209926) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(-0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.306863) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-2.0200502) q[0];
rz(2.8858378) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(2.2671525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16337285) q[0];
sx q[0];
rz(-1.1436497) q[0];
sx q[0];
rz(2.8073729) q[0];
rz(-pi) q[1];
rz(-1.3358243) q[2];
sx q[2];
rz(-0.58832303) q[2];
sx q[2];
rz(-0.36662835) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16972152) q[1];
sx q[1];
rz(-1.1217146) q[1];
sx q[1];
rz(2.4596618) q[1];
x q[2];
rz(3.087895) q[3];
sx q[3];
rz(-1.3018381) q[3];
sx q[3];
rz(-0.69124903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(-0.43593105) q[2];
rz(2.46051) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(-2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23713672) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(1.0748192) q[0];
rz(-2.3020321) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(-0.39594617) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0229189) q[0];
sx q[0];
rz(-2.4791105) q[0];
sx q[0];
rz(-0.89612095) q[0];
rz(-pi) q[1];
rz(-1.5244353) q[2];
sx q[2];
rz(-2.8472387) q[2];
sx q[2];
rz(0.72325318) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46519854) q[1];
sx q[1];
rz(-0.68323831) q[1];
sx q[1];
rz(0.015124358) q[1];
rz(-pi) q[2];
rz(-2.8860502) q[3];
sx q[3];
rz(-1.5934172) q[3];
sx q[3];
rz(0.25657755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6039156) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(2.9023857) q[2];
rz(-3.0662597) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.4199715) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(-2.3216632) q[0];
rz(0.48768249) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(-2.908169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.061302) q[0];
sx q[0];
rz(-1.6675321) q[0];
sx q[0];
rz(3.0808099) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64220631) q[2];
sx q[2];
rz(-1.3807447) q[2];
sx q[2];
rz(-1.1914636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0049131752) q[1];
sx q[1];
rz(-0.79489743) q[1];
sx q[1];
rz(-0.48308785) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78419533) q[3];
sx q[3];
rz(-1.2894221) q[3];
sx q[3];
rz(-2.1649862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0908115) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(0.74742571) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(2.8994765) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500047) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(-3.0773556) q[0];
rz(2.1977987) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(-0.76104004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2634537) q[0];
sx q[0];
rz(-1.8315151) q[0];
sx q[0];
rz(3.1104452) q[0];
rz(-pi) q[1];
rz(2.9391187) q[2];
sx q[2];
rz(-1.379181) q[2];
sx q[2];
rz(-1.7970049) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3984822) q[1];
sx q[1];
rz(-2.8201436) q[1];
sx q[1];
rz(-0.076996315) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8499591) q[3];
sx q[3];
rz(-1.2984635) q[3];
sx q[3];
rz(-0.84701049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3349907) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(3.0495194) q[2];
rz(-2.4798685) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(-1.7592643) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6732366) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-0.026542149) q[0];
rz(2.2684855) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(0.32593265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8317141) q[0];
sx q[0];
rz(-1.1967812) q[0];
sx q[0];
rz(-0.563234) q[0];
x q[1];
rz(-1.2665777) q[2];
sx q[2];
rz(-1.1433257) q[2];
sx q[2];
rz(1.3104591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95820108) q[1];
sx q[1];
rz(-0.32918731) q[1];
sx q[1];
rz(2.7126461) q[1];
x q[2];
rz(1.5876706) q[3];
sx q[3];
rz(-1.4697187) q[3];
sx q[3];
rz(2.8942787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59763336) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(1.4298965) q[3];
sx q[3];
rz(-1.1681898) q[3];
sx q[3];
rz(-0.54106075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-0.72189271) q[0];
rz(-1.7294653) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(0.11925764) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3463466) q[0];
sx q[0];
rz(-2.3756785) q[0];
sx q[0];
rz(2.7555097) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7330405) q[2];
sx q[2];
rz(-1.4622697) q[2];
sx q[2];
rz(-2.721399) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5057999) q[1];
sx q[1];
rz(-2.0675106) q[1];
sx q[1];
rz(-2.990681) q[1];
rz(-pi) q[2];
rz(0.42640949) q[3];
sx q[3];
rz(-0.75424131) q[3];
sx q[3];
rz(-2.756556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(-1.7822441) q[2];
rz(2.3826777) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(-2.6045077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83157241) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(1.0634364) q[0];
rz(-2.8670782) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(-0.88561052) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6016156) q[0];
sx q[0];
rz(-2.5936539) q[0];
sx q[0];
rz(-2.6849296) q[0];
rz(-pi) q[1];
x q[1];
rz(0.015651264) q[2];
sx q[2];
rz(-0.99132292) q[2];
sx q[2];
rz(-1.5660812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9319716) q[1];
sx q[1];
rz(-1.069427) q[1];
sx q[1];
rz(2.6905836) q[1];
rz(0.86730154) q[3];
sx q[3];
rz(-2.2059545) q[3];
sx q[3];
rz(2.9150073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(2.0098861) q[2];
rz(-1.0845832) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(-1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30329147) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(-2.0595179) q[0];
rz(1.2754296) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(1.1358322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1866859) q[0];
sx q[0];
rz(-2.1011155) q[0];
sx q[0];
rz(-2.9493939) q[0];
rz(-1.5485974) q[2];
sx q[2];
rz(-2.0054842) q[2];
sx q[2];
rz(-1.7322025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.218704) q[1];
sx q[1];
rz(-1.6195546) q[1];
sx q[1];
rz(-0.4263652) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7157852) q[3];
sx q[3];
rz(-0.96447456) q[3];
sx q[3];
rz(3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0231126) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(-0.87289587) q[2];
rz(0.84351271) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89235598) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(-1.0247914) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(-1.1970253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7322757) q[0];
sx q[0];
rz(-1.1336375) q[0];
sx q[0];
rz(0.93426312) q[0];
rz(-pi) q[1];
rz(-2.0335774) q[2];
sx q[2];
rz(-0.26460755) q[2];
sx q[2];
rz(0.77361425) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7750818) q[1];
sx q[1];
rz(-1.0771891) q[1];
sx q[1];
rz(2.2713186) q[1];
rz(-pi) q[2];
rz(-0.61693807) q[3];
sx q[3];
rz(-2.5411798) q[3];
sx q[3];
rz(-2.2850349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3060351) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(0.79997921) q[2];
rz(1.9647313) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(-2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70893127) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(-2.6196383) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(-2.7675046) q[2];
sx q[2];
rz(-1.903152) q[2];
sx q[2];
rz(-2.8852035) q[2];
rz(2.2189191) q[3];
sx q[3];
rz(-1.0395944) q[3];
sx q[3];
rz(1.1208054) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
