OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0786809) q[0];
sx q[0];
rz(-0.39830387) q[0];
sx q[0];
rz(2.8226573) q[0];
rz(-0.4660663) q[1];
sx q[1];
rz(-1.708344) q[1];
sx q[1];
rz(0.27369174) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39616991) q[0];
sx q[0];
rz(-1.92116) q[0];
sx q[0];
rz(0.78846447) q[0];
rz(-pi) q[1];
rz(0.75184043) q[2];
sx q[2];
rz(-0.32773563) q[2];
sx q[2];
rz(-2.73955) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0658354) q[1];
sx q[1];
rz(-0.21327107) q[1];
sx q[1];
rz(0.98827459) q[1];
rz(-2.3832604) q[3];
sx q[3];
rz(-1.8852295) q[3];
sx q[3];
rz(0.093859501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.29229257) q[2];
sx q[2];
rz(-2.2165074) q[2];
sx q[2];
rz(-1.7577457) q[2];
rz(0.84550953) q[3];
sx q[3];
rz(-3.129831) q[3];
sx q[3];
rz(0.99172926) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94075769) q[0];
sx q[0];
rz(-2.2469914) q[0];
sx q[0];
rz(1.0351329) q[0];
rz(1.9738522) q[1];
sx q[1];
rz(-2.9393241) q[1];
sx q[1];
rz(2.367173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37062936) q[0];
sx q[0];
rz(-1.9815195) q[0];
sx q[0];
rz(-1.0324917) q[0];
rz(-pi) q[1];
rz(2.0032194) q[2];
sx q[2];
rz(-0.70821643) q[2];
sx q[2];
rz(-2.0067635) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.32774156) q[1];
sx q[1];
rz(-2.6151513) q[1];
sx q[1];
rz(1.3415643) q[1];
rz(-0.033942698) q[3];
sx q[3];
rz(-1.5205624) q[3];
sx q[3];
rz(1.7173345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8059798) q[2];
sx q[2];
rz(-3.0642982) q[2];
sx q[2];
rz(0.94266164) q[2];
rz(-3.0521657) q[3];
sx q[3];
rz(-0.663203) q[3];
sx q[3];
rz(-2.8360143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.5494004) q[0];
sx q[0];
rz(-2.5553199) q[0];
sx q[0];
rz(0.19202448) q[0];
rz(2.5281455) q[1];
sx q[1];
rz(-0.44812056) q[1];
sx q[1];
rz(3.1268069) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0917182) q[0];
sx q[0];
rz(-0.8862859) q[0];
sx q[0];
rz(0.52264638) q[0];
rz(-pi) q[1];
rz(2.4528265) q[2];
sx q[2];
rz(-1.3564132) q[2];
sx q[2];
rz(2.2685879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5303684) q[1];
sx q[1];
rz(-2.1448128) q[1];
sx q[1];
rz(2.8635129) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1802894) q[3];
sx q[3];
rz(-0.064924463) q[3];
sx q[3];
rz(-2.8466259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52745596) q[2];
sx q[2];
rz(-1.5946031) q[2];
sx q[2];
rz(-0.0097489348) q[2];
rz(2.5629432) q[3];
sx q[3];
rz(-1.0520244) q[3];
sx q[3];
rz(-2.3600522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(3.0233199) q[0];
sx q[0];
rz(-2.3319722) q[0];
sx q[0];
rz(-0.6672346) q[0];
rz(2.1369333) q[1];
sx q[1];
rz(-2.9861082) q[1];
sx q[1];
rz(-1.7612693) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28128703) q[0];
sx q[0];
rz(-1.5800887) q[0];
sx q[0];
rz(-1.0929266) q[0];
rz(-0.3227206) q[2];
sx q[2];
rz(-1.1185794) q[2];
sx q[2];
rz(-1.9911901) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9050153) q[1];
sx q[1];
rz(-1.6029753) q[1];
sx q[1];
rz(-2.6479812) q[1];
rz(-1.0466688) q[3];
sx q[3];
rz(-1.7780684) q[3];
sx q[3];
rz(2.2818499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37463793) q[2];
sx q[2];
rz(-1.1344323) q[2];
sx q[2];
rz(-3.0881506) q[2];
rz(-2.5080569) q[3];
sx q[3];
rz(-0.41826785) q[3];
sx q[3];
rz(0.0088648349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87090129) q[0];
sx q[0];
rz(-2.5636261) q[0];
sx q[0];
rz(-2.0722678) q[0];
rz(-1.2879734) q[1];
sx q[1];
rz(-3.0777212) q[1];
sx q[1];
rz(0.091648253) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6955513) q[0];
sx q[0];
rz(-1.5484838) q[0];
sx q[0];
rz(3.1299921) q[0];
rz(-pi) q[1];
rz(0.32817082) q[2];
sx q[2];
rz(-2.0053021) q[2];
sx q[2];
rz(-0.19319867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5875467) q[1];
sx q[1];
rz(-2.5623706) q[1];
sx q[1];
rz(-2.278797) q[1];
rz(1.5318457) q[3];
sx q[3];
rz(-2.340303) q[3];
sx q[3];
rz(-1.0391446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.84534591) q[2];
sx q[2];
rz(-2.3643934) q[2];
sx q[2];
rz(0.11543154) q[2];
rz(-0.7655862) q[3];
sx q[3];
rz(-1.4976394) q[3];
sx q[3];
rz(2.7287741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63824832) q[0];
sx q[0];
rz(-0.64953506) q[0];
sx q[0];
rz(0.58393884) q[0];
rz(3.0931547) q[1];
sx q[1];
rz(-2.9190639) q[1];
sx q[1];
rz(-2.4979874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2023425) q[0];
sx q[0];
rz(-1.5233375) q[0];
sx q[0];
rz(1.2623293) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2320643) q[2];
sx q[2];
rz(-1.9478746) q[2];
sx q[2];
rz(-0.23185767) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4847624) q[1];
sx q[1];
rz(-1.7028013) q[1];
sx q[1];
rz(2.2836192) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77966431) q[3];
sx q[3];
rz(-1.4309959) q[3];
sx q[3];
rz(-1.9817022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.863997) q[2];
sx q[2];
rz(-2.3433351) q[2];
sx q[2];
rz(-2.4994728) q[2];
rz(0.13907214) q[3];
sx q[3];
rz(-0.13834794) q[3];
sx q[3];
rz(1.6819491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794401) q[0];
sx q[0];
rz(-0.42069778) q[0];
sx q[0];
rz(0.62533373) q[0];
rz(-2.9026237) q[1];
sx q[1];
rz(-2.9049554) q[1];
sx q[1];
rz(1.3561358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0129798) q[0];
sx q[0];
rz(-1.5389568) q[0];
sx q[0];
rz(1.6252993) q[0];
rz(1.8618334) q[2];
sx q[2];
rz(-1.9328397) q[2];
sx q[2];
rz(-0.14823981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6125638) q[1];
sx q[1];
rz(-2.7269438) q[1];
sx q[1];
rz(-2.0842444) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5111968) q[3];
sx q[3];
rz(-0.26132628) q[3];
sx q[3];
rz(-3.0898409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8465392) q[2];
sx q[2];
rz(-1.2650547) q[2];
sx q[2];
rz(-0.96578252) q[2];
rz(-0.24671181) q[3];
sx q[3];
rz(-1.8762981) q[3];
sx q[3];
rz(1.3707967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63916373) q[0];
sx q[0];
rz(-0.40488365) q[0];
sx q[0];
rz(3.0054481) q[0];
rz(-2.4038521) q[1];
sx q[1];
rz(-2.9164011) q[1];
sx q[1];
rz(-1.9053316) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5571283) q[0];
sx q[0];
rz(-2.4274913) q[0];
sx q[0];
rz(-1.0034998) q[0];
rz(-pi) q[1];
rz(-0.18119104) q[2];
sx q[2];
rz(-2.199442) q[2];
sx q[2];
rz(-1.218007) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.511179) q[1];
sx q[1];
rz(-1.8466338) q[1];
sx q[1];
rz(2.7008369) q[1];
x q[2];
rz(0.32847877) q[3];
sx q[3];
rz(-0.94624472) q[3];
sx q[3];
rz(0.49027944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4514734) q[2];
sx q[2];
rz(-1.5855007) q[2];
sx q[2];
rz(-0.95009178) q[2];
rz(2.7443366) q[3];
sx q[3];
rz(-2.6451431) q[3];
sx q[3];
rz(0.44601405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5815247) q[0];
sx q[0];
rz(-0.6811322) q[0];
sx q[0];
rz(-0.11823046) q[0];
rz(-1.0826348) q[1];
sx q[1];
rz(-1.2018459) q[1];
sx q[1];
rz(1.3523678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0948321) q[0];
sx q[0];
rz(-2.7343395) q[0];
sx q[0];
rz(-2.880275) q[0];
rz(-pi) q[1];
rz(0.1583516) q[2];
sx q[2];
rz(-0.54979804) q[2];
sx q[2];
rz(-2.9937772) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4115826) q[1];
sx q[1];
rz(-1.4435539) q[1];
sx q[1];
rz(-2.0557269) q[1];
x q[2];
rz(2.8297573) q[3];
sx q[3];
rz(-1.1118299) q[3];
sx q[3];
rz(2.3941657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3236986) q[2];
sx q[2];
rz(-2.4511621) q[2];
sx q[2];
rz(0.77198088) q[2];
rz(3.0585994) q[3];
sx q[3];
rz(-1.8738184) q[3];
sx q[3];
rz(2.6193589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25987396) q[0];
sx q[0];
rz(-2.3840955) q[0];
sx q[0];
rz(0.74257332) q[0];
rz(-1.2965797) q[1];
sx q[1];
rz(-0.80755889) q[1];
sx q[1];
rz(0.47320941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98223393) q[0];
sx q[0];
rz(-1.8372941) q[0];
sx q[0];
rz(-3.052211) q[0];
rz(-pi) q[1];
rz(2.8933018) q[2];
sx q[2];
rz(-1.3220698) q[2];
sx q[2];
rz(1.7493713) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4025583) q[1];
sx q[1];
rz(-2.0328882) q[1];
sx q[1];
rz(-2.7527304) q[1];
x q[2];
rz(0.27373154) q[3];
sx q[3];
rz(-0.76430862) q[3];
sx q[3];
rz(2.4814062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5903198) q[2];
sx q[2];
rz(-1.4164305) q[2];
sx q[2];
rz(1.9080706) q[2];
rz(0.34520712) q[3];
sx q[3];
rz(-2.749741) q[3];
sx q[3];
rz(0.83693081) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.681916) q[0];
sx q[0];
rz(-1.3688594) q[0];
sx q[0];
rz(-1.3874227) q[0];
rz(3.0826898) q[1];
sx q[1];
rz(-1.420493) q[1];
sx q[1];
rz(-2.4891985) q[1];
rz(-2.7564213) q[2];
sx q[2];
rz(-1.9829911) q[2];
sx q[2];
rz(0.34172716) q[2];
rz(0.98109365) q[3];
sx q[3];
rz(-2.4193939) q[3];
sx q[3];
rz(-1.9388225) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
