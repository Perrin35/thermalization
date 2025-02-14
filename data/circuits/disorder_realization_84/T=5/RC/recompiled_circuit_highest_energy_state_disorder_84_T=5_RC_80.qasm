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
rz(-0.33997384) q[0];
sx q[0];
rz(-0.61859328) q[0];
sx q[0];
rz(-3.0292188) q[0];
rz(-2.1984341) q[1];
sx q[1];
rz(-1.4961493) q[1];
sx q[1];
rz(-1.563031) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9364632) q[0];
sx q[0];
rz(-2.056707) q[0];
sx q[0];
rz(-1.0690752) q[0];
x q[1];
rz(-1.4077827) q[2];
sx q[2];
rz(-2.2634677) q[2];
sx q[2];
rz(2.7977347) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1670363) q[1];
sx q[1];
rz(-0.99460973) q[1];
sx q[1];
rz(-0.78414708) q[1];
rz(1.1279757) q[3];
sx q[3];
rz(-0.72924727) q[3];
sx q[3];
rz(-3.0506718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9872687) q[2];
sx q[2];
rz(-1.4877886) q[2];
sx q[2];
rz(-0.38779116) q[2];
rz(0.44529861) q[3];
sx q[3];
rz(-0.53578353) q[3];
sx q[3];
rz(2.9545412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8086789) q[0];
sx q[0];
rz(-2.6412922) q[0];
sx q[0];
rz(-2.152541) q[0];
rz(-1.5679081) q[1];
sx q[1];
rz(-0.75850073) q[1];
sx q[1];
rz(3.0737976) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7967927) q[0];
sx q[0];
rz(-1.4541309) q[0];
sx q[0];
rz(-0.33504265) q[0];
rz(-pi) q[1];
rz(1.2607093) q[2];
sx q[2];
rz(-1.3011622) q[2];
sx q[2];
rz(-0.53169429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2941563) q[1];
sx q[1];
rz(-1.7136317) q[1];
sx q[1];
rz(0.60522176) q[1];
rz(-pi) q[2];
rz(1.7420473) q[3];
sx q[3];
rz(-0.6500385) q[3];
sx q[3];
rz(2.5105944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13061358) q[2];
sx q[2];
rz(-2.4108672) q[2];
sx q[2];
rz(-2.0558426) q[2];
rz(-2.7340414) q[3];
sx q[3];
rz(-1.6435888) q[3];
sx q[3];
rz(1.7727859) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2456335) q[0];
sx q[0];
rz(-2.2838554) q[0];
sx q[0];
rz(-1.6813543) q[0];
rz(-1.1145837) q[1];
sx q[1];
rz(-2.3268301) q[1];
sx q[1];
rz(0.97273716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25165194) q[0];
sx q[0];
rz(-0.06315843) q[0];
sx q[0];
rz(-0.09357293) q[0];
rz(-pi) q[1];
rz(-0.26370987) q[2];
sx q[2];
rz(-2.5463841) q[2];
sx q[2];
rz(3.120228) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3537843) q[1];
sx q[1];
rz(-2.5671509) q[1];
sx q[1];
rz(-0.30409388) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17334799) q[3];
sx q[3];
rz(-0.9066994) q[3];
sx q[3];
rz(-0.085395902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14791791) q[2];
sx q[2];
rz(-2.6213054) q[2];
sx q[2];
rz(-0.49368039) q[2];
rz(-2.4364021) q[3];
sx q[3];
rz(-0.54809904) q[3];
sx q[3];
rz(1.7459474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.482835) q[0];
sx q[0];
rz(-0.63084298) q[0];
sx q[0];
rz(-3.0870497) q[0];
rz(1.4788117) q[1];
sx q[1];
rz(-2.3334267) q[1];
sx q[1];
rz(0.85173839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211713) q[0];
sx q[0];
rz(-1.2454438) q[0];
sx q[0];
rz(-2.6634548) q[0];
x q[1];
rz(-0.24864218) q[2];
sx q[2];
rz(-0.40462769) q[2];
sx q[2];
rz(2.8554631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6229565) q[1];
sx q[1];
rz(-1.1723435) q[1];
sx q[1];
rz(-1.4116628) q[1];
rz(-pi) q[2];
rz(-0.69290956) q[3];
sx q[3];
rz(-1.1947144) q[3];
sx q[3];
rz(-2.2364837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8405174) q[2];
sx q[2];
rz(-3.0121351) q[2];
sx q[2];
rz(-1.2811309) q[2];
rz(1.2766131) q[3];
sx q[3];
rz(-1.8254779) q[3];
sx q[3];
rz(-2.2006456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5285444) q[0];
sx q[0];
rz(-3.0597882) q[0];
sx q[0];
rz(1.4962037) q[0];
rz(1.243535) q[1];
sx q[1];
rz(-1.4356109) q[1];
sx q[1];
rz(-3.1028455) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.782283) q[0];
sx q[0];
rz(-1.4371431) q[0];
sx q[0];
rz(1.8275285) q[0];
rz(0.51182039) q[2];
sx q[2];
rz(-1.9241779) q[2];
sx q[2];
rz(1.3450373) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.52415614) q[1];
sx q[1];
rz(-1.8548522) q[1];
sx q[1];
rz(-1.3139105) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0465602) q[3];
sx q[3];
rz(-1.8390313) q[3];
sx q[3];
rz(-0.07650693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99899989) q[2];
sx q[2];
rz(-0.6730538) q[2];
sx q[2];
rz(-2.1633945) q[2];
rz(3.0603958) q[3];
sx q[3];
rz(-1.8742163) q[3];
sx q[3];
rz(0.80011884) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4303495) q[0];
sx q[0];
rz(-2.5177903) q[0];
sx q[0];
rz(1.8728363) q[0];
rz(0.87780344) q[1];
sx q[1];
rz(-1.9453847) q[1];
sx q[1];
rz(0.21003221) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97392573) q[0];
sx q[0];
rz(-0.96441482) q[0];
sx q[0];
rz(0.5180562) q[0];
rz(2.4792489) q[2];
sx q[2];
rz(-0.46842945) q[2];
sx q[2];
rz(-1.989502) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91466838) q[1];
sx q[1];
rz(-1.4836425) q[1];
sx q[1];
rz(2.8512849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29573567) q[3];
sx q[3];
rz(-1.9882759) q[3];
sx q[3];
rz(-1.118286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.36748019) q[2];
sx q[2];
rz(-1.8364649) q[2];
sx q[2];
rz(1.0779862) q[2];
rz(0.20400253) q[3];
sx q[3];
rz(-1.9670468) q[3];
sx q[3];
rz(1.3998122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3486977) q[0];
sx q[0];
rz(-2.5070511) q[0];
sx q[0];
rz(-3.0358411) q[0];
rz(-1.8271242) q[1];
sx q[1];
rz(-2.1038838) q[1];
sx q[1];
rz(2.7076142) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8815241) q[0];
sx q[0];
rz(-0.48598841) q[0];
sx q[0];
rz(1.6218779) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9527516) q[2];
sx q[2];
rz(-1.8318118) q[2];
sx q[2];
rz(-0.37303916) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8934007) q[1];
sx q[1];
rz(-1.6576644) q[1];
sx q[1];
rz(0.76201622) q[1];
rz(-pi) q[2];
rz(2.8777255) q[3];
sx q[3];
rz(-1.646606) q[3];
sx q[3];
rz(-2.3621462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5964261) q[2];
sx q[2];
rz(-0.20238987) q[2];
sx q[2];
rz(1.8983967) q[2];
rz(3.0498114) q[3];
sx q[3];
rz(-1.6538849) q[3];
sx q[3];
rz(1.53246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5246326) q[0];
sx q[0];
rz(-1.7175104) q[0];
sx q[0];
rz(2.4272163) q[0];
rz(2.3081035) q[1];
sx q[1];
rz(-2.1099213) q[1];
sx q[1];
rz(0.74976903) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0548328) q[0];
sx q[0];
rz(-1.0386837) q[0];
sx q[0];
rz(-0.52901308) q[0];
x q[1];
rz(-1.1718458) q[2];
sx q[2];
rz(-1.0064501) q[2];
sx q[2];
rz(2.7442748) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6121889) q[1];
sx q[1];
rz(-2.2893808) q[1];
sx q[1];
rz(2.8118333) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8213098) q[3];
sx q[3];
rz(-1.6289454) q[3];
sx q[3];
rz(0.15722903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3124776) q[2];
sx q[2];
rz(-1.1822327) q[2];
sx q[2];
rz(1.0324837) q[2];
rz(1.9898344) q[3];
sx q[3];
rz(-2.0572898) q[3];
sx q[3];
rz(1.1758218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51079232) q[0];
sx q[0];
rz(-1.3277338) q[0];
sx q[0];
rz(3.0603141) q[0];
rz(-2.6649113) q[1];
sx q[1];
rz(-1.3414693) q[1];
sx q[1];
rz(0.016990677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1611902) q[0];
sx q[0];
rz(-2.542164) q[0];
sx q[0];
rz(-0.73841788) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65036739) q[2];
sx q[2];
rz(-2.0135437) q[2];
sx q[2];
rz(2.3897417) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5754328) q[1];
sx q[1];
rz(-1.448307) q[1];
sx q[1];
rz(3.0585101) q[1];
rz(-3.105427) q[3];
sx q[3];
rz(-1.6828575) q[3];
sx q[3];
rz(-1.6074635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2396669) q[2];
sx q[2];
rz(-1.6096175) q[2];
sx q[2];
rz(1.7473183) q[2];
rz(1.4782921) q[3];
sx q[3];
rz(-2.6152488) q[3];
sx q[3];
rz(-2.8315376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33206853) q[0];
sx q[0];
rz(-2.9991751) q[0];
sx q[0];
rz(-0.19752565) q[0];
rz(-1.7212414) q[1];
sx q[1];
rz(-1.0528126) q[1];
sx q[1];
rz(-1.24409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8007788) q[0];
sx q[0];
rz(-1.4288623) q[0];
sx q[0];
rz(0.63951389) q[0];
rz(-2.317987) q[2];
sx q[2];
rz(-1.428229) q[2];
sx q[2];
rz(0.31566516) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1411775) q[1];
sx q[1];
rz(-1.1338935) q[1];
sx q[1];
rz(-1.8975329) q[1];
rz(-3.1116034) q[3];
sx q[3];
rz(-0.99195671) q[3];
sx q[3];
rz(2.8473583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.55599219) q[2];
sx q[2];
rz(-0.25871325) q[2];
sx q[2];
rz(-0.82536215) q[2];
rz(-0.9591006) q[3];
sx q[3];
rz(-0.93950713) q[3];
sx q[3];
rz(1.3033029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.30967228) q[0];
sx q[0];
rz(-0.29255022) q[0];
sx q[0];
rz(-1.1462611) q[0];
rz(-1.4347026) q[1];
sx q[1];
rz(-2.2271894) q[1];
sx q[1];
rz(0.31417876) q[1];
rz(-0.5774647) q[2];
sx q[2];
rz(-2.2784735) q[2];
sx q[2];
rz(1.4934595) q[2];
rz(0.17735986) q[3];
sx q[3];
rz(-1.0724417) q[3];
sx q[3];
rz(0.90870212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
