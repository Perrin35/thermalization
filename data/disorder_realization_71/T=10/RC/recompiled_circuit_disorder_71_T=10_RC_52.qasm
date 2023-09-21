OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(-2.3117476) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(-0.87632626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60183817) q[0];
sx q[0];
rz(-1.1798501) q[0];
sx q[0];
rz(-1.3215617) q[0];
rz(-1.8560156) q[2];
sx q[2];
rz(-0.60456317) q[2];
sx q[2];
rz(2.4100458) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6688924) q[1];
sx q[1];
rz(-1.2408222) q[1];
sx q[1];
rz(-0.015923576) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34378864) q[3];
sx q[3];
rz(-0.4292092) q[3];
sx q[3];
rz(-2.6478298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9709388) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(-2.412964) q[2];
rz(-2.6206) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(-2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347297) q[0];
sx q[0];
rz(-1.9711718) q[0];
sx q[0];
rz(1.1215425) q[0];
rz(-2.8858378) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(0.87444011) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16337285) q[0];
sx q[0];
rz(-1.997943) q[0];
sx q[0];
rz(2.8073729) q[0];
rz(2.9875056) q[2];
sx q[2];
rz(-2.140897) q[2];
sx q[2];
rz(0.64683435) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.892131) q[1];
sx q[1];
rz(-2.3453237) q[1];
sx q[1];
rz(-2.4888121) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3784834) q[3];
sx q[3];
rz(-2.8674539) q[3];
sx q[3];
rz(0.4916693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4008537) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(2.7056616) q[2];
rz(-2.46051) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(0.40288231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9044559) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(-2.0667734) q[0];
rz(0.83956051) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(2.7456465) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11055761) q[0];
sx q[0];
rz(-1.9651411) q[0];
sx q[0];
rz(-1.0236543) q[0];
x q[1];
rz(-3.1275438) q[2];
sx q[2];
rz(-1.2767681) q[2];
sx q[2];
rz(0.77169466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0242651) q[1];
sx q[1];
rz(-1.5612484) q[1];
sx q[1];
rz(2.4584103) q[1];
x q[2];
rz(-1.5941761) q[3];
sx q[3];
rz(-1.3153207) q[3];
sx q[3];
rz(1.3083096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6039156) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(-2.9023857) q[2];
rz(-0.075332969) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(1.8384365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72162119) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(-2.3216632) q[0];
rz(-2.6539102) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(-0.23342361) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963835) q[0];
sx q[0];
rz(-1.5102981) q[0];
sx q[0];
rz(1.66771) q[0];
rz(0.31077023) q[2];
sx q[2];
rz(-2.4756873) q[2];
sx q[2];
rz(0.13194612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9277716) q[1];
sx q[1];
rz(-1.908761) q[1];
sx q[1];
rz(-2.4073699) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-2.0908115) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(2.3941669) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(-2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915879) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(3.0773556) q[0];
rz(2.1977987) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-2.3805526) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8781389) q[0];
sx q[0];
rz(-1.8315151) q[0];
sx q[0];
rz(-0.031147416) q[0];
rz(-pi) q[1];
rz(0.76743482) q[2];
sx q[2];
rz(-0.27786294) q[2];
sx q[2];
rz(0.9741191) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.099247301) q[1];
sx q[1];
rz(-1.5464916) q[1];
sx q[1];
rz(-2.821032) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8588572) q[3];
sx q[3];
rz(-1.3021819) q[3];
sx q[3];
rz(-2.3408567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3349907) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(3.0495194) q[2];
rz(0.66172415) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(-1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46835607) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(3.1150505) q[0];
rz(-0.87310711) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(-0.32593265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78544261) q[0];
sx q[0];
rz(-2.4768562) q[0];
sx q[0];
rz(-0.63389969) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58165254) q[2];
sx q[2];
rz(-0.51917167) q[2];
sx q[2];
rz(-0.66228629) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1833916) q[1];
sx q[1];
rz(-2.8124053) q[1];
sx q[1];
rz(-2.7126461) q[1];
rz(0.10109191) q[3];
sx q[3];
rz(-1.5540082) q[3];
sx q[3];
rz(-1.3217795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.59763336) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(-1.4298965) q[3];
sx q[3];
rz(-1.1681898) q[3];
sx q[3];
rz(0.54106075) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8687246) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-0.72189271) q[0];
rz(-1.7294653) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(3.022335) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3463466) q[0];
sx q[0];
rz(-0.7659142) q[0];
sx q[0];
rz(-0.38608293) q[0];
rz(-pi) q[1];
rz(-0.97738816) q[2];
sx q[2];
rz(-2.9466629) q[2];
sx q[2];
rz(-2.5755663) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1489361) q[1];
sx q[1];
rz(-1.7033556) q[1];
sx q[1];
rz(-1.0692783) q[1];
rz(-pi) q[2];
rz(2.4339606) q[3];
sx q[3];
rz(-1.8579357) q[3];
sx q[3];
rz(-1.5054782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44935903) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(-1.3593486) q[2];
rz(-0.75891495) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(-0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.83157241) q[0];
sx q[0];
rz(-0.65905237) q[0];
sx q[0];
rz(2.0781562) q[0];
rz(2.8670782) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(0.88561052) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6016156) q[0];
sx q[0];
rz(-0.5479387) q[0];
sx q[0];
rz(2.6849296) q[0];
rz(-pi) q[1];
rz(3.1259414) q[2];
sx q[2];
rz(-0.99132292) q[2];
sx q[2];
rz(-1.5755115) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13247989) q[1];
sx q[1];
rz(-1.9630034) q[1];
sx q[1];
rz(-1.023804) q[1];
x q[2];
rz(0.76836821) q[3];
sx q[3];
rz(-1.0230912) q[3];
sx q[3];
rz(-0.87793575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72835913) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(2.0098861) q[2];
rz(1.0845832) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(-1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(1.0820748) q[0];
rz(1.8661631) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(2.0057604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58763114) q[0];
sx q[0];
rz(-0.56092867) q[0];
sx q[0];
rz(1.2558623) q[0];
rz(-1.5485974) q[2];
sx q[2];
rz(-2.0054842) q[2];
sx q[2];
rz(-1.7322025) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.75479924) q[1];
sx q[1];
rz(-2.7126185) q[1];
sx q[1];
rz(-3.0241443) q[1];
rz(-2.2215861) q[3];
sx q[3];
rz(-1.917106) q[3];
sx q[3];
rz(-1.3045834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0231126) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(2.2980799) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(-1.2783485) q[0];
rz(1.0247914) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(-1.1970253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.409317) q[0];
sx q[0];
rz(-1.1336375) q[0];
sx q[0];
rz(2.2073295) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0212101) q[2];
sx q[2];
rz(-1.3345846) q[2];
sx q[2];
rz(1.2506968) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4492606) q[1];
sx q[1];
rz(-2.3093866) q[1];
sx q[1];
rz(-0.87528054) q[1];
x q[2];
rz(0.61693807) q[3];
sx q[3];
rz(-2.5411798) q[3];
sx q[3];
rz(2.2850349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3060351) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(2.3416134) q[2];
rz(-1.9647313) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(0.95361382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70893127) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(2.6196383) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(0.37408806) q[2];
sx q[2];
rz(-1.903152) q[2];
sx q[2];
rz(-2.8852035) q[2];
rz(-2.5064777) q[3];
sx q[3];
rz(-2.1182346) q[3];
sx q[3];
rz(2.3253141) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
