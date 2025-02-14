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
rz(1.3311812) q[0];
sx q[0];
rz(-1.5067195) q[0];
sx q[0];
rz(0.17568406) q[0];
rz(-2.076258) q[1];
sx q[1];
rz(-1.330436) q[1];
sx q[1];
rz(-2.4863844) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86533812) q[0];
sx q[0];
rz(-2.0575581) q[0];
sx q[0];
rz(-2.7810762) q[0];
rz(-pi) q[1];
rz(1.6263481) q[2];
sx q[2];
rz(-1.6695938) q[2];
sx q[2];
rz(-2.9227119) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.44819278) q[1];
sx q[1];
rz(-1.498135) q[1];
sx q[1];
rz(-0.10814114) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1813222) q[3];
sx q[3];
rz(-1.3283848) q[3];
sx q[3];
rz(-1.2038972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4483999) q[2];
sx q[2];
rz(-3.0312263) q[2];
sx q[2];
rz(-2.8006862) q[2];
rz(2.36813) q[3];
sx q[3];
rz(-2.1229459) q[3];
sx q[3];
rz(-0.10507467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4648436) q[0];
sx q[0];
rz(-2.8571547) q[0];
sx q[0];
rz(-0.25962096) q[0];
rz(-2.3530841) q[1];
sx q[1];
rz(-0.85616833) q[1];
sx q[1];
rz(-1.2451008) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21331025) q[0];
sx q[0];
rz(-1.6716372) q[0];
sx q[0];
rz(-0.95243246) q[0];
x q[1];
rz(-1.4711667) q[2];
sx q[2];
rz(-1.8848231) q[2];
sx q[2];
rz(0.043836029) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4739405) q[1];
sx q[1];
rz(-0.86939916) q[1];
sx q[1];
rz(-0.46469537) q[1];
x q[2];
rz(2.8471242) q[3];
sx q[3];
rz(-2.614902) q[3];
sx q[3];
rz(-1.1931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9022687) q[2];
sx q[2];
rz(-2.6332492) q[2];
sx q[2];
rz(2.2321205) q[2];
rz(-2.4746312) q[3];
sx q[3];
rz(-1.704155) q[3];
sx q[3];
rz(0.10233574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42720902) q[0];
sx q[0];
rz(-2.7294071) q[0];
sx q[0];
rz(-1.3877731) q[0];
rz(2.1448403) q[1];
sx q[1];
rz(-0.69147384) q[1];
sx q[1];
rz(2.6549285) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.023611) q[0];
sx q[0];
rz(-1.1661464) q[0];
sx q[0];
rz(2.9756474) q[0];
rz(-2.2816554) q[2];
sx q[2];
rz(-0.87274466) q[2];
sx q[2];
rz(2.6280478) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8139517) q[1];
sx q[1];
rz(-2.608229) q[1];
sx q[1];
rz(2.9017828) q[1];
x q[2];
rz(-0.14778845) q[3];
sx q[3];
rz(-2.1923965) q[3];
sx q[3];
rz(-0.10336598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3643058) q[2];
sx q[2];
rz(-1.2591668) q[2];
sx q[2];
rz(-2.6510009) q[2];
rz(2.3250438) q[3];
sx q[3];
rz(-2.3994763) q[3];
sx q[3];
rz(2.1013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33646026) q[0];
sx q[0];
rz(-1.8959624) q[0];
sx q[0];
rz(-0.98139393) q[0];
rz(-2.7384752) q[1];
sx q[1];
rz(-0.50794452) q[1];
sx q[1];
rz(-2.6037762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92141672) q[0];
sx q[0];
rz(-2.8517637) q[0];
sx q[0];
rz(-1.0224065) q[0];
x q[1];
rz(1.4188624) q[2];
sx q[2];
rz(-2.4151994) q[2];
sx q[2];
rz(1.6223729) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0742948) q[1];
sx q[1];
rz(-0.63772574) q[1];
sx q[1];
rz(0.89115759) q[1];
rz(-pi) q[2];
rz(-0.89610164) q[3];
sx q[3];
rz(-1.3235914) q[3];
sx q[3];
rz(-0.02442115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8889403) q[2];
sx q[2];
rz(-2.6344968) q[2];
sx q[2];
rz(-2.0895152) q[2];
rz(-0.55225152) q[3];
sx q[3];
rz(-1.4058607) q[3];
sx q[3];
rz(1.9210057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.030516142) q[0];
sx q[0];
rz(-1.80856) q[0];
sx q[0];
rz(-0.037574969) q[0];
rz(-2.3887718) q[1];
sx q[1];
rz(-1.5767187) q[1];
sx q[1];
rz(-3.0573696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5355229) q[0];
sx q[0];
rz(-1.3227934) q[0];
sx q[0];
rz(-2.1169784) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4677643) q[2];
sx q[2];
rz(-2.6884656) q[2];
sx q[2];
rz(1.782531) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3854691) q[1];
sx q[1];
rz(-1.1471355) q[1];
sx q[1];
rz(3.0961354) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0952206) q[3];
sx q[3];
rz(-1.7313469) q[3];
sx q[3];
rz(-0.35453803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5273744) q[2];
sx q[2];
rz(-2.1671593) q[2];
sx q[2];
rz(-0.43240377) q[2];
rz(-1.6480986) q[3];
sx q[3];
rz(-1.3727539) q[3];
sx q[3];
rz(-2.4026332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768141) q[0];
sx q[0];
rz(-1.6141163) q[0];
sx q[0];
rz(0.76171184) q[0];
rz(1.0234045) q[1];
sx q[1];
rz(-0.49097148) q[1];
sx q[1];
rz(2.8010119) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7135431) q[0];
sx q[0];
rz(-0.78589702) q[0];
sx q[0];
rz(-0.12554306) q[0];
x q[1];
rz(-1.4144142) q[2];
sx q[2];
rz(-2.1264653) q[2];
sx q[2];
rz(-1.2909691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2557166) q[1];
sx q[1];
rz(-1.9105366) q[1];
sx q[1];
rz(-0.28743108) q[1];
x q[2];
rz(-0.16400985) q[3];
sx q[3];
rz(-2.384438) q[3];
sx q[3];
rz(-1.3379768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9684101) q[2];
sx q[2];
rz(-2.0081655) q[2];
sx q[2];
rz(-0.75562149) q[2];
rz(-0.38718811) q[3];
sx q[3];
rz(-1.7500992) q[3];
sx q[3];
rz(0.79425991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8887535) q[0];
sx q[0];
rz(-3.0536953) q[0];
sx q[0];
rz(1.9320236) q[0];
rz(1.5225438) q[1];
sx q[1];
rz(-0.77170283) q[1];
sx q[1];
rz(-1.3873772) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98859859) q[0];
sx q[0];
rz(-1.8805046) q[0];
sx q[0];
rz(1.9266846) q[0];
x q[1];
rz(2.8262994) q[2];
sx q[2];
rz(-2.1209426) q[2];
sx q[2];
rz(1.2619051) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55303326) q[1];
sx q[1];
rz(-2.297595) q[1];
sx q[1];
rz(2.5735106) q[1];
rz(-pi) q[2];
rz(-1.8903743) q[3];
sx q[3];
rz(-2.0154023) q[3];
sx q[3];
rz(2.3426659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9428955) q[2];
sx q[2];
rz(-0.2936475) q[2];
sx q[2];
rz(0.4100619) q[2];
rz(2.6041218) q[3];
sx q[3];
rz(-2.0987857) q[3];
sx q[3];
rz(-1.0679831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-0.10802565) q[0];
sx q[0];
rz(-1.2355868) q[0];
sx q[0];
rz(-2.1349452) q[0];
rz(2.0058696) q[1];
sx q[1];
rz(-0.4907116) q[1];
sx q[1];
rz(2.8699285) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1170809) q[0];
sx q[0];
rz(-0.70444626) q[0];
sx q[0];
rz(-2.4570877) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4681787) q[2];
sx q[2];
rz(-2.1072763) q[2];
sx q[2];
rz(-1.7458037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0783196) q[1];
sx q[1];
rz(-2.1542179) q[1];
sx q[1];
rz(2.2607445) q[1];
rz(-1.7944505) q[3];
sx q[3];
rz(-0.8394548) q[3];
sx q[3];
rz(-0.37089965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9378308) q[2];
sx q[2];
rz(-0.93272847) q[2];
sx q[2];
rz(-2.0358098) q[2];
rz(1.743099) q[3];
sx q[3];
rz(-2.1574557) q[3];
sx q[3];
rz(-2.2891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49649134) q[0];
sx q[0];
rz(-2.190525) q[0];
sx q[0];
rz(2.7880461) q[0];
rz(-1.4191267) q[1];
sx q[1];
rz(-1.5357693) q[1];
sx q[1];
rz(-0.062573418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.94189) q[0];
sx q[0];
rz(-1.0257922) q[0];
sx q[0];
rz(2.1502994) q[0];
rz(-pi) q[1];
rz(-2.563971) q[2];
sx q[2];
rz(-2.2934539) q[2];
sx q[2];
rz(-1.8236782) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3741604) q[1];
sx q[1];
rz(-1.8783775) q[1];
sx q[1];
rz(0.49616431) q[1];
x q[2];
rz(-1.4776286) q[3];
sx q[3];
rz(-2.2312632) q[3];
sx q[3];
rz(2.4371393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8211557) q[2];
sx q[2];
rz(-1.374036) q[2];
sx q[2];
rz(2.721411) q[2];
rz(0.34475103) q[3];
sx q[3];
rz(-0.77778608) q[3];
sx q[3];
rz(2.2043118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47278136) q[0];
sx q[0];
rz(-2.9720699) q[0];
sx q[0];
rz(1.8602759) q[0];
rz(1.5538813) q[1];
sx q[1];
rz(-0.80931598) q[1];
sx q[1];
rz(2.4344427) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2408694) q[0];
sx q[0];
rz(-1.5587731) q[0];
sx q[0];
rz(2.2012086) q[0];
x q[1];
rz(2.1936839) q[2];
sx q[2];
rz(-2.4595408) q[2];
sx q[2];
rz(-0.94833224) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7950629) q[1];
sx q[1];
rz(-1.7938385) q[1];
sx q[1];
rz(-2.8233379) q[1];
rz(-pi) q[2];
rz(2.1955802) q[3];
sx q[3];
rz(-1.0407699) q[3];
sx q[3];
rz(-1.8998787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3215434) q[2];
sx q[2];
rz(-2.0115435) q[2];
sx q[2];
rz(0.496544) q[2];
rz(1.8494362) q[3];
sx q[3];
rz(-2.8437331) q[3];
sx q[3];
rz(2.7636102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.9847066) q[0];
sx q[0];
rz(-0.26837415) q[0];
sx q[0];
rz(-1.0347086) q[0];
rz(-2.5490419) q[1];
sx q[1];
rz(-2.4609346) q[1];
sx q[1];
rz(1.5998283) q[1];
rz(0.040475673) q[2];
sx q[2];
rz(-1.802609) q[2];
sx q[2];
rz(2.6229057) q[2];
rz(1.7962528) q[3];
sx q[3];
rz(-0.97889789) q[3];
sx q[3];
rz(-1.2759664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
