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
rz(-0.77684075) q[0];
sx q[0];
rz(2.2651894) q[0];
sx q[0];
rz(9.6776008) q[0];
rz(1.1480992) q[1];
sx q[1];
rz(-2.2263081) q[1];
sx q[1];
rz(-0.80438703) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641805) q[0];
sx q[0];
rz(-1.2578811) q[0];
sx q[0];
rz(1.1699647) q[0];
rz(-1.8041274) q[2];
sx q[2];
rz(-1.1510282) q[2];
sx q[2];
rz(2.6884318) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.731718) q[1];
sx q[1];
rz(-1.9148286) q[1];
sx q[1];
rz(-2.7469977) q[1];
rz(-0.31320235) q[3];
sx q[3];
rz(-0.36010183) q[3];
sx q[3];
rz(-1.859377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.09314166) q[2];
sx q[2];
rz(-0.96258771) q[2];
sx q[2];
rz(0.87240458) q[2];
rz(2.8060272) q[3];
sx q[3];
rz(-1.395697) q[3];
sx q[3];
rz(0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63808477) q[0];
sx q[0];
rz(-2.8758949) q[0];
sx q[0];
rz(0.22915325) q[0];
rz(0.19042641) q[1];
sx q[1];
rz(-1.3399905) q[1];
sx q[1];
rz(0.71358877) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.021727) q[0];
sx q[0];
rz(-1.3029126) q[0];
sx q[0];
rz(-2.2557206) q[0];
x q[1];
rz(-2.607478) q[2];
sx q[2];
rz(-0.91098173) q[2];
sx q[2];
rz(2.5995863) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4137607) q[1];
sx q[1];
rz(-1.276386) q[1];
sx q[1];
rz(2.4319699) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9715367) q[3];
sx q[3];
rz(-0.77203686) q[3];
sx q[3];
rz(2.7626729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68636346) q[2];
sx q[2];
rz(-0.31193048) q[2];
sx q[2];
rz(2.6658106) q[2];
rz(-2.0717715) q[3];
sx q[3];
rz(-1.0082303) q[3];
sx q[3];
rz(-0.76655918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0692724) q[0];
sx q[0];
rz(-1.9444436) q[0];
sx q[0];
rz(-1.9092165) q[0];
rz(2.6909289) q[1];
sx q[1];
rz(-1.181299) q[1];
sx q[1];
rz(2.252069) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7018676) q[0];
sx q[0];
rz(-2.151818) q[0];
sx q[0];
rz(2.0608611) q[0];
rz(0.81368229) q[2];
sx q[2];
rz(-2.431834) q[2];
sx q[2];
rz(-2.8738632) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3308308) q[1];
sx q[1];
rz(-2.2137801) q[1];
sx q[1];
rz(-0.26442702) q[1];
rz(-1.9185478) q[3];
sx q[3];
rz(-1.9257716) q[3];
sx q[3];
rz(2.6048911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4543317) q[2];
sx q[2];
rz(-0.4250409) q[2];
sx q[2];
rz(-1.2588151) q[2];
rz(1.9076294) q[3];
sx q[3];
rz(-2.0089269) q[3];
sx q[3];
rz(-0.0080571938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4312129) q[0];
sx q[0];
rz(-0.90573913) q[0];
sx q[0];
rz(-1.6718965) q[0];
rz(-1.1471033) q[1];
sx q[1];
rz(-2.1397782) q[1];
sx q[1];
rz(2.240644) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16580627) q[0];
sx q[0];
rz(-2.8167509) q[0];
sx q[0];
rz(2.8715517) q[0];
x q[1];
rz(-0.9587165) q[2];
sx q[2];
rz(-1.0462073) q[2];
sx q[2];
rz(-0.9768578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9228153) q[1];
sx q[1];
rz(-0.73000693) q[1];
sx q[1];
rz(-0.14132146) q[1];
x q[2];
rz(1.4114831) q[3];
sx q[3];
rz(-0.84419227) q[3];
sx q[3];
rz(-2.5558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6441696) q[2];
sx q[2];
rz(-1.3966509) q[2];
sx q[2];
rz(-0.72745848) q[2];
rz(0.71508956) q[3];
sx q[3];
rz(-2.6810724) q[3];
sx q[3];
rz(0.49811825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46399507) q[0];
sx q[0];
rz(-1.7514739) q[0];
sx q[0];
rz(-2.4053307) q[0];
rz(-2.5212506) q[1];
sx q[1];
rz(-2.5963929) q[1];
sx q[1];
rz(-1.5249407) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10337457) q[0];
sx q[0];
rz(-3.0854825) q[0];
sx q[0];
rz(2.9655064) q[0];
rz(-pi) q[1];
rz(2.1651046) q[2];
sx q[2];
rz(-0.39820489) q[2];
sx q[2];
rz(-1.2334241) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8680537) q[1];
sx q[1];
rz(-2.2333849) q[1];
sx q[1];
rz(0.48133793) q[1];
rz(2.3297263) q[3];
sx q[3];
rz(-0.97211876) q[3];
sx q[3];
rz(-0.52385274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5357431) q[2];
sx q[2];
rz(-0.77631408) q[2];
sx q[2];
rz(2.1898451) q[2];
rz(0.4176628) q[3];
sx q[3];
rz(-0.2499191) q[3];
sx q[3];
rz(1.9865659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34053892) q[0];
sx q[0];
rz(-1.8908353) q[0];
sx q[0];
rz(0.75403768) q[0];
rz(-2.2484089) q[1];
sx q[1];
rz(-2.1008396) q[1];
sx q[1];
rz(0.16597861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5170068) q[0];
sx q[0];
rz(-0.69147325) q[0];
sx q[0];
rz(0.16384478) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7506261) q[2];
sx q[2];
rz(-2.0702614) q[2];
sx q[2];
rz(-1.7714372) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41478911) q[1];
sx q[1];
rz(-0.94894743) q[1];
sx q[1];
rz(1.0110823) q[1];
x q[2];
rz(-1.2955722) q[3];
sx q[3];
rz(-1.4365734) q[3];
sx q[3];
rz(1.8695267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37504998) q[2];
sx q[2];
rz(-1.1376209) q[2];
sx q[2];
rz(0.005391187) q[2];
rz(-3.052875) q[3];
sx q[3];
rz(-1.5327449) q[3];
sx q[3];
rz(-2.3458792) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2583112) q[0];
sx q[0];
rz(-1.2092051) q[0];
sx q[0];
rz(-0.2051556) q[0];
rz(0.908665) q[1];
sx q[1];
rz(-2.2712207) q[1];
sx q[1];
rz(0.044513449) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1703029) q[0];
sx q[0];
rz(-1.7006386) q[0];
sx q[0];
rz(1.4128039) q[0];
rz(-1.3088063) q[2];
sx q[2];
rz(-0.14523187) q[2];
sx q[2];
rz(1.1373718) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14340885) q[1];
sx q[1];
rz(-2.1381731) q[1];
sx q[1];
rz(-0.25556775) q[1];
rz(-pi) q[2];
rz(1.4511075) q[3];
sx q[3];
rz(-1.5933678) q[3];
sx q[3];
rz(2.165497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32555106) q[2];
sx q[2];
rz(-0.52246919) q[2];
sx q[2];
rz(0.52118707) q[2];
rz(2.9442287) q[3];
sx q[3];
rz(-1.3857931) q[3];
sx q[3];
rz(-0.50281966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407783) q[0];
sx q[0];
rz(-0.1939119) q[0];
sx q[0];
rz(-0.25522301) q[0];
rz(0.1420282) q[1];
sx q[1];
rz(-1.4165001) q[1];
sx q[1];
rz(2.7878888) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3981374) q[0];
sx q[0];
rz(-2.6530735) q[0];
sx q[0];
rz(2.3228881) q[0];
x q[1];
rz(2.1944517) q[2];
sx q[2];
rz(-1.1532239) q[2];
sx q[2];
rz(0.13969487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8579944) q[1];
sx q[1];
rz(-1.2322752) q[1];
sx q[1];
rz(2.0320862) q[1];
rz(0.74905101) q[3];
sx q[3];
rz(-1.1370249) q[3];
sx q[3];
rz(-1.6925136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8346617) q[2];
sx q[2];
rz(-2.8262409) q[2];
sx q[2];
rz(1.6174512) q[2];
rz(2.2751685) q[3];
sx q[3];
rz(-1.6774991) q[3];
sx q[3];
rz(0.58969897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6397112) q[0];
sx q[0];
rz(-2.9852133) q[0];
sx q[0];
rz(-1.7891275) q[0];
rz(-1.9068708) q[1];
sx q[1];
rz(-0.23209485) q[1];
sx q[1];
rz(0.51188767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5049625) q[0];
sx q[0];
rz(-1.413336) q[0];
sx q[0];
rz(3.0877293) q[0];
x q[1];
rz(2.9779766) q[2];
sx q[2];
rz(-0.77218845) q[2];
sx q[2];
rz(-0.26547394) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7831969) q[1];
sx q[1];
rz(-1.8640717) q[1];
sx q[1];
rz(-2.6986928) q[1];
rz(2.1549822) q[3];
sx q[3];
rz(-1.8711149) q[3];
sx q[3];
rz(0.14726135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1341683) q[2];
sx q[2];
rz(-1.8210501) q[2];
sx q[2];
rz(2.7443938) q[2];
rz(-2.126179) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(-0.5526244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39887244) q[0];
sx q[0];
rz(-1.2074559) q[0];
sx q[0];
rz(0.43750986) q[0];
rz(0.73879009) q[1];
sx q[1];
rz(-2.3414325) q[1];
sx q[1];
rz(1.4791666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9112944) q[0];
sx q[0];
rz(-2.0574967) q[0];
sx q[0];
rz(-1.507797) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9301038) q[2];
sx q[2];
rz(-1.6793161) q[2];
sx q[2];
rz(-1.4695652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91875118) q[1];
sx q[1];
rz(-1.6770308) q[1];
sx q[1];
rz(-1.6332455) q[1];
rz(2.2774901) q[3];
sx q[3];
rz(-2.2165944) q[3];
sx q[3];
rz(1.0458664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9208357) q[2];
sx q[2];
rz(-1.7831384) q[2];
sx q[2];
rz(-2.9662568) q[2];
rz(-1.0026503) q[3];
sx q[3];
rz(-0.23313871) q[3];
sx q[3];
rz(1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6709082) q[0];
sx q[0];
rz(-1.6841472) q[0];
sx q[0];
rz(2.0367784) q[0];
rz(-0.15057527) q[1];
sx q[1];
rz(-2.2655948) q[1];
sx q[1];
rz(1.2949952) q[1];
rz(1.5631093) q[2];
sx q[2];
rz(-1.3725217) q[2];
sx q[2];
rz(1.3524806) q[2];
rz(-1.8035268) q[3];
sx q[3];
rz(-1.3265811) q[3];
sx q[3];
rz(-1.4221232) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
