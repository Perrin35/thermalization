OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.75582957) q[0];
sx q[0];
rz(-1.4094149) q[0];
sx q[0];
rz(-0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298676) q[0];
sx q[0];
rz(-1.6352788) q[0];
sx q[0];
rz(0.026260016) q[0];
x q[1];
rz(2.4550081) q[2];
sx q[2];
rz(-1.540544) q[2];
sx q[2];
rz(-0.1498915) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20323682) q[1];
sx q[1];
rz(-2.6563546) q[1];
sx q[1];
rz(-0.39802246) q[1];
rz(2.9362039) q[3];
sx q[3];
rz(-2.5377512) q[3];
sx q[3];
rz(2.3103034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0191779) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(-1.4953556) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0533957) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(0.85900599) q[0];
rz(-2.7711218) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(-1.6765615) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252936) q[0];
sx q[0];
rz(-1.3536317) q[0];
sx q[0];
rz(0.3152245) q[0];
rz(-pi) q[1];
rz(-1.2334521) q[2];
sx q[2];
rz(-2.0806899) q[2];
sx q[2];
rz(-0.22491977) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4620004) q[1];
sx q[1];
rz(-0.781171) q[1];
sx q[1];
rz(0.52266927) q[1];
rz(0.46918418) q[3];
sx q[3];
rz(-0.57933925) q[3];
sx q[3];
rz(-2.7964696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7039965) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-2.583288) q[2];
rz(-1.0393418) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(-0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52477437) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-2.9779789) q[0];
rz(-0.71584654) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(1.7680426) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87953075) q[0];
sx q[0];
rz(-2.6143605) q[0];
sx q[0];
rz(3.0920045) q[0];
rz(-1.9203556) q[2];
sx q[2];
rz(-2.3974843) q[2];
sx q[2];
rz(-1.0457872) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92257567) q[1];
sx q[1];
rz(-2.1891928) q[1];
sx q[1];
rz(0.69505691) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3397066) q[3];
sx q[3];
rz(-2.826414) q[3];
sx q[3];
rz(1.2847881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13087656) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(1.019657) q[2];
rz(2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(-0.61292928) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402128) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(0.83909488) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(2.531321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3195517) q[0];
sx q[0];
rz(-0.10986957) q[0];
sx q[0];
rz(-1.6739474) q[0];
x q[1];
rz(1.1876112) q[2];
sx q[2];
rz(-1.2387347) q[2];
sx q[2];
rz(-0.039475723) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3210541) q[1];
sx q[1];
rz(-0.74756261) q[1];
sx q[1];
rz(1.9588406) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88704349) q[3];
sx q[3];
rz(-2.4409557) q[3];
sx q[3];
rz(-2.5168583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6976167) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(-2.234263) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6326555) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(1.6089815) q[0];
rz(2.4081047) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(1.3132494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816313) q[0];
sx q[0];
rz(-0.93802035) q[0];
sx q[0];
rz(2.7347793) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1057642) q[2];
sx q[2];
rz(-1.5140859) q[2];
sx q[2];
rz(2.5113475) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.207706) q[1];
sx q[1];
rz(-1.943214) q[1];
sx q[1];
rz(-1.1697342) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5562708) q[3];
sx q[3];
rz(-2.1748741) q[3];
sx q[3];
rz(0.6795336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0017172) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(2.4772947) q[2];
rz(-0.70513606) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(3.0378708) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076684549) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(-0.054811906) q[0];
rz(-2.1272155) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(2.6584113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7532363) q[0];
sx q[0];
rz(-1.6325132) q[0];
sx q[0];
rz(1.389099) q[0];
x q[1];
rz(-1.6364355) q[2];
sx q[2];
rz(-2.6780431) q[2];
sx q[2];
rz(2.133873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3997765) q[1];
sx q[1];
rz(-2.2175334) q[1];
sx q[1];
rz(0.6439376) q[1];
rz(-pi) q[2];
rz(-1.2383078) q[3];
sx q[3];
rz(-1.4629435) q[3];
sx q[3];
rz(0.40963848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9873535) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(-1.9412458) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744556) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(-0.95463395) q[0];
rz(2.5700991) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(2.8894997) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1958256) q[0];
sx q[0];
rz(-0.44647631) q[0];
sx q[0];
rz(-2.8004942) q[0];
x q[1];
rz(2.3977631) q[2];
sx q[2];
rz(-0.93223909) q[2];
sx q[2];
rz(0.078787412) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8512307) q[1];
sx q[1];
rz(-2.2154659) q[1];
sx q[1];
rz(-2.6840997) q[1];
rz(2.0310568) q[3];
sx q[3];
rz(-1.2292976) q[3];
sx q[3];
rz(-0.5865435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38348848) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(-0.90448109) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7082108) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(-1.7277539) q[0];
rz(0.66043234) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(1.7699014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6721281) q[0];
sx q[0];
rz(-0.95333316) q[0];
sx q[0];
rz(-3.0332964) q[0];
rz(-pi) q[1];
rz(1.3604836) q[2];
sx q[2];
rz(-0.82092199) q[2];
sx q[2];
rz(-2.6471777) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5270556) q[1];
sx q[1];
rz(-1.0053047) q[1];
sx q[1];
rz(1.2177699) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4189818) q[3];
sx q[3];
rz(-1.2137128) q[3];
sx q[3];
rz(0.15541542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82683212) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(0.39789847) q[2];
rz(0.49063101) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9234377) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(2.9072705) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(0.27059069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71352495) q[0];
sx q[0];
rz(-1.7881835) q[0];
sx q[0];
rz(-2.9146951) q[0];
rz(-pi) q[1];
rz(1.9174689) q[2];
sx q[2];
rz(-1.3895831) q[2];
sx q[2];
rz(2.5172362) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2724243) q[1];
sx q[1];
rz(-0.13726343) q[1];
sx q[1];
rz(-0.024172524) q[1];
rz(-pi) q[2];
rz(-1.0911646) q[3];
sx q[3];
rz(-1.661146) q[3];
sx q[3];
rz(0.8271715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33971912) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(1.489893) q[2];
rz(0.70288944) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(-1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(0.15429601) q[0];
rz(2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-2.399209) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48001227) q[0];
sx q[0];
rz(-2.1614657) q[0];
sx q[0];
rz(-2.3175879) q[0];
rz(-pi) q[1];
rz(2.7102091) q[2];
sx q[2];
rz(-2.1200392) q[2];
sx q[2];
rz(-1.0647578) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.91126373) q[1];
sx q[1];
rz(-1.2250369) q[1];
sx q[1];
rz(1.5826349) q[1];
rz(-pi) q[2];
rz(-0.14245716) q[3];
sx q[3];
rz(-1.5420621) q[3];
sx q[3];
rz(-1.8491668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2877038) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5819736) q[2];
rz(-0.21073267) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(-0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29006526) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(-1.7998981) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(1.2269536) q[2];
sx q[2];
rz(-2.362102) q[2];
sx q[2];
rz(-1.4951928) q[2];
rz(-0.7704173) q[3];
sx q[3];
rz(-2.9478248) q[3];
sx q[3];
rz(-0.41268681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
