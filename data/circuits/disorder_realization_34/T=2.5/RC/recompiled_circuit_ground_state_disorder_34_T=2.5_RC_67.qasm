OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.9095705) q[0];
sx q[0];
rz(2.4973305) q[0];
sx q[0];
rz(7.147026) q[0];
rz(4.7148352) q[1];
sx q[1];
rz(5.0532053) q[1];
sx q[1];
rz(2.4911575) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0599974) q[0];
sx q[0];
rz(-1.0733111) q[0];
sx q[0];
rz(1.803969) q[0];
rz(-pi) q[1];
rz(-3.0502599) q[2];
sx q[2];
rz(-2.8300081) q[2];
sx q[2];
rz(0.1908737) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0697223) q[1];
sx q[1];
rz(-2.156667) q[1];
sx q[1];
rz(-0.61572325) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4635383) q[3];
sx q[3];
rz(-1.4441731) q[3];
sx q[3];
rz(2.7789314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9170561) q[2];
sx q[2];
rz(-1.7660331) q[2];
sx q[2];
rz(1.7898233) q[2];
rz(0.81870643) q[3];
sx q[3];
rz(-1.6823781) q[3];
sx q[3];
rz(-1.5514577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.8390035) q[0];
sx q[0];
rz(-0.67794696) q[0];
sx q[0];
rz(2.5658521) q[0];
rz(-2.2585244) q[1];
sx q[1];
rz(-1.6022316) q[1];
sx q[1];
rz(-1.8227122) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3991845) q[0];
sx q[0];
rz(-2.3517646) q[0];
sx q[0];
rz(2.7862048) q[0];
rz(-pi) q[1];
rz(3.1167745) q[2];
sx q[2];
rz(-0.78486004) q[2];
sx q[2];
rz(-1.0944686) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3471862) q[1];
sx q[1];
rz(-2.3270895) q[1];
sx q[1];
rz(-2.0655355) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6448375) q[3];
sx q[3];
rz(-0.36128929) q[3];
sx q[3];
rz(-2.3827117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5004702) q[2];
sx q[2];
rz(-1.0176071) q[2];
sx q[2];
rz(-2.9024331) q[2];
rz(0.33254361) q[3];
sx q[3];
rz(-1.6859237) q[3];
sx q[3];
rz(-1.0913764) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34514937) q[0];
sx q[0];
rz(-2.394016) q[0];
sx q[0];
rz(0.24459608) q[0];
rz(2.7469475) q[1];
sx q[1];
rz(-0.60401812) q[1];
sx q[1];
rz(1.9371202) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3217425) q[0];
sx q[0];
rz(-1.5816147) q[0];
sx q[0];
rz(-0.011011684) q[0];
x q[1];
rz(2.7543805) q[2];
sx q[2];
rz(-2.0287598) q[2];
sx q[2];
rz(2.2004623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.50379163) q[1];
sx q[1];
rz(-2.0855013) q[1];
sx q[1];
rz(0.33333244) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7359636) q[3];
sx q[3];
rz(-2.0361216) q[3];
sx q[3];
rz(0.80287186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.34387732) q[2];
sx q[2];
rz(-0.98261967) q[2];
sx q[2];
rz(-0.24125153) q[2];
rz(1.7734211) q[3];
sx q[3];
rz(-1.3894812) q[3];
sx q[3];
rz(-2.1696137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5695246) q[0];
sx q[0];
rz(-1.854874) q[0];
sx q[0];
rz(-2.0401814) q[0];
rz(1.4934348) q[1];
sx q[1];
rz(-0.28683528) q[1];
sx q[1];
rz(-2.1410087) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7582486) q[0];
sx q[0];
rz(-2.7186518) q[0];
sx q[0];
rz(-1.3995885) q[0];
rz(0.0098044013) q[2];
sx q[2];
rz(-1.0145309) q[2];
sx q[2];
rz(-0.74273587) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8800671) q[1];
sx q[1];
rz(-1.6717922) q[1];
sx q[1];
rz(-1.5965881) q[1];
rz(-pi) q[2];
rz(-2.9977418) q[3];
sx q[3];
rz(-1.9884304) q[3];
sx q[3];
rz(-0.9506433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7857886) q[2];
sx q[2];
rz(-1.641909) q[2];
sx q[2];
rz(-0.22264063) q[2];
rz(1.6629793) q[3];
sx q[3];
rz(-0.40545774) q[3];
sx q[3];
rz(-1.8805898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.8292238) q[0];
sx q[0];
rz(-2.0033422) q[0];
sx q[0];
rz(0.19947048) q[0];
rz(1.5169187) q[1];
sx q[1];
rz(-1.805504) q[1];
sx q[1];
rz(-2.6768501) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6049395) q[0];
sx q[0];
rz(-2.8230328) q[0];
sx q[0];
rz(2.0037988) q[0];
rz(-1.6301441) q[2];
sx q[2];
rz(-1.3497102) q[2];
sx q[2];
rz(-0.57002744) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3252651) q[1];
sx q[1];
rz(-1.2817329) q[1];
sx q[1];
rz(-1.8613893) q[1];
rz(-1.5681276) q[3];
sx q[3];
rz(-0.87335372) q[3];
sx q[3];
rz(1.3476882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0988079) q[2];
sx q[2];
rz(-1.837919) q[2];
sx q[2];
rz(2.888077) q[2];
rz(0.86067307) q[3];
sx q[3];
rz(-2.6200675) q[3];
sx q[3];
rz(-1.1389987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0694224) q[0];
sx q[0];
rz(-1.6149898) q[0];
sx q[0];
rz(1.7813064) q[0];
rz(2.5379429) q[1];
sx q[1];
rz(-0.79704469) q[1];
sx q[1];
rz(0.48702494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0031155) q[0];
sx q[0];
rz(-0.91775187) q[0];
sx q[0];
rz(-1.4794635) q[0];
rz(-pi) q[1];
rz(1.4757882) q[2];
sx q[2];
rz(-0.83432799) q[2];
sx q[2];
rz(0.67259865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.071515171) q[1];
sx q[1];
rz(-1.2116451) q[1];
sx q[1];
rz(-0.17724322) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0661845) q[3];
sx q[3];
rz(-0.45201354) q[3];
sx q[3];
rz(1.3936123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6765678) q[2];
sx q[2];
rz(-2.1239231) q[2];
sx q[2];
rz(0.49099311) q[2];
rz(-0.53389126) q[3];
sx q[3];
rz(-0.54755727) q[3];
sx q[3];
rz(0.045031358) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.295306) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(1.2714161) q[0];
rz(2.2170587) q[1];
sx q[1];
rz(-1.144616) q[1];
sx q[1];
rz(-2.1163993) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3640396) q[0];
sx q[0];
rz(-1.784783) q[0];
sx q[0];
rz(0.63978934) q[0];
rz(-pi) q[1];
rz(-0.94550262) q[2];
sx q[2];
rz(-1.7472113) q[2];
sx q[2];
rz(0.27798619) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5793663) q[1];
sx q[1];
rz(-2.0831897) q[1];
sx q[1];
rz(-2.775869) q[1];
rz(-1.8162433) q[3];
sx q[3];
rz(-2.2625828) q[3];
sx q[3];
rz(-0.46190767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9292494) q[2];
sx q[2];
rz(-0.29444567) q[2];
sx q[2];
rz(2.219131) q[2];
rz(0.38659066) q[3];
sx q[3];
rz(-1.2873193) q[3];
sx q[3];
rz(-1.5132025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59231049) q[0];
sx q[0];
rz(-3.1070502) q[0];
sx q[0];
rz(0.70796815) q[0];
rz(0.51900807) q[1];
sx q[1];
rz(-0.52878562) q[1];
sx q[1];
rz(-2.4755898) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9922657) q[0];
sx q[0];
rz(-2.104787) q[0];
sx q[0];
rz(-0.1073246) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0254499) q[2];
sx q[2];
rz(-1.2525038) q[2];
sx q[2];
rz(1.4333985) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.076326536) q[1];
sx q[1];
rz(-0.90172036) q[1];
sx q[1];
rz(1.5764023) q[1];
rz(-pi) q[2];
rz(1.610393) q[3];
sx q[3];
rz(-0.61568135) q[3];
sx q[3];
rz(-2.1510268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4697326) q[2];
sx q[2];
rz(-1.1598776) q[2];
sx q[2];
rz(1.0603909) q[2];
rz(1.7139858) q[3];
sx q[3];
rz(-1.5771882) q[3];
sx q[3];
rz(2.3744627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1308744) q[0];
sx q[0];
rz(-2.3227782) q[0];
sx q[0];
rz(-2.4406216) q[0];
rz(-0.88841191) q[1];
sx q[1];
rz(-0.41524926) q[1];
sx q[1];
rz(1.433724) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5007223) q[0];
sx q[0];
rz(-1.6789915) q[0];
sx q[0];
rz(2.6263424) q[0];
rz(0.40590717) q[2];
sx q[2];
rz(-1.5957615) q[2];
sx q[2];
rz(-1.7737103) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4605324) q[1];
sx q[1];
rz(-1.5523445) q[1];
sx q[1];
rz(-3.1204719) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69207003) q[3];
sx q[3];
rz(-1.0502397) q[3];
sx q[3];
rz(-0.51560054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99011699) q[2];
sx q[2];
rz(-1.2274123) q[2];
sx q[2];
rz(1.4768451) q[2];
rz(-2.9448523) q[3];
sx q[3];
rz(-2.0048095) q[3];
sx q[3];
rz(0.93973947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6075341) q[0];
sx q[0];
rz(-1.5628096) q[0];
sx q[0];
rz(-2.0205355) q[0];
rz(-0.47814449) q[1];
sx q[1];
rz(-2.4604764) q[1];
sx q[1];
rz(-0.71472439) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9909458) q[0];
sx q[0];
rz(-2.0523221) q[0];
sx q[0];
rz(-1.7981547) q[0];
rz(-0.98866971) q[2];
sx q[2];
rz(-2.5091362) q[2];
sx q[2];
rz(-1.1256696) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0352382) q[1];
sx q[1];
rz(-2.0788621) q[1];
sx q[1];
rz(2.6324327) q[1];
rz(-pi) q[2];
rz(-0.86831324) q[3];
sx q[3];
rz(-0.83783093) q[3];
sx q[3];
rz(0.7759076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5152682) q[2];
sx q[2];
rz(-1.3268665) q[2];
sx q[2];
rz(1.2088306) q[2];
rz(-1.6995466) q[3];
sx q[3];
rz(-2.2550826) q[3];
sx q[3];
rz(3.076021) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6616853) q[0];
sx q[0];
rz(-1.8330782) q[0];
sx q[0];
rz(-0.080396419) q[0];
rz(-2.5936364) q[1];
sx q[1];
rz(-2.4520271) q[1];
sx q[1];
rz(1.3507623) q[1];
rz(1.1016079) q[2];
sx q[2];
rz(-1.3528578) q[2];
sx q[2];
rz(2.8982671) q[2];
rz(-1.9962541) q[3];
sx q[3];
rz(-1.1444935) q[3];
sx q[3];
rz(-1.7589699) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
