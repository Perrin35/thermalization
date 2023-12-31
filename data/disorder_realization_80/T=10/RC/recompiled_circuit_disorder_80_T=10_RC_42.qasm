OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7497082) q[0];
sx q[0];
rz(-2.9449129) q[0];
sx q[0];
rz(-1.1893907) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(0.37766159) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.057214) q[0];
sx q[0];
rz(-0.14176653) q[0];
sx q[0];
rz(-1.0300107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2871735) q[2];
sx q[2];
rz(-1.5429153) q[2];
sx q[2];
rz(-1.4717799) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1805815) q[1];
sx q[1];
rz(-2.2716224) q[1];
sx q[1];
rz(-2.5518774) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8336485) q[3];
sx q[3];
rz(-0.43711284) q[3];
sx q[3];
rz(-1.2446158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.010014023) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(-0.67260355) q[2];
rz(2.9721695) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23068962) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(-0.22856523) q[0];
rz(-2.9810492) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(-0.28796089) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9790968) q[0];
sx q[0];
rz(-1.7485577) q[0];
sx q[0];
rz(-1.3773247) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1638853) q[2];
sx q[2];
rz(-1.2427254) q[2];
sx q[2];
rz(-1.6239945) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60625633) q[1];
sx q[1];
rz(-1.1621981) q[1];
sx q[1];
rz(2.4819083) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6091299) q[3];
sx q[3];
rz(-1.9209071) q[3];
sx q[3];
rz(-0.73327524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(1.6681558) q[2];
rz(-0.93747059) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(-0.37500769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8163452) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(2.8741799) q[0];
rz(1.422241) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(-0.95169383) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0591473) q[0];
sx q[0];
rz(-1.379231) q[0];
sx q[0];
rz(-3.0491769) q[0];
x q[1];
rz(-2.9159413) q[2];
sx q[2];
rz(-1.0661085) q[2];
sx q[2];
rz(-0.80863189) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.051085) q[1];
sx q[1];
rz(-1.2625492) q[1];
sx q[1];
rz(-2.3329263) q[1];
rz(-pi) q[2];
rz(0.53593105) q[3];
sx q[3];
rz(-2.1246353) q[3];
sx q[3];
rz(1.904408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13016985) q[2];
sx q[2];
rz(-1.6458076) q[2];
sx q[2];
rz(0.34417957) q[2];
rz(-2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0043871) q[0];
sx q[0];
rz(-1.4973649) q[0];
sx q[0];
rz(1.8970867) q[0];
rz(0.1291153) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(-2.7688162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6308206) q[0];
sx q[0];
rz(-1.456506) q[0];
sx q[0];
rz(2.3833016) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4332306) q[2];
sx q[2];
rz(-0.82651143) q[2];
sx q[2];
rz(2.7162958) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72153889) q[1];
sx q[1];
rz(-2.5905847) q[1];
sx q[1];
rz(1.0934456) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9863015) q[3];
sx q[3];
rz(-2.5875475) q[3];
sx q[3];
rz(-0.98741764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.63561511) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(-2.2367031) q[2];
rz(2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(-0.99159616) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47914094) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(-1.8537846) q[0];
rz(-1.4783391) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(-0.53422654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984765) q[0];
sx q[0];
rz(-1.2623708) q[0];
sx q[0];
rz(-0.31750676) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0787557) q[2];
sx q[2];
rz(-1.8740219) q[2];
sx q[2];
rz(-0.505503) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2248762) q[1];
sx q[1];
rz(-1.6460878) q[1];
sx q[1];
rz(-3.0776943) q[1];
x q[2];
rz(0.71877919) q[3];
sx q[3];
rz(-1.6703509) q[3];
sx q[3];
rz(-2.0841141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(-0.14349288) q[2];
rz(-1.3714553) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(-2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7025529) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-0.23705661) q[0];
rz(1.2409695) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(-1.3751078) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8538044) q[0];
sx q[0];
rz(-1.8332991) q[0];
sx q[0];
rz(0.35065035) q[0];
x q[1];
rz(0.74890045) q[2];
sx q[2];
rz(-1.5178711) q[2];
sx q[2];
rz(2.5089094) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68723893) q[1];
sx q[1];
rz(-2.7308186) q[1];
sx q[1];
rz(-0.06128581) q[1];
x q[2];
rz(-1.1477908) q[3];
sx q[3];
rz(-2.0241963) q[3];
sx q[3];
rz(0.51844937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53283006) q[2];
sx q[2];
rz(-0.71981788) q[2];
sx q[2];
rz(0.14870816) q[2];
rz(0.016629774) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(-2.9490411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6454813) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(2.2684229) q[0];
rz(-2.219615) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(-3.057664) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5859563) q[0];
sx q[0];
rz(-1.3480098) q[0];
sx q[0];
rz(2.1369834) q[0];
rz(-pi) q[1];
rz(-0.53972466) q[2];
sx q[2];
rz(-2.735609) q[2];
sx q[2];
rz(-0.27137953) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13087194) q[1];
sx q[1];
rz(-2.2166703) q[1];
sx q[1];
rz(2.3400397) q[1];
rz(-pi) q[2];
rz(-0.84079068) q[3];
sx q[3];
rz(-1.6222686) q[3];
sx q[3];
rz(-1.7712902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7579047) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(-2.596358) q[2];
rz(-2.7111354) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(-0.30495131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51171821) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(1.2114725) q[0];
rz(-2.1684872) q[1];
sx q[1];
rz(-0.5935697) q[1];
sx q[1];
rz(-0.94582311) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4771172) q[0];
sx q[0];
rz(-0.23453377) q[0];
sx q[0];
rz(-0.3420591) q[0];
rz(-pi) q[1];
rz(2.2875167) q[2];
sx q[2];
rz(-1.0208566) q[2];
sx q[2];
rz(-0.16888176) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1434053) q[1];
sx q[1];
rz(-0.37969509) q[1];
sx q[1];
rz(-3.0642964) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010057851) q[3];
sx q[3];
rz(-2.4754482) q[3];
sx q[3];
rz(-0.87107623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2716081) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(2.8611709) q[2];
rz(2.5366606) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(-2.159507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9873001) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(-2.5262685) q[0];
rz(-0.92957169) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(2.5659134) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79918039) q[0];
sx q[0];
rz(-0.77501446) q[0];
sx q[0];
rz(-0.70223017) q[0];
rz(-pi) q[1];
rz(-2.0909537) q[2];
sx q[2];
rz(-0.22503223) q[2];
sx q[2];
rz(-2.2573543) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3153783) q[1];
sx q[1];
rz(-1.678102) q[1];
sx q[1];
rz(-0.025749287) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9541295) q[3];
sx q[3];
rz(-1.0378569) q[3];
sx q[3];
rz(0.77123469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3118887) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(3.1196307) q[2];
rz(-2.9711376) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(-2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72572529) q[0];
sx q[0];
rz(-0.9265582) q[0];
sx q[0];
rz(-2.5073994) q[0];
rz(-0.028907396) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(-0.96910563) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4934851) q[0];
sx q[0];
rz(-0.095795184) q[0];
sx q[0];
rz(-2.2936054) q[0];
rz(-pi) q[1];
rz(-1.3766039) q[2];
sx q[2];
rz(-2.3079434) q[2];
sx q[2];
rz(3.0498691) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5574054) q[1];
sx q[1];
rz(-2.4483042) q[1];
sx q[1];
rz(0.28179817) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5892331) q[3];
sx q[3];
rz(-0.5061572) q[3];
sx q[3];
rz(1.5137223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6897631) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(-0.36007145) q[2];
rz(1.5277956) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(-0.56267363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2071028) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(-3.0974401) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(2.8051878) q[2];
sx q[2];
rz(-1.207926) q[2];
sx q[2];
rz(1.5628857) q[2];
rz(2.5526657) q[3];
sx q[3];
rz(-2.6116461) q[3];
sx q[3];
rz(-2.6245821) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
