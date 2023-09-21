OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9956545) q[0];
sx q[0];
rz(-0.50322682) q[0];
sx q[0];
rz(2.4174262) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(0.78483265) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0137579) q[0];
sx q[0];
rz(-0.36359596) q[0];
sx q[0];
rz(-2.512393) q[0];
rz(-pi) q[1];
rz(-1.6720812) q[2];
sx q[2];
rz(-1.1541608) q[2];
sx q[2];
rz(3.0243304) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29014978) q[1];
sx q[1];
rz(-0.57934299) q[1];
sx q[1];
rz(1.1555374) q[1];
rz(-1.647244) q[3];
sx q[3];
rz(-0.41562286) q[3];
sx q[3];
rz(-1.3941744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.589754) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(3.0736249) q[2];
rz(-0.12456482) q[3];
sx q[3];
rz(-2.8187276) q[3];
sx q[3];
rz(-1.7547866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.92000604) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(0.24366972) q[0];
rz(-0.63175732) q[1];
sx q[1];
rz(-1.4032204) q[1];
sx q[1];
rz(-1.7858645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6512017) q[0];
sx q[0];
rz(-0.25679195) q[0];
sx q[0];
rz(1.4635565) q[0];
x q[1];
rz(1.182105) q[2];
sx q[2];
rz(-1.2388065) q[2];
sx q[2];
rz(1.7322844) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4888549) q[1];
sx q[1];
rz(-1.4421717) q[1];
sx q[1];
rz(2.0205523) q[1];
x q[2];
rz(-0.59216604) q[3];
sx q[3];
rz(-1.6093996) q[3];
sx q[3];
rz(-1.6174699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-0.9920384) q[2];
sx q[2];
rz(0.24965723) q[2];
rz(-2.6349973) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(0.33199582) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8963985) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(-2.2431592) q[0];
rz(1.8067182) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(-1.867884) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1153591) q[0];
sx q[0];
rz(-0.44239487) q[0];
sx q[0];
rz(-2.1917403) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31295915) q[2];
sx q[2];
rz(-2.8319781) q[2];
sx q[2];
rz(1.6329873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85325235) q[1];
sx q[1];
rz(-0.77008343) q[1];
sx q[1];
rz(-2.6283162) q[1];
x q[2];
rz(-2.0687194) q[3];
sx q[3];
rz(-0.91091279) q[3];
sx q[3];
rz(-2.4592196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0818103) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(-0.95345062) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(-0.22687337) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2801441) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(2.8934073) q[0];
rz(-1.0379627) q[1];
sx q[1];
rz(-1.1231517) q[1];
sx q[1];
rz(-0.074137069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.603133) q[0];
sx q[0];
rz(-2.5939301) q[0];
sx q[0];
rz(1.0663701) q[0];
x q[1];
rz(-0.89216994) q[2];
sx q[2];
rz(-1.8461508) q[2];
sx q[2];
rz(-0.25472578) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84872765) q[1];
sx q[1];
rz(-1.8855727) q[1];
sx q[1];
rz(0.12401144) q[1];
x q[2];
rz(-1.7449964) q[3];
sx q[3];
rz(-1.8074236) q[3];
sx q[3];
rz(0.06121204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(3.1029491) q[2];
rz(-2.1679227) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.5357635) q[0];
sx q[0];
rz(-1.3624396) q[0];
rz(-0.81659395) q[1];
sx q[1];
rz(-1.8530308) q[1];
sx q[1];
rz(-1.978925) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1122738) q[0];
sx q[0];
rz(-0.11538878) q[0];
sx q[0];
rz(-1.3089887) q[0];
rz(0.16935279) q[2];
sx q[2];
rz(-1.8893818) q[2];
sx q[2];
rz(2.7404075) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34449023) q[1];
sx q[1];
rz(-1.3791729) q[1];
sx q[1];
rz(1.2002798) q[1];
x q[2];
rz(-2.707162) q[3];
sx q[3];
rz(-1.0940922) q[3];
sx q[3];
rz(1.3694976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(2.999372) q[2];
rz(2.2375315) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(-0.18946762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11480039) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(1.6739155) q[0];
rz(-0.57178512) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(0.30803672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338617) q[0];
sx q[0];
rz(-2.9836285) q[0];
sx q[0];
rz(0.64361848) q[0];
rz(-pi) q[1];
rz(-0.96254827) q[2];
sx q[2];
rz(-2.0934009) q[2];
sx q[2];
rz(0.42524291) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8923556) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(1.874079) q[1];
rz(2.8940053) q[3];
sx q[3];
rz(-0.35841225) q[3];
sx q[3];
rz(-0.7022411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77928153) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(-1.9936838) q[2];
rz(0.71427304) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66184735) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(3.0116144) q[0];
rz(3.1107483) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(2.470509) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6993461) q[0];
sx q[0];
rz(-1.5171577) q[0];
sx q[0];
rz(1.5357114) q[0];
x q[1];
rz(-3.1110711) q[2];
sx q[2];
rz(-1.3866716) q[2];
sx q[2];
rz(0.39779278) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56747251) q[1];
sx q[1];
rz(-1.6735055) q[1];
sx q[1];
rz(1.1979539) q[1];
x q[2];
rz(-3.0942261) q[3];
sx q[3];
rz(-2.2455375) q[3];
sx q[3];
rz(2.588152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.122763) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(-0.028586483) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(-1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(2.7291765) q[0];
rz(-1.6917797) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.9746045) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3827688) q[0];
sx q[0];
rz(-2.1817657) q[0];
sx q[0];
rz(0.2610892) q[0];
x q[1];
rz(0.27600482) q[2];
sx q[2];
rz(-2.2609684) q[2];
sx q[2];
rz(-1.3921757) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3762714) q[1];
sx q[1];
rz(-1.5404535) q[1];
sx q[1];
rz(-0.18860753) q[1];
x q[2];
rz(1.0409045) q[3];
sx q[3];
rz(-1.7003254) q[3];
sx q[3];
rz(1.1216175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2010487) q[2];
sx q[2];
rz(-0.87934914) q[2];
sx q[2];
rz(-0.18903014) q[2];
rz(0.14686251) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(-1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1259595) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(-2.2145859) q[0];
rz(1.3828297) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(1.6519201) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531567) q[0];
sx q[0];
rz(-2.2373767) q[0];
sx q[0];
rz(0.96320926) q[0];
rz(1.4719047) q[2];
sx q[2];
rz(-0.47404587) q[2];
sx q[2];
rz(-0.74741077) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8330194) q[1];
sx q[1];
rz(-0.25990572) q[1];
sx q[1];
rz(2.412917) q[1];
rz(-pi) q[2];
rz(-1.0350111) q[3];
sx q[3];
rz(-2.0097369) q[3];
sx q[3];
rz(-0.5700232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5513409) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(-1.7112973) q[2];
rz(-2.5643505) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(1.3841217) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5883314) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(0.28840315) q[0];
rz(-0.53238955) q[1];
sx q[1];
rz(-2.6817697) q[1];
sx q[1];
rz(0.14702252) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9592181) q[0];
sx q[0];
rz(-1.1897414) q[0];
sx q[0];
rz(-2.1931838) q[0];
x q[1];
rz(-0.23994259) q[2];
sx q[2];
rz(-1.502617) q[2];
sx q[2];
rz(-2.6477674) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30565572) q[1];
sx q[1];
rz(-2.5218997) q[1];
sx q[1];
rz(-1.5020919) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7266502) q[3];
sx q[3];
rz(-0.999756) q[3];
sx q[3];
rz(-1.6243638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0599351) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(-0.74404136) q[2];
rz(-0.75731164) q[3];
sx q[3];
rz(-1.363874) q[3];
sx q[3];
rz(-1.3967167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.025678) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(-0.81746447) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(-0.13851891) q[2];
sx q[2];
rz(-0.45590966) q[2];
sx q[2];
rz(1.6675303) q[2];
rz(2.1327303) q[3];
sx q[3];
rz(-1.4603793) q[3];
sx q[3];
rz(0.59752656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];