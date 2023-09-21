OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(-2.6383658) q[0];
sx q[0];
rz(0.72416645) q[0];
rz(6.9231482) q[1];
sx q[1];
rz(5.7531113) q[1];
sx q[1];
rz(2.35676) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35207507) q[0];
sx q[0];
rz(-1.2791469) q[0];
sx q[0];
rz(1.350499) q[0];
rz(-pi) q[1];
rz(0.41854026) q[2];
sx q[2];
rz(-1.4782018) q[2];
sx q[2];
rz(1.4946403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5074871) q[1];
sx q[1];
rz(-1.3480942) q[1];
sx q[1];
rz(2.1102064) q[1];
rz(-pi) q[2];
rz(-3.1078994) q[3];
sx q[3];
rz(-1.9851306) q[3];
sx q[3];
rz(-1.8309483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(0.067967728) q[2];
rz(-3.0170278) q[3];
sx q[3];
rz(-2.8187276) q[3];
sx q[3];
rz(-1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92000604) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(-0.24366972) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.4032204) q[1];
sx q[1];
rz(-1.3557281) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5403554) q[0];
sx q[0];
rz(-1.3155126) q[0];
sx q[0];
rz(0.028098696) q[0];
rz(0.83265702) q[2];
sx q[2];
rz(-2.6359733) q[2];
sx q[2];
rz(-2.6308699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4888549) q[1];
sx q[1];
rz(-1.4421717) q[1];
sx q[1];
rz(1.1210404) q[1];
x q[2];
rz(-1.5242819) q[3];
sx q[3];
rz(-0.97913137) q[3];
sx q[3];
rz(-3.068963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0624861) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(-0.24965723) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24519414) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(0.89843345) q[0];
rz(-1.3348745) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(1.867884) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0262336) q[0];
sx q[0];
rz(-2.6991978) q[0];
sx q[0];
rz(-2.1917403) q[0];
rz(-0.31295915) q[2];
sx q[2];
rz(-2.8319781) q[2];
sx q[2];
rz(-1.5086053) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
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
rz(pi/2) q[1];
sx q[1];
rz(1.0597824) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(-0.95345062) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8614486) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(-0.24818534) q[0];
rz(-2.10363) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(-0.074137069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5384597) q[0];
sx q[0];
rz(-2.5939301) q[0];
sx q[0];
rz(-2.0752226) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.147869) q[2];
sx q[2];
rz(-2.4175156) q[2];
sx q[2];
rz(0.99087447) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76064199) q[1];
sx q[1];
rz(-1.4529072) q[1];
sx q[1];
rz(-1.2537434) q[1];
x q[2];
rz(-1.3965963) q[3];
sx q[3];
rz(-1.8074236) q[3];
sx q[3];
rz(-0.06121204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8903824) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-3.1029491) q[2];
rz(-0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(-0.27004778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.5357635) q[0];
sx q[0];
rz(1.779153) q[0];
rz(2.3249987) q[1];
sx q[1];
rz(-1.8530308) q[1];
sx q[1];
rz(-1.978925) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029318854) q[0];
sx q[0];
rz(-0.11538878) q[0];
sx q[0];
rz(-1.3089887) q[0];
rz(-pi) q[1];
rz(-1.0983724) q[2];
sx q[2];
rz(-2.782151) q[2];
sx q[2];
rz(-0.098509468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77046493) q[1];
sx q[1];
rz(-0.41509291) q[1];
sx q[1];
rz(-1.0789372) q[1];
rz(2.707162) q[3];
sx q[3];
rz(-2.0475004) q[3];
sx q[3];
rz(-1.7720951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6267307) q[2];
sx q[2];
rz(-1.0802439) q[2];
sx q[2];
rz(-2.999372) q[2];
rz(-0.90406117) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(-0.18946762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11480039) q[0];
sx q[0];
rz(-0.27652201) q[0];
sx q[0];
rz(-1.6739155) q[0];
rz(-0.57178512) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(2.8335559) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338617) q[0];
sx q[0];
rz(-2.9836285) q[0];
sx q[0];
rz(0.64361848) q[0];
rz(-pi) q[1];
rz(-2.5295837) q[2];
sx q[2];
rz(-1.0527805) q[2];
sx q[2];
rz(-0.81105622) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.088965485) q[1];
sx q[1];
rz(-1.7672156) q[1];
sx q[1];
rz(-2.2599225) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4792535) q[3];
sx q[3];
rz(-1.9178101) q[3];
sx q[3];
rz(-0.43859827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3623111) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(-1.1479088) q[2];
rz(0.71427304) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(-2.4961297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66184735) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(3.0116144) q[0];
rz(0.030844363) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(-0.67108363) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13043159) q[0];
sx q[0];
rz(-1.5357619) q[0];
sx q[0];
rz(0.053671562) q[0];
rz(-pi) q[1];
rz(1.7332156) q[2];
sx q[2];
rz(-0.186609) q[2];
sx q[2];
rz(-0.23255541) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.259686) q[1];
sx q[1];
rz(-2.7555008) q[1];
sx q[1];
rz(-1.8465471) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5116793) q[3];
sx q[3];
rz(-2.4654508) q[3];
sx q[3];
rz(2.6638871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0188296) q[2];
sx q[2];
rz(-1.5315703) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(-0.028586483) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1639444) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(-2.7291765) q[0];
rz(-1.6917797) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.9746045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7588239) q[0];
sx q[0];
rz(-0.95982691) q[0];
sx q[0];
rz(-2.8805034) q[0];
x q[1];
rz(-0.86161676) q[2];
sx q[2];
rz(-1.3590727) q[2];
sx q[2];
rz(0.35702969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1046909) q[1];
sx q[1];
rz(-2.9505886) q[1];
sx q[1];
rz(-2.9810993) q[1];
rz(-pi) q[2];
rz(-2.9917631) q[3];
sx q[3];
rz(-2.095788) q[3];
sx q[3];
rz(-0.52469745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-0.18903014) q[2];
rz(0.14686251) q[3];
sx q[3];
rz(-0.18460128) q[3];
sx q[3];
rz(1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(0.92700672) q[0];
rz(1.3828297) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(-1.6519201) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2531567) q[0];
sx q[0];
rz(-2.2373767) q[0];
sx q[0];
rz(-0.96320926) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4719047) q[2];
sx q[2];
rz(-2.6675468) q[2];
sx q[2];
rz(0.74741077) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3085732) q[1];
sx q[1];
rz(-0.25990572) q[1];
sx q[1];
rz(0.72867568) q[1];
rz(-pi) q[2];
rz(0.49976607) q[3];
sx q[3];
rz(-1.0904113) q[3];
sx q[3];
rz(-1.8936611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5513409) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(1.4302953) q[2];
rz(-2.5643505) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(-1.757471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5883314) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(-2.8531895) q[0];
rz(-2.6092031) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(-2.9945701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4924016) q[0];
sx q[0];
rz(-2.1426139) q[0];
sx q[0];
rz(2.6834821) q[0];
x q[1];
rz(0.23994259) q[2];
sx q[2];
rz(-1.6389756) q[2];
sx q[2];
rz(0.49382526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22132561) q[1];
sx q[1];
rz(-2.1888071) q[1];
sx q[1];
rz(3.0926535) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7266502) q[3];
sx q[3];
rz(-2.1418367) q[3];
sx q[3];
rz(1.5172289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0816575) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(-0.74404136) q[2];
rz(-0.75731164) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(-1.7448759) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.025678) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(-0.81746447) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(3.0030737) q[2];
sx q[2];
rz(-0.45590966) q[2];
sx q[2];
rz(1.6675303) q[2];
rz(1.3656473) q[3];
sx q[3];
rz(-2.57006) q[3];
sx q[3];
rz(-0.80001696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
