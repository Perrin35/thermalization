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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7895176) q[0];
sx q[0];
rz(-1.8624458) q[0];
sx q[0];
rz(1.7910936) q[0];
x q[1];
rz(1.6720812) q[2];
sx q[2];
rz(-1.1541608) q[2];
sx q[2];
rz(0.11726221) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29014978) q[1];
sx q[1];
rz(-0.57934299) q[1];
sx q[1];
rz(1.9860553) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1562528) q[3];
sx q[3];
rz(-1.539955) q[3];
sx q[3];
rz(-2.8950092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.589754) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-0.067967728) q[2];
rz(3.0170278) q[3];
sx q[3];
rz(-2.8187276) q[3];
sx q[3];
rz(1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.92000604) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(2.8979229) q[0];
rz(0.63175732) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(1.3557281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5403554) q[0];
sx q[0];
rz(-1.8260801) q[0];
sx q[0];
rz(-3.113494) q[0];
rz(-0.35661125) q[2];
sx q[2];
rz(-1.2043673) q[2];
sx q[2];
rz(-0.29417843) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7999468) q[1];
sx q[1];
rz(-0.46657545) q[1];
sx q[1];
rz(-1.2816309) q[1];
x q[2];
rz(1.5242819) q[3];
sx q[3];
rz(-2.1624613) q[3];
sx q[3];
rz(0.072629645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(-0.24965723) q[2];
rz(-0.50659531) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(-2.8095968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24519414) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(2.2431592) q[0];
rz(1.3348745) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(-1.2737087) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0231409) q[0];
sx q[0];
rz(-1.8225192) q[0];
sx q[0];
rz(1.9385507) q[0];
rz(-pi) q[1];
rz(-1.6689698) q[2];
sx q[2];
rz(-1.2766826) q[2];
sx q[2];
rz(1.3054747) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0395567) q[1];
sx q[1];
rz(-1.2219056) q[1];
sx q[1];
rz(0.7015014) q[1];
rz(-pi) q[2];
rz(-2.0687194) q[3];
sx q[3];
rz(-2.2306799) q[3];
sx q[3];
rz(-0.68237309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0818103) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(2.188142) q[2];
rz(-3.1070784) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(-0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2801441) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(-2.8934073) q[0];
rz(-2.10363) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(-0.074137069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1124681) q[0];
sx q[0];
rz(-1.0974786) q[0];
sx q[0];
rz(-2.8549457) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2494227) q[2];
sx q[2];
rz(-1.2954419) q[2];
sx q[2];
rz(-2.8868669) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46596913) q[1];
sx q[1];
rz(-2.8040261) q[1];
sx q[1];
rz(1.9338495) q[1];
rz(1.3965963) q[3];
sx q[3];
rz(-1.3341691) q[3];
sx q[3];
rz(-0.06121204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8903824) q[2];
sx q[2];
rz(-0.40428287) q[2];
sx q[2];
rz(3.1029491) q[2];
rz(2.1679227) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(-2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(-1.779153) q[0];
rz(0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(-1.978925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8602596) q[0];
sx q[0];
rz(-1.6006002) q[0];
sx q[0];
rz(-1.4593065) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0432203) q[2];
sx q[2];
rz(-2.782151) q[2];
sx q[2];
rz(-3.0430832) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8414383) q[1];
sx q[1];
rz(-1.9342124) q[1];
sx q[1];
rz(-2.9363948) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.707162) q[3];
sx q[3];
rz(-1.0940922) q[3];
sx q[3];
rz(1.3694976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5148619) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(2.999372) q[2];
rz(2.2375315) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(-2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11480039) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.4676771) q[0];
rz(-0.57178512) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(-2.8335559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7361569) q[0];
sx q[0];
rz(-1.476256) q[0];
sx q[0];
rz(-0.12673881) q[0];
rz(-0.96254827) q[2];
sx q[2];
rz(-2.0934009) q[2];
sx q[2];
rz(-2.7163497) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8923556) q[1];
sx q[1];
rz(-0.71214572) q[1];
sx q[1];
rz(-1.874079) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8940053) q[3];
sx q[3];
rz(-0.35841225) q[3];
sx q[3];
rz(0.7022411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77928153) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(1.9936838) q[2];
rz(-0.71427304) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66184735) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(-0.1299783) q[0];
rz(3.1107483) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(2.470509) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6993461) q[0];
sx q[0];
rz(-1.5171577) q[0];
sx q[0];
rz(-1.6058812) q[0];
x q[1];
rz(1.408377) q[2];
sx q[2];
rz(-2.9549837) q[2];
sx q[2];
rz(2.9090372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8819067) q[1];
sx q[1];
rz(-2.7555008) q[1];
sx q[1];
rz(-1.2950456) q[1];
x q[2];
rz(3.0942261) q[3];
sx q[3];
rz(-0.8960552) q[3];
sx q[3];
rz(-0.55344068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0188296) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(2.5637131) q[2];
rz(3.1130062) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(-1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9776483) q[0];
sx q[0];
rz(-2.2387235) q[0];
sx q[0];
rz(-0.41241616) q[0];
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
rz(0.8186504) q[0];
sx q[0];
rz(-0.65781051) q[0];
sx q[0];
rz(-1.2176745) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8655878) q[2];
sx q[2];
rz(-2.2609684) q[2];
sx q[2];
rz(1.3921757) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.941277) q[1];
sx q[1];
rz(-1.759316) q[1];
sx q[1];
rz(1.6016866) q[1];
x q[2];
rz(-1.8230209) q[3];
sx q[3];
rz(-0.54402292) q[3];
sx q[3];
rz(0.23214425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2010487) q[2];
sx q[2];
rz(-0.87934914) q[2];
sx q[2];
rz(2.9525625) q[2];
rz(2.9947301) q[3];
sx q[3];
rz(-0.18460128) q[3];
sx q[3];
rz(1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(2.2145859) q[0];
rz(-1.758763) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(-1.4896726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2531567) q[0];
sx q[0];
rz(-2.2373767) q[0];
sx q[0];
rz(2.1783834) q[0];
rz(-pi) q[1];
rz(-1.4719047) q[2];
sx q[2];
rz(-2.6675468) q[2];
sx q[2];
rz(2.3941819) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8330194) q[1];
sx q[1];
rz(-2.8816869) q[1];
sx q[1];
rz(0.72867568) q[1];
rz(-pi) q[2];
rz(-2.1065815) q[3];
sx q[3];
rz(-1.1318558) q[3];
sx q[3];
rz(-0.5700232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5902517) q[2];
sx q[2];
rz(-0.44289032) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5532613) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(-0.28840315) q[0];
rz(2.6092031) q[1];
sx q[1];
rz(-2.6817697) q[1];
sx q[1];
rz(-2.9945701) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9592181) q[0];
sx q[0];
rz(-1.1897414) q[0];
sx q[0];
rz(2.1931838) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23994259) q[2];
sx q[2];
rz(-1.502617) q[2];
sx q[2];
rz(0.49382526) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.920267) q[1];
sx q[1];
rz(-0.95278554) q[1];
sx q[1];
rz(-3.0926535) q[1];
rz(-pi) q[2];
rz(1.4149425) q[3];
sx q[3];
rz(-0.999756) q[3];
sx q[3];
rz(-1.6243638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0816575) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(-2.3975513) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(1.1159146) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(0.81746447) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(-3.0030737) q[2];
sx q[2];
rz(-2.685683) q[2];
sx q[2];
rz(-1.4740623) q[2];
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