OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95076686) q[0];
sx q[0];
rz(1.289225) q[0];
sx q[0];
rz(9.3318648) q[0];
rz(3.0463123) q[1];
sx q[1];
rz(-2.4089101) q[1];
sx q[1];
rz(-1.9021775) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6042955) q[0];
sx q[0];
rz(-1.3193466) q[0];
sx q[0];
rz(-0.23668134) q[0];
rz(0.26387604) q[2];
sx q[2];
rz(-1.7603163) q[2];
sx q[2];
rz(2.5352728) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.8653633) q[1];
sx q[1];
rz(-1.5412093) q[1];
sx q[1];
rz(-0.33025708) q[1];
x q[2];
rz(0.57115023) q[3];
sx q[3];
rz(-1.1451654) q[3];
sx q[3];
rz(-1.7895123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4702845) q[2];
sx q[2];
rz(-1.4802063) q[2];
sx q[2];
rz(-1.7926463) q[2];
rz(0.74364439) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(2.9644137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0002366) q[0];
sx q[0];
rz(-1.5970705) q[0];
sx q[0];
rz(0.29362383) q[0];
rz(-2.1108421) q[1];
sx q[1];
rz(-1.3615969) q[1];
sx q[1];
rz(0.34293276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2156649) q[0];
sx q[0];
rz(-2.5747888) q[0];
sx q[0];
rz(1.180992) q[0];
x q[1];
rz(1.8937102) q[2];
sx q[2];
rz(-2.7512449) q[2];
sx q[2];
rz(0.54803145) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3174177) q[1];
sx q[1];
rz(-1.3567827) q[1];
sx q[1];
rz(2.432968) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9092122) q[3];
sx q[3];
rz(-0.65753257) q[3];
sx q[3];
rz(-0.38207182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0127516) q[2];
sx q[2];
rz(-1.3920471) q[2];
sx q[2];
rz(-2.3823605) q[2];
rz(-3.1208842) q[3];
sx q[3];
rz(-1.8755251) q[3];
sx q[3];
rz(0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52221209) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(-3.1371327) q[0];
rz(-2.4941173) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(-0.78537816) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68092184) q[0];
sx q[0];
rz(-2.9630532) q[0];
sx q[0];
rz(2.6723249) q[0];
x q[1];
rz(2.9663229) q[2];
sx q[2];
rz(-1.880155) q[2];
sx q[2];
rz(0.64495211) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8080146) q[1];
sx q[1];
rz(-1.6393264) q[1];
sx q[1];
rz(1.0295111) q[1];
x q[2];
rz(-2.4786199) q[3];
sx q[3];
rz(-2.1372652) q[3];
sx q[3];
rz(-2.278355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7711827) q[2];
sx q[2];
rz(-0.9202756) q[2];
sx q[2];
rz(1.3578337) q[2];
rz(0.34058288) q[3];
sx q[3];
rz(-0.95652306) q[3];
sx q[3];
rz(-2.3497439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99712813) q[0];
sx q[0];
rz(-1.0839394) q[0];
sx q[0];
rz(0.88515627) q[0];
rz(-2.3947233) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(-1.0964099) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089450739) q[0];
sx q[0];
rz(-0.47596395) q[0];
sx q[0];
rz(-2.3030998) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.030641067) q[2];
sx q[2];
rz(-1.5113792) q[2];
sx q[2];
rz(-0.63842809) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4064286) q[1];
sx q[1];
rz(-2.6044835) q[1];
sx q[1];
rz(-2.6178611) q[1];
rz(0.55620749) q[3];
sx q[3];
rz(-0.81327754) q[3];
sx q[3];
rz(0.63268703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5235644) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(-0.024624126) q[2];
rz(-0.63557449) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(-2.2583101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6898952) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(-0.46689335) q[0];
rz(-3.0310071) q[1];
sx q[1];
rz(-0.46375912) q[1];
sx q[1];
rz(2.0707524) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12884451) q[0];
sx q[0];
rz(-1.8086047) q[0];
sx q[0];
rz(-1.1046011) q[0];
x q[1];
rz(-2.6314965) q[2];
sx q[2];
rz(-0.3613216) q[2];
sx q[2];
rz(-2.9053743) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.22510281) q[1];
sx q[1];
rz(-1.9472633) q[1];
sx q[1];
rz(2.589499) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70160206) q[3];
sx q[3];
rz(-1.0402586) q[3];
sx q[3];
rz(-0.70487937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0232627) q[2];
sx q[2];
rz(-1.7795965) q[2];
sx q[2];
rz(-0.85232097) q[2];
rz(3.1039589) q[3];
sx q[3];
rz(-1.9084385) q[3];
sx q[3];
rz(-0.5948624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0090050176) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(1.8527385) q[0];
rz(1.6291078) q[1];
sx q[1];
rz(-2.0178724) q[1];
sx q[1];
rz(1.9992794) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2347533) q[0];
sx q[0];
rz(-1.2182353) q[0];
sx q[0];
rz(1.4666739) q[0];
rz(0.062131957) q[2];
sx q[2];
rz(-1.3634063) q[2];
sx q[2];
rz(2.1211989) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9276322) q[1];
sx q[1];
rz(-2.7134656) q[1];
sx q[1];
rz(2.5090748) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6352245) q[3];
sx q[3];
rz(-0.35249235) q[3];
sx q[3];
rz(1.5086205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(0.12040559) q[2];
rz(1.1558007) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8026546) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(1.4539723) q[0];
rz(1.5974207) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(-0.45101756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2149725) q[0];
sx q[0];
rz(-0.30719137) q[0];
sx q[0];
rz(0.19812576) q[0];
x q[1];
rz(-1.2043549) q[2];
sx q[2];
rz(-1.3533354) q[2];
sx q[2];
rz(1.5254453) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.800182) q[1];
sx q[1];
rz(-1.4497611) q[1];
sx q[1];
rz(2.7054663) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0093098442) q[3];
sx q[3];
rz(-2.3311989) q[3];
sx q[3];
rz(1.1610462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.94577998) q[2];
sx q[2];
rz(-1.7273629) q[2];
sx q[2];
rz(1.7808524) q[2];
rz(-2.1360548) q[3];
sx q[3];
rz(-1.1755627) q[3];
sx q[3];
rz(1.3895234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239531) q[0];
sx q[0];
rz(-0.12549505) q[0];
sx q[0];
rz(-0.70607591) q[0];
rz(-1.5977244) q[1];
sx q[1];
rz(-1.4223301) q[1];
sx q[1];
rz(0.85618883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5498915) q[0];
sx q[0];
rz(-0.020792637) q[0];
sx q[0];
rz(1.2077622) q[0];
rz(-1.1107619) q[2];
sx q[2];
rz(-2.3383814) q[2];
sx q[2];
rz(0.75424657) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44196196) q[1];
sx q[1];
rz(-0.61304997) q[1];
sx q[1];
rz(2.7837201) q[1];
rz(1.3230349) q[3];
sx q[3];
rz(-2.6399586) q[3];
sx q[3];
rz(-2.786123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84264821) q[2];
sx q[2];
rz(-1.7405258) q[2];
sx q[2];
rz(-2.974158) q[2];
rz(-1.5873448) q[3];
sx q[3];
rz(-0.7730248) q[3];
sx q[3];
rz(2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7635968) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(2.7698621) q[0];
rz(-1.4338214) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(0.30002123) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8512631) q[0];
sx q[0];
rz(-1.6870878) q[0];
sx q[0];
rz(2.854748) q[0];
x q[1];
rz(1.6113564) q[2];
sx q[2];
rz(-2.1475361) q[2];
sx q[2];
rz(-0.64740136) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91050657) q[1];
sx q[1];
rz(-0.66201895) q[1];
sx q[1];
rz(-0.20259095) q[1];
rz(1.7172708) q[3];
sx q[3];
rz(-2.4333242) q[3];
sx q[3];
rz(-2.012501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6649449) q[2];
sx q[2];
rz(-0.46755329) q[2];
sx q[2];
rz(-2.7139968) q[2];
rz(1.5852196) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(-0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72117358) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(0.46947259) q[0];
rz(2.2162614) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(-3.1184149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1184922) q[0];
sx q[0];
rz(-1.7351302) q[0];
sx q[0];
rz(0.53043764) q[0];
rz(0.59892861) q[2];
sx q[2];
rz(-0.99236503) q[2];
sx q[2];
rz(2.3430062) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35649037) q[1];
sx q[1];
rz(-1.39686) q[1];
sx q[1];
rz(-0.16696232) q[1];
rz(-0.84081991) q[3];
sx q[3];
rz(-1.429002) q[3];
sx q[3];
rz(-2.8805582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8699441) q[2];
sx q[2];
rz(-1.2056489) q[2];
sx q[2];
rz(0.49087697) q[2];
rz(1.0902181) q[3];
sx q[3];
rz(-2.1722983) q[3];
sx q[3];
rz(-0.31392613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.28521095) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(-0.27676997) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(-2.7748952) q[2];
sx q[2];
rz(-2.4187805) q[2];
sx q[2];
rz(0.4772966) q[2];
rz(-0.033944081) q[3];
sx q[3];
rz(-0.53474075) q[3];
sx q[3];
rz(-0.026215601) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
