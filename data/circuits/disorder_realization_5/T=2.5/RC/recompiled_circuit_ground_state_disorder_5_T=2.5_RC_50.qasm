OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.759999) q[0];
sx q[0];
rz(2.969709) q[0];
sx q[0];
rz(10.039731) q[0];
rz(1.3568658) q[1];
sx q[1];
rz(4.7525726) q[1];
sx q[1];
rz(9.5819028) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2113094) q[0];
sx q[0];
rz(-2.0727949) q[0];
sx q[0];
rz(-0.94663488) q[0];
rz(-0.31034361) q[2];
sx q[2];
rz(-2.9039798) q[2];
sx q[2];
rz(-1.0050736) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9188668) q[1];
sx q[1];
rz(-0.91321975) q[1];
sx q[1];
rz(0.01941733) q[1];
rz(-1.3376921) q[3];
sx q[3];
rz(-0.86458428) q[3];
sx q[3];
rz(-1.2561089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60363808) q[2];
sx q[2];
rz(-0.31542888) q[2];
sx q[2];
rz(-1.8559378) q[2];
rz(-1.8042608) q[3];
sx q[3];
rz(-0.95271102) q[3];
sx q[3];
rz(-3.0281236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86785698) q[0];
sx q[0];
rz(-1.8524167) q[0];
sx q[0];
rz(0.56647545) q[0];
rz(-1.4631588) q[1];
sx q[1];
rz(-0.17371829) q[1];
sx q[1];
rz(2.4991551) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3644407) q[0];
sx q[0];
rz(-2.6760347) q[0];
sx q[0];
rz(-1.6861395) q[0];
rz(-0.42936705) q[2];
sx q[2];
rz(-0.68898669) q[2];
sx q[2];
rz(-1.8935204) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3088544) q[1];
sx q[1];
rz(-0.843261) q[1];
sx q[1];
rz(0.067752167) q[1];
rz(-2.4107433) q[3];
sx q[3];
rz(-1.0789144) q[3];
sx q[3];
rz(-3.0036894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2598339) q[2];
sx q[2];
rz(-0.79462785) q[2];
sx q[2];
rz(2.3074522) q[2];
rz(0.086890876) q[3];
sx q[3];
rz(-2.0565242) q[3];
sx q[3];
rz(2.0104008) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6777545) q[0];
sx q[0];
rz(-1.0967655) q[0];
sx q[0];
rz(-2.0820397) q[0];
rz(0.73486596) q[1];
sx q[1];
rz(-3.1262472) q[1];
sx q[1];
rz(2.9954092) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8627166) q[0];
sx q[0];
rz(-1.4116316) q[0];
sx q[0];
rz(1.0907286) q[0];
rz(0.0038613316) q[2];
sx q[2];
rz(-1.6950102) q[2];
sx q[2];
rz(0.079380097) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9707036) q[1];
sx q[1];
rz(-1.7196989) q[1];
sx q[1];
rz(-1.7335684) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22821088) q[3];
sx q[3];
rz(-1.2973229) q[3];
sx q[3];
rz(0.21450689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.369578) q[2];
sx q[2];
rz(-0.39083189) q[2];
sx q[2];
rz(2.3005627) q[2];
rz(3.0817025) q[3];
sx q[3];
rz(-1.4845029) q[3];
sx q[3];
rz(0.64391518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.197072) q[0];
sx q[0];
rz(-0.50839013) q[0];
sx q[0];
rz(-0.12666853) q[0];
rz(-1.4177167) q[1];
sx q[1];
rz(-0.020514943) q[1];
sx q[1];
rz(0.98974481) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5602172) q[0];
sx q[0];
rz(-2.2994083) q[0];
sx q[0];
rz(-1.5422264) q[0];
rz(-pi) q[1];
rz(1.1737972) q[2];
sx q[2];
rz(-1.4922895) q[2];
sx q[2];
rz(2.5369413) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.77074487) q[1];
sx q[1];
rz(-2.641692) q[1];
sx q[1];
rz(-1.1558394) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5775388) q[3];
sx q[3];
rz(-0.28279916) q[3];
sx q[3];
rz(1.1081616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7030455) q[2];
sx q[2];
rz(-2.788471) q[2];
sx q[2];
rz(2.9959196) q[2];
rz(0.4438256) q[3];
sx q[3];
rz(-1.5747109) q[3];
sx q[3];
rz(-1.1183967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038427453) q[0];
sx q[0];
rz(-1.9157836) q[0];
sx q[0];
rz(-2.3577754) q[0];
rz(-1.8812284) q[1];
sx q[1];
rz(-0.0030219373) q[1];
sx q[1];
rz(-2.9154215) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0167671) q[0];
sx q[0];
rz(-1.1065605) q[0];
sx q[0];
rz(2.8874257) q[0];
rz(-pi) q[1];
rz(1.1532835) q[2];
sx q[2];
rz(-1.8476356) q[2];
sx q[2];
rz(-1.4275488) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60701398) q[1];
sx q[1];
rz(-1.6456938) q[1];
sx q[1];
rz(-2.043488) q[1];
rz(-pi) q[2];
x q[2];
rz(2.466684) q[3];
sx q[3];
rz(-2.2332472) q[3];
sx q[3];
rz(-1.1250594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0164531) q[2];
sx q[2];
rz(-0.21802248) q[2];
sx q[2];
rz(1.7832635) q[2];
rz(-0.45768091) q[3];
sx q[3];
rz(-1.4384559) q[3];
sx q[3];
rz(-2.5911736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.054475527) q[0];
sx q[0];
rz(-3.1396907) q[0];
sx q[0];
rz(-3.0892293) q[0];
rz(-0.26854435) q[1];
sx q[1];
rz(-1.9895376) q[1];
sx q[1];
rz(2.9103738) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5729882) q[0];
sx q[0];
rz(-1.3954961) q[0];
sx q[0];
rz(0.25956395) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74502142) q[2];
sx q[2];
rz(-1.3437437) q[2];
sx q[2];
rz(-2.678427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0760804) q[1];
sx q[1];
rz(-1.8022746) q[1];
sx q[1];
rz(0.3409909) q[1];
rz(2.1778132) q[3];
sx q[3];
rz(-1.7208817) q[3];
sx q[3];
rz(-0.039473783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82983595) q[2];
sx q[2];
rz(-0.51218963) q[2];
sx q[2];
rz(2.3066985) q[2];
rz(-0.087415047) q[3];
sx q[3];
rz(-2.0525457) q[3];
sx q[3];
rz(-2.1277229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3111303) q[0];
sx q[0];
rz(-2.1489547) q[0];
sx q[0];
rz(2.2845238) q[0];
rz(2.3509707) q[1];
sx q[1];
rz(-3.1309083) q[1];
sx q[1];
rz(-2.0649921) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.650855) q[0];
sx q[0];
rz(-1.8371353) q[0];
sx q[0];
rz(-1.8835212) q[0];
rz(-pi) q[1];
rz(2.358165) q[2];
sx q[2];
rz(-1.9840709) q[2];
sx q[2];
rz(0.30488067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69546135) q[1];
sx q[1];
rz(-0.82099719) q[1];
sx q[1];
rz(1.9638792) q[1];
rz(-pi) q[2];
rz(-0.3420426) q[3];
sx q[3];
rz(-1.6989467) q[3];
sx q[3];
rz(-0.32610937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2662346) q[2];
sx q[2];
rz(-0.42576063) q[2];
sx q[2];
rz(2.2597964) q[2];
rz(-1.4074696) q[3];
sx q[3];
rz(-2.2018645) q[3];
sx q[3];
rz(2.5772429) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8565916) q[0];
sx q[0];
rz(-3.1398616) q[0];
sx q[0];
rz(0.28975394) q[0];
rz(1.9081135) q[1];
sx q[1];
rz(-1.8461485) q[1];
sx q[1];
rz(-0.27898702) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4623337) q[0];
sx q[0];
rz(-1.7391608) q[0];
sx q[0];
rz(3.116045) q[0];
rz(-pi) q[1];
rz(-0.60207587) q[2];
sx q[2];
rz(-0.47594949) q[2];
sx q[2];
rz(2.0747831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28767868) q[1];
sx q[1];
rz(-0.33386896) q[1];
sx q[1];
rz(-1.3592657) q[1];
rz(-2.2905228) q[3];
sx q[3];
rz(-1.2382231) q[3];
sx q[3];
rz(0.50902396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4616619) q[2];
sx q[2];
rz(-1.8671904) q[2];
sx q[2];
rz(-2.8604841) q[2];
rz(-0.57349652) q[3];
sx q[3];
rz(-2.1888013) q[3];
sx q[3];
rz(-0.41962418) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9004274) q[0];
sx q[0];
rz(-0.92210162) q[0];
sx q[0];
rz(1.4379372) q[0];
rz(2.9519713) q[1];
sx q[1];
rz(-0.01556839) q[1];
sx q[1];
rz(-1.9059034) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8164809) q[0];
sx q[0];
rz(-2.6345523) q[0];
sx q[0];
rz(1.965824) q[0];
rz(-pi) q[1];
rz(-0.81053712) q[2];
sx q[2];
rz(-0.97522465) q[2];
sx q[2];
rz(-0.11071225) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37285629) q[1];
sx q[1];
rz(-2.8017178) q[1];
sx q[1];
rz(1.0699141) q[1];
x q[2];
rz(0.017589324) q[3];
sx q[3];
rz(-1.5776199) q[3];
sx q[3];
rz(-1.9987035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8427061) q[2];
sx q[2];
rz(-1.3624531) q[2];
sx q[2];
rz(-2.5617808) q[2];
rz(1.7719571) q[3];
sx q[3];
rz(-2.7176791) q[3];
sx q[3];
rz(-0.4637318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6746826) q[0];
sx q[0];
rz(-1.5002102) q[0];
sx q[0];
rz(1.4379733) q[0];
rz(-1.9101608) q[1];
sx q[1];
rz(-0.91130251) q[1];
sx q[1];
rz(1.4968754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5176277) q[0];
sx q[0];
rz(-1.3056439) q[0];
sx q[0];
rz(0.82268663) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6193498) q[2];
sx q[2];
rz(-1.2127611) q[2];
sx q[2];
rz(-2.4286662) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0857049) q[1];
sx q[1];
rz(-3.1236557) q[1];
sx q[1];
rz(-2.9824801) q[1];
rz(-1.0375828) q[3];
sx q[3];
rz(-0.54736558) q[3];
sx q[3];
rz(-1.411388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15184312) q[2];
sx q[2];
rz(-1.5785549) q[2];
sx q[2];
rz(2.6452276) q[2];
rz(2.1967891) q[3];
sx q[3];
rz(-0.0050408575) q[3];
sx q[3];
rz(2.4387824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6481358) q[0];
sx q[0];
rz(-1.428816) q[0];
sx q[0];
rz(-1.1270123) q[0];
rz(1.5588749) q[1];
sx q[1];
rz(-0.3180779) q[1];
sx q[1];
rz(0.19837468) q[1];
rz(1.3140903) q[2];
sx q[2];
rz(-1.5460925) q[2];
sx q[2];
rz(1.8800541) q[2];
rz(-2.5924223) q[3];
sx q[3];
rz(-2.5097889) q[3];
sx q[3];
rz(-0.24181152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
