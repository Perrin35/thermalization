OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4164299) q[0];
sx q[0];
rz(-0.13983146) q[0];
sx q[0];
rz(-2.5319985) q[0];
rz(-2.4729589) q[1];
sx q[1];
rz(-0.86548391) q[1];
sx q[1];
rz(-3.0545711) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7323157) q[0];
sx q[0];
rz(-1.9415932) q[0];
sx q[0];
rz(2.2229574) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4423223) q[2];
sx q[2];
rz(-0.80768425) q[2];
sx q[2];
rz(1.3070004) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.062292369) q[1];
sx q[1];
rz(-1.0907409) q[1];
sx q[1];
rz(2.9855707) q[1];
rz(-pi) q[2];
rz(0.41214715) q[3];
sx q[3];
rz(-1.5282359) q[3];
sx q[3];
rz(1.7736848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13880754) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(-0.37386093) q[2];
rz(-2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(-2.9132304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574629) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.194838) q[0];
rz(-3.0589814) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(3.1412178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9449687) q[0];
sx q[0];
rz(-1.225291) q[0];
sx q[0];
rz(-2.0161122) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8325009) q[2];
sx q[2];
rz(-1.3771025) q[2];
sx q[2];
rz(-1.2145834) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2692711) q[1];
sx q[1];
rz(-1.806013) q[1];
sx q[1];
rz(-0.61551952) q[1];
x q[2];
rz(-2.681148) q[3];
sx q[3];
rz(-1.6201971) q[3];
sx q[3];
rz(1.7668262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3327545) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(2.9648798) q[2];
rz(-0.79408944) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394543) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(-0.40476558) q[0];
rz(1.8602712) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.3084897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8264016) q[0];
sx q[0];
rz(-2.0668525) q[0];
sx q[0];
rz(1.0415003) q[0];
rz(-pi) q[1];
rz(-0.53738014) q[2];
sx q[2];
rz(-1.3668622) q[2];
sx q[2];
rz(-2.1539719) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.96924671) q[1];
sx q[1];
rz(-1.3174651) q[1];
sx q[1];
rz(-2.2971056) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4806467) q[3];
sx q[3];
rz(-1.1336859) q[3];
sx q[3];
rz(-2.3428832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9222766) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(-3.1211839) q[2];
rz(-1.071788) q[3];
sx q[3];
rz(-1.1276378) q[3];
sx q[3];
rz(0.66334692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(0.91598696) q[0];
rz(-0.46332106) q[1];
sx q[1];
rz(-1.0659734) q[1];
sx q[1];
rz(2.0844918) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0650478) q[0];
sx q[0];
rz(-1.1335982) q[0];
sx q[0];
rz(1.5989499) q[0];
rz(-pi) q[1];
rz(0.82614233) q[2];
sx q[2];
rz(-2.38378) q[2];
sx q[2];
rz(0.10105029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.033356655) q[1];
sx q[1];
rz(-1.557096) q[1];
sx q[1];
rz(1.800888) q[1];
x q[2];
rz(2.9017157) q[3];
sx q[3];
rz(-0.7719709) q[3];
sx q[3];
rz(-0.039226942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2477734) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(2.7988953) q[2];
rz(1.6977067) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(2.4263583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3559568) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(-3.0551531) q[0];
rz(-1.7516288) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(2.6729029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5422573) q[0];
sx q[0];
rz(-2.242803) q[0];
sx q[0];
rz(0.50962781) q[0];
rz(-pi) q[1];
rz(-1.7519978) q[2];
sx q[2];
rz(-1.8743519) q[2];
sx q[2];
rz(0.62523491) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7680109) q[1];
sx q[1];
rz(-1.7005159) q[1];
sx q[1];
rz(2.9380161) q[1];
rz(-pi) q[2];
rz(-0.080805578) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.1420574) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(-2.2772677) q[2];
rz(-0.21720973) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7140759) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(-0.19700225) q[0];
rz(-1.3621832) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(2.8053455) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2339736) q[0];
sx q[0];
rz(-0.95171463) q[0];
sx q[0];
rz(-1.0236077) q[0];
x q[1];
rz(2.5480812) q[2];
sx q[2];
rz(-1.0700018) q[2];
sx q[2];
rz(1.3154495) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0834004) q[1];
sx q[1];
rz(-0.5914878) q[1];
sx q[1];
rz(-3.0565492) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69643173) q[3];
sx q[3];
rz(-1.7034354) q[3];
sx q[3];
rz(-1.2030676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(0.84632787) q[2];
rz(-1.9019295) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(-3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3939312) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(0.25892648) q[0];
rz(-1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-2.0475725) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0897113) q[0];
sx q[0];
rz(-1.9117172) q[0];
sx q[0];
rz(2.0345576) q[0];
x q[1];
rz(-1.3482434) q[2];
sx q[2];
rz(-1.5828653) q[2];
sx q[2];
rz(1.0782858) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6106231) q[1];
sx q[1];
rz(-1.5986643) q[1];
sx q[1];
rz(1.665297) q[1];
rz(-1.5488946) q[3];
sx q[3];
rz(-1.4727482) q[3];
sx q[3];
rz(2.736562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9672433) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(-2.8209177) q[2];
rz(-2.705412) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(-2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86673474) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(-2.136769) q[0];
rz(-2.5514305) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(0.80387962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.852254) q[0];
sx q[0];
rz(-0.58710557) q[0];
sx q[0];
rz(3.0703074) q[0];
rz(0.55401037) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(-2.6487034) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.67919532) q[1];
sx q[1];
rz(-3.052366) q[1];
sx q[1];
rz(1.6855082) q[1];
rz(-0.6635267) q[3];
sx q[3];
rz(-1.6420206) q[3];
sx q[3];
rz(-1.0458667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8089495) q[2];
sx q[2];
rz(-2.2968764) q[2];
sx q[2];
rz(-2.1311029) q[2];
rz(0.60339749) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.5378961) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(0.4075152) q[0];
rz(-2.852476) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(-0.75072748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6347046) q[0];
sx q[0];
rz(-1.2872211) q[0];
sx q[0];
rz(0.30446913) q[0];
x q[1];
rz(2.8640792) q[2];
sx q[2];
rz(-0.29007426) q[2];
sx q[2];
rz(-1.2558503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99414413) q[1];
sx q[1];
rz(-1.3683043) q[1];
sx q[1];
rz(1.2575498) q[1];
rz(-0.084089355) q[3];
sx q[3];
rz(-0.30232271) q[3];
sx q[3];
rz(-1.7526527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0728545) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(1.9753974) q[2];
rz(1.3646305) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(-3.1197746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5584548) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(1.3903842) q[0];
rz(-2.8109) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.6814544) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9598823) q[0];
sx q[0];
rz(-1.0707756) q[0];
sx q[0];
rz(-0.48258968) q[0];
rz(-0.023256217) q[2];
sx q[2];
rz(-2.0746982) q[2];
sx q[2];
rz(2.0397253) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32669386) q[1];
sx q[1];
rz(-1.1732475) q[1];
sx q[1];
rz(1.3709929) q[1];
rz(-2.4246033) q[3];
sx q[3];
rz(-0.90962142) q[3];
sx q[3];
rz(-2.9644074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.34974393) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(-1.6798518) q[2];
rz(1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(-2.6859443) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5158952) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(1.760578) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(1.21576) q[2];
sx q[2];
rz(-2.3264865) q[2];
sx q[2];
rz(-2.9441499) q[2];
rz(0.47386668) q[3];
sx q[3];
rz(-2.4110473) q[3];
sx q[3];
rz(-1.9163781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
