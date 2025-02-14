OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2616413) q[0];
sx q[0];
rz(-0.14544848) q[0];
sx q[0];
rz(0.26785904) q[0];
rz(1.574006) q[1];
sx q[1];
rz(-2.9730453) q[1];
sx q[1];
rz(2.5634917) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12500873) q[0];
sx q[0];
rz(-1.8468231) q[0];
sx q[0];
rz(2.2651423) q[0];
x q[1];
rz(-0.16139754) q[2];
sx q[2];
rz(-2.4366852) q[2];
sx q[2];
rz(-2.7053331) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0511623) q[1];
sx q[1];
rz(-1.9774861) q[1];
sx q[1];
rz(2.6654408) q[1];
x q[2];
rz(-1.3760482) q[3];
sx q[3];
rz(-0.71310242) q[3];
sx q[3];
rz(-2.3793067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1987004) q[2];
sx q[2];
rz(-2.7304724) q[2];
sx q[2];
rz(1.6655507) q[2];
rz(0.57925159) q[3];
sx q[3];
rz(-1.9871291) q[3];
sx q[3];
rz(2.4756685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2838374) q[0];
sx q[0];
rz(-1.0455766) q[0];
sx q[0];
rz(0.055140821) q[0];
rz(0.91446963) q[1];
sx q[1];
rz(-1.6998467) q[1];
sx q[1];
rz(1.9257911) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037579868) q[0];
sx q[0];
rz(-1.5079588) q[0];
sx q[0];
rz(1.5893905) q[0];
rz(-pi) q[1];
rz(0.74389768) q[2];
sx q[2];
rz(-0.66989952) q[2];
sx q[2];
rz(3.0238341) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2760073) q[1];
sx q[1];
rz(-2.7424208) q[1];
sx q[1];
rz(0.047476032) q[1];
rz(-1.7082105) q[3];
sx q[3];
rz(-2.0169037) q[3];
sx q[3];
rz(2.844939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9332283) q[2];
sx q[2];
rz(-2.6231982) q[2];
sx q[2];
rz(2.5210157) q[2];
rz(1.4536475) q[3];
sx q[3];
rz(-2.1098638) q[3];
sx q[3];
rz(1.0769963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2700014) q[0];
sx q[0];
rz(-0.33137614) q[0];
sx q[0];
rz(1.510386) q[0];
rz(2.467678) q[1];
sx q[1];
rz(-1.8464512) q[1];
sx q[1];
rz(-1.6253701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0533901) q[0];
sx q[0];
rz(-0.53117472) q[0];
sx q[0];
rz(-1.6199528) q[0];
rz(-0.38336945) q[2];
sx q[2];
rz(-0.96770937) q[2];
sx q[2];
rz(2.4281339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.060884692) q[1];
sx q[1];
rz(-1.1509893) q[1];
sx q[1];
rz(1.3221198) q[1];
rz(-pi) q[2];
rz(0.45511873) q[3];
sx q[3];
rz(-0.60448217) q[3];
sx q[3];
rz(-0.68463078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7565833) q[2];
sx q[2];
rz(-1.5998806) q[2];
sx q[2];
rz(2.4760683) q[2];
rz(-2.3982128) q[3];
sx q[3];
rz(-2.0205108) q[3];
sx q[3];
rz(-0.22751787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7706364) q[0];
sx q[0];
rz(-2.0891068) q[0];
sx q[0];
rz(-1.8413405) q[0];
rz(-3.0216253) q[1];
sx q[1];
rz(-2.0610466) q[1];
sx q[1];
rz(1.0467451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9366953) q[0];
sx q[0];
rz(-3.129382) q[0];
sx q[0];
rz(1.3365251) q[0];
x q[1];
rz(-0.36645269) q[2];
sx q[2];
rz(-2.2546446) q[2];
sx q[2];
rz(-0.89095913) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0485301) q[1];
sx q[1];
rz(-2.0566062) q[1];
sx q[1];
rz(0.062315224) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98174121) q[3];
sx q[3];
rz(-1.9509754) q[3];
sx q[3];
rz(-1.5693992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7946502) q[2];
sx q[2];
rz(-2.0796937) q[2];
sx q[2];
rz(0.77155716) q[2];
rz(2.0906585) q[3];
sx q[3];
rz(-1.0335048) q[3];
sx q[3];
rz(-1.9482025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0070888) q[0];
sx q[0];
rz(-1.6935885) q[0];
sx q[0];
rz(0.80528468) q[0];
rz(1.443642) q[1];
sx q[1];
rz(-2.3256358) q[1];
sx q[1];
rz(3.1051292) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48591025) q[0];
sx q[0];
rz(-0.56942372) q[0];
sx q[0];
rz(1.6605366) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1389937) q[2];
sx q[2];
rz(-1.6284962) q[2];
sx q[2];
rz(-0.94896832) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.441022) q[1];
sx q[1];
rz(-2.6614958) q[1];
sx q[1];
rz(-1.6861077) q[1];
x q[2];
rz(2.8920435) q[3];
sx q[3];
rz(-0.59246906) q[3];
sx q[3];
rz(-1.0147926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61520758) q[2];
sx q[2];
rz(-2.0168763) q[2];
sx q[2];
rz(0.57662326) q[2];
rz(-1.5455101) q[3];
sx q[3];
rz(-2.2871064) q[3];
sx q[3];
rz(2.5647148) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1103519) q[0];
sx q[0];
rz(-1.4875655) q[0];
sx q[0];
rz(-2.5417969) q[0];
rz(0.80232969) q[1];
sx q[1];
rz(-1.0917412) q[1];
sx q[1];
rz(-0.64782992) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9947056) q[0];
sx q[0];
rz(-2.0118666) q[0];
sx q[0];
rz(0.85710454) q[0];
x q[1];
rz(1.3607499) q[2];
sx q[2];
rz(-1.577652) q[2];
sx q[2];
rz(3.1355646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.036025612) q[1];
sx q[1];
rz(-2.470825) q[1];
sx q[1];
rz(0.48700602) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2180473) q[3];
sx q[3];
rz(-2.1622373) q[3];
sx q[3];
rz(-0.54852911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.99737793) q[2];
sx q[2];
rz(-0.37313676) q[2];
sx q[2];
rz(2.0443661) q[2];
rz(1.0002452) q[3];
sx q[3];
rz(-0.92366832) q[3];
sx q[3];
rz(1.3057115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10709396) q[0];
sx q[0];
rz(-1.9871563) q[0];
sx q[0];
rz(-0.19919285) q[0];
rz(-1.2983324) q[1];
sx q[1];
rz(-0.36111626) q[1];
sx q[1];
rz(2.4995506) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10810971) q[0];
sx q[0];
rz(-2.3831522) q[0];
sx q[0];
rz(-1.2940501) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0986311) q[2];
sx q[2];
rz(-1.9342074) q[2];
sx q[2];
rz(0.052415457) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1470829) q[1];
sx q[1];
rz(-1.6509984) q[1];
sx q[1];
rz(2.9417646) q[1];
rz(-1.1690833) q[3];
sx q[3];
rz(-0.48334941) q[3];
sx q[3];
rz(-0.95248896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.22884998) q[2];
sx q[2];
rz(-1.0656837) q[2];
sx q[2];
rz(1.9672811) q[2];
rz(-2.2187388) q[3];
sx q[3];
rz(-2.703981) q[3];
sx q[3];
rz(-2.9722884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4787503) q[0];
sx q[0];
rz(-1.4816544) q[0];
sx q[0];
rz(-2.1424275) q[0];
rz(2.303458) q[1];
sx q[1];
rz(-2.2331388) q[1];
sx q[1];
rz(2.1160486) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9860172) q[0];
sx q[0];
rz(-2.3087569) q[0];
sx q[0];
rz(0.9132333) q[0];
x q[1];
rz(1.8286669) q[2];
sx q[2];
rz(-1.317655) q[2];
sx q[2];
rz(-2.2035901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57069381) q[1];
sx q[1];
rz(-2.0624196) q[1];
sx q[1];
rz(2.6925283) q[1];
rz(0.59896627) q[3];
sx q[3];
rz(-2.3271797) q[3];
sx q[3];
rz(-1.6391476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.020236882) q[2];
sx q[2];
rz(-1.3573703) q[2];
sx q[2];
rz(2.9244002) q[2];
rz(-1.1413261) q[3];
sx q[3];
rz(-2.7666028) q[3];
sx q[3];
rz(-2.2264437) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4311669) q[0];
sx q[0];
rz(-1.2942261) q[0];
sx q[0];
rz(-3.1173832) q[0];
rz(1.0912033) q[1];
sx q[1];
rz(-2.0228491) q[1];
sx q[1];
rz(-2.3275163) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41939999) q[0];
sx q[0];
rz(-2.2792313) q[0];
sx q[0];
rz(-1.6838151) q[0];
rz(0.54583728) q[2];
sx q[2];
rz(-1.0190735) q[2];
sx q[2];
rz(1.1941225) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7660038) q[1];
sx q[1];
rz(-0.75889041) q[1];
sx q[1];
rz(-1.5726456) q[1];
rz(-pi) q[2];
rz(1.7992448) q[3];
sx q[3];
rz(-1.486393) q[3];
sx q[3];
rz(-0.24385246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39390627) q[2];
sx q[2];
rz(-0.86204356) q[2];
sx q[2];
rz(-1.5300592) q[2];
rz(-2.090442) q[3];
sx q[3];
rz(-1.3408778) q[3];
sx q[3];
rz(-0.10147258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6018588) q[0];
sx q[0];
rz(-2.1506385) q[0];
sx q[0];
rz(-2.6644326) q[0];
rz(1.5513783) q[1];
sx q[1];
rz(-1.4789707) q[1];
sx q[1];
rz(-1.8345376) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4068749) q[0];
sx q[0];
rz(-0.85556036) q[0];
sx q[0];
rz(0.060013219) q[0];
x q[1];
rz(-0.81002323) q[2];
sx q[2];
rz(-0.56213435) q[2];
sx q[2];
rz(0.13371828) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9809147) q[1];
sx q[1];
rz(-1.0458993) q[1];
sx q[1];
rz(2.4452433) q[1];
x q[2];
rz(-0.99193345) q[3];
sx q[3];
rz(-2.4483878) q[3];
sx q[3];
rz(-2.7424911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9755379) q[2];
sx q[2];
rz(-2.1259978) q[2];
sx q[2];
rz(-2.3892152) q[2];
rz(-3.0289529) q[3];
sx q[3];
rz(-2.967716) q[3];
sx q[3];
rz(1.1627452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1174711) q[0];
sx q[0];
rz(-1.385067) q[0];
sx q[0];
rz(-2.0140482) q[0];
rz(-0.37359259) q[1];
sx q[1];
rz(-1.5126546) q[1];
sx q[1];
rz(1.0027813) q[1];
rz(-1.1524947) q[2];
sx q[2];
rz(-1.9807182) q[2];
sx q[2];
rz(2.9302927) q[2];
rz(-1.7711025) q[3];
sx q[3];
rz(-2.3080993) q[3];
sx q[3];
rz(-0.75626683) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
