OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.36800185) q[0];
sx q[0];
rz(-0.79080963) q[0];
sx q[0];
rz(0.33413449) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96555644) q[0];
sx q[0];
rz(-2.3210225) q[0];
sx q[0];
rz(-1.6198938) q[0];
x q[1];
rz(-1.5300418) q[2];
sx q[2];
rz(-2.0353122) q[2];
sx q[2];
rz(-1.0711311) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76518607) q[1];
sx q[1];
rz(-1.1131439) q[1];
sx q[1];
rz(-0.54995579) q[1];
rz(-2.9524607) q[3];
sx q[3];
rz(-1.3526219) q[3];
sx q[3];
rz(-1.8060341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.618764) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(-2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537271) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(2.7541449) q[0];
rz(-0.93915835) q[1];
sx q[1];
rz(-2.1444131) q[1];
sx q[1];
rz(-1.739025) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99968602) q[0];
sx q[0];
rz(-1.5667331) q[0];
sx q[0];
rz(-0.02948972) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3520794) q[2];
sx q[2];
rz(-2.1195076) q[2];
sx q[2];
rz(0.95900853) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.034033876) q[1];
sx q[1];
rz(-0.65980136) q[1];
sx q[1];
rz(-2.4504513) q[1];
rz(-pi) q[2];
rz(3.0121147) q[3];
sx q[3];
rz(-2.5054512) q[3];
sx q[3];
rz(0.57487088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42276057) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(2.823901) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(-0.40107045) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6221878) q[0];
sx q[0];
rz(-1.6739968) q[0];
sx q[0];
rz(-1.897057) q[0];
rz(-pi) q[1];
rz(-1.3813854) q[2];
sx q[2];
rz(-0.93511287) q[2];
sx q[2];
rz(2.1798101) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33441015) q[1];
sx q[1];
rz(-1.5644329) q[1];
sx q[1];
rz(2.4114354) q[1];
rz(-pi) q[2];
rz(1.1776351) q[3];
sx q[3];
rz(-1.7522246) q[3];
sx q[3];
rz(-1.7593311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59427375) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(2.1650971) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(0.78021375) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926086) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(1.6216888) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(-2.8881883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4356416) q[0];
sx q[0];
rz(-0.55103978) q[0];
sx q[0];
rz(-1.7680697) q[0];
x q[1];
rz(1.0850111) q[2];
sx q[2];
rz(-1.5554785) q[2];
sx q[2];
rz(-2.3682396) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.036382347) q[1];
sx q[1];
rz(-1.6670334) q[1];
sx q[1];
rz(0.63016816) q[1];
rz(-2.4921791) q[3];
sx q[3];
rz(-1.6522437) q[3];
sx q[3];
rz(0.34929517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(-2.3542662) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(-1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(-0.45809349) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23819085) q[0];
sx q[0];
rz(-1.5731249) q[0];
sx q[0];
rz(-0.44263126) q[0];
rz(-pi) q[1];
rz(0.96180054) q[2];
sx q[2];
rz(-0.70448175) q[2];
sx q[2];
rz(1.5722164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0702857) q[1];
sx q[1];
rz(-0.83155337) q[1];
sx q[1];
rz(-2.9830095) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3985653) q[3];
sx q[3];
rz(-0.50695626) q[3];
sx q[3];
rz(0.34190049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1725585) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(-2.9120973) q[2];
rz(3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5620419) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(-0.91947412) q[0];
rz(-0.062285034) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(1.2671635) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125268) q[0];
sx q[0];
rz(-2.292284) q[0];
sx q[0];
rz(0.52961403) q[0];
x q[1];
rz(2.9847449) q[2];
sx q[2];
rz(-2.3914797) q[2];
sx q[2];
rz(1.4175121) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3082723) q[1];
sx q[1];
rz(-2.7566524) q[1];
sx q[1];
rz(-2.1182548) q[1];
rz(-2.0165765) q[3];
sx q[3];
rz(-0.52281724) q[3];
sx q[3];
rz(-2.0857874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1798114) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-2.5578257) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0631183) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(2.069058) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(2.0297208) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3355013) q[0];
sx q[0];
rz(-1.4772692) q[0];
sx q[0];
rz(2.0557687) q[0];
rz(-pi) q[1];
rz(1.9260336) q[2];
sx q[2];
rz(-1.0676427) q[2];
sx q[2];
rz(-0.1728729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.48549451) q[1];
sx q[1];
rz(-1.6472677) q[1];
sx q[1];
rz(2.8511726) q[1];
x q[2];
rz(-2.2124347) q[3];
sx q[3];
rz(-0.60141364) q[3];
sx q[3];
rz(1.8638368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(1.1676577) q[2];
rz(1.6052823) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.4060098) q[3];
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
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(0.73079601) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(2.3866167) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98159957) q[0];
sx q[0];
rz(-2.2845075) q[0];
sx q[0];
rz(0.1937565) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72889363) q[2];
sx q[2];
rz(-1.1155827) q[2];
sx q[2];
rz(-0.58065562) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4432115) q[1];
sx q[1];
rz(-1.4581919) q[1];
sx q[1];
rz(2.5618782) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8927535) q[3];
sx q[3];
rz(-2.4340981) q[3];
sx q[3];
rz(-1.0196109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76453152) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(0.63684741) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4144142) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(2.0027347) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(0.019502217) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6396128) q[0];
sx q[0];
rz(-1.4907955) q[0];
sx q[0];
rz(1.3638858) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9901864) q[2];
sx q[2];
rz(-1.7689686) q[2];
sx q[2];
rz(2.1422447) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1972678) q[1];
sx q[1];
rz(-1.7515344) q[1];
sx q[1];
rz(0.9428057) q[1];
rz(-pi) q[2];
rz(-1.8832302) q[3];
sx q[3];
rz(-1.4133246) q[3];
sx q[3];
rz(-2.0563994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1372244) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(-0.84890378) q[2];
rz(2.7539339) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(-1.7548521) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.9932995) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0120221) q[0];
sx q[0];
rz(-0.14229933) q[0];
sx q[0];
rz(2.9905031) q[0];
rz(-1.7933153) q[2];
sx q[2];
rz(-2.1902124) q[2];
sx q[2];
rz(-0.36048181) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1710098) q[1];
sx q[1];
rz(-2.2323771) q[1];
sx q[1];
rz(-2.9243102) q[1];
x q[2];
rz(0.34928068) q[3];
sx q[3];
rz(-1.5216856) q[3];
sx q[3];
rz(-1.9625361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-0.36515507) q[2];
rz(3.0129516) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(-0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(3.0317595) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-2.1784492) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(-1.1744432) q[2];
sx q[2];
rz(-0.94559961) q[2];
sx q[2];
rz(-0.0018975817) q[2];
rz(0.036988463) q[3];
sx q[3];
rz(-1.0109517) q[3];
sx q[3];
rz(3.0299822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
