OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(-1.5664772) q[0];
sx q[0];
rz(-1.1251261) q[0];
rz(1.0220802) q[1];
sx q[1];
rz(-0.66749579) q[1];
sx q[1];
rz(-2.3281085) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4018702) q[0];
sx q[0];
rz(-0.39416322) q[0];
sx q[0];
rz(-1.9456844) q[0];
x q[1];
rz(-1.9314693) q[2];
sx q[2];
rz(-1.1481592) q[2];
sx q[2];
rz(-2.8354743) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58373815) q[1];
sx q[1];
rz(-2.0101133) q[1];
sx q[1];
rz(-1.735461) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0510819) q[3];
sx q[3];
rz(-1.7337243) q[3];
sx q[3];
rz(0.75608692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6525314) q[2];
sx q[2];
rz(-1.4514613) q[2];
sx q[2];
rz(-2.2129464) q[2];
rz(1.5993902) q[3];
sx q[3];
rz(-1.3336811) q[3];
sx q[3];
rz(-1.9699875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9376675) q[0];
sx q[0];
rz(-1.3805905) q[0];
sx q[0];
rz(-3.0145338) q[0];
rz(2.1584885) q[1];
sx q[1];
rz(-1.7763205) q[1];
sx q[1];
rz(-2.3703221) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5708904) q[0];
sx q[0];
rz(-1.9237907) q[0];
sx q[0];
rz(-2.3285026) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7423058) q[2];
sx q[2];
rz(-2.7993188) q[2];
sx q[2];
rz(-2.8245087) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82133085) q[1];
sx q[1];
rz(-1.2148569) q[1];
sx q[1];
rz(-1.9411646) q[1];
rz(-0.57621376) q[3];
sx q[3];
rz(-2.0348843) q[3];
sx q[3];
rz(2.4903542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9942921) q[2];
sx q[2];
rz(-1.9376126) q[2];
sx q[2];
rz(0.0083943923) q[2];
rz(2.4781135) q[3];
sx q[3];
rz(-1.2365664) q[3];
sx q[3];
rz(-2.8765163) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10107772) q[0];
sx q[0];
rz(-0.82656693) q[0];
sx q[0];
rz(-2.7000632) q[0];
rz(2.1532374) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(-3.006014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044047318) q[0];
sx q[0];
rz(-1.5870461) q[0];
sx q[0];
rz(1.5784997) q[0];
rz(-pi) q[1];
rz(2.709948) q[2];
sx q[2];
rz(-2.3767391) q[2];
sx q[2];
rz(2.5333196) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24568711) q[1];
sx q[1];
rz(-2.062633) q[1];
sx q[1];
rz(0.40426429) q[1];
rz(-pi) q[2];
rz(-1.7884004) q[3];
sx q[3];
rz(-0.74624589) q[3];
sx q[3];
rz(1.7640132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93379891) q[2];
sx q[2];
rz(-1.7094882) q[2];
sx q[2];
rz(3.1033893) q[2];
rz(0.52538747) q[3];
sx q[3];
rz(-0.62756687) q[3];
sx q[3];
rz(0.6853404) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66048375) q[0];
sx q[0];
rz(-0.79367343) q[0];
sx q[0];
rz(-2.5307181) q[0];
rz(1.3350217) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(0.23922051) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6253875) q[0];
sx q[0];
rz(-0.62407485) q[0];
sx q[0];
rz(1.1136057) q[0];
rz(-2.6692713) q[2];
sx q[2];
rz(-1.2178253) q[2];
sx q[2];
rz(-1.8582839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.661631) q[1];
sx q[1];
rz(-1.1203655) q[1];
sx q[1];
rz(1.9305139) q[1];
x q[2];
rz(-1.7710822) q[3];
sx q[3];
rz(-1.1173751) q[3];
sx q[3];
rz(-0.1399006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7188344) q[2];
sx q[2];
rz(-1.4449291) q[2];
sx q[2];
rz(-3.139843) q[2];
rz(2.8478029) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(2.8747115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2413498) q[0];
sx q[0];
rz(-1.7647864) q[0];
sx q[0];
rz(-0.84306651) q[0];
rz(-0.33755606) q[1];
sx q[1];
rz(-1.866021) q[1];
sx q[1];
rz(-1.5379803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0564954) q[0];
sx q[0];
rz(-0.81696327) q[0];
sx q[0];
rz(0.2635862) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29422167) q[2];
sx q[2];
rz(-0.18202848) q[2];
sx q[2];
rz(0.95532571) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2258721) q[1];
sx q[1];
rz(-0.32494007) q[1];
sx q[1];
rz(-2.48808) q[1];
x q[2];
rz(2.9318453) q[3];
sx q[3];
rz(-1.4919859) q[3];
sx q[3];
rz(-0.94890734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7097077) q[2];
sx q[2];
rz(-2.0543435) q[2];
sx q[2];
rz(-1.0464) q[2];
rz(2.8228068) q[3];
sx q[3];
rz(-2.4304978) q[3];
sx q[3];
rz(-0.54767245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043561291) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(0.64796722) q[0];
rz(-1.3994392) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(-1.3290149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15547046) q[0];
sx q[0];
rz(-0.24674812) q[0];
sx q[0];
rz(1.107649) q[0];
rz(0.29000303) q[2];
sx q[2];
rz(-1.776374) q[2];
sx q[2];
rz(-2.9221591) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4635515) q[1];
sx q[1];
rz(-1.8811418) q[1];
sx q[1];
rz(0.84872679) q[1];
rz(-1.9372609) q[3];
sx q[3];
rz(-2.0559337) q[3];
sx q[3];
rz(1.7462891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7602188) q[2];
sx q[2];
rz(-1.9332644) q[2];
sx q[2];
rz(-1.5488497) q[2];
rz(-3.0717487) q[3];
sx q[3];
rz(-1.8826238) q[3];
sx q[3];
rz(2.2729592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1182564) q[0];
sx q[0];
rz(-1.2160439) q[0];
sx q[0];
rz(0.94183952) q[0];
rz(-0.16009227) q[1];
sx q[1];
rz(-1.6611049) q[1];
sx q[1];
rz(-0.19217415) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036445905) q[0];
sx q[0];
rz(-0.22686401) q[0];
sx q[0];
rz(1.731621) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1015119) q[2];
sx q[2];
rz(-1.049384) q[2];
sx q[2];
rz(-1.5837976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.46017473) q[1];
sx q[1];
rz(-2.3263086) q[1];
sx q[1];
rz(2.6667206) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5302079) q[3];
sx q[3];
rz(-2.6436664) q[3];
sx q[3];
rz(2.9663309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.75752246) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(-2.7080217) q[2];
rz(1.5571669) q[3];
sx q[3];
rz(-1.6470563) q[3];
sx q[3];
rz(2.7092194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(1.7973416) q[0];
rz(1.9380219) q[1];
sx q[1];
rz(-1.1886339) q[1];
sx q[1];
rz(1.7787836) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8289611) q[0];
sx q[0];
rz(-1.5927218) q[0];
sx q[0];
rz(1.6005105) q[0];
rz(-pi) q[1];
rz(2.9672181) q[2];
sx q[2];
rz(-2.3078636) q[2];
sx q[2];
rz(-2.9437609) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4688022) q[1];
sx q[1];
rz(-1.4458477) q[1];
sx q[1];
rz(2.9114257) q[1];
rz(-2.8239488) q[3];
sx q[3];
rz(-2.2322725) q[3];
sx q[3];
rz(0.7484439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5726996) q[2];
sx q[2];
rz(-2.3685444) q[2];
sx q[2];
rz(2.1577238) q[2];
rz(-0.24108663) q[3];
sx q[3];
rz(-0.9328931) q[3];
sx q[3];
rz(-2.4510395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.6702061) q[0];
sx q[0];
rz(-0.79600483) q[0];
sx q[0];
rz(-0.27467003) q[0];
rz(1.341691) q[1];
sx q[1];
rz(-0.62961737) q[1];
sx q[1];
rz(-2.0955657) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2316596) q[0];
sx q[0];
rz(-0.69485661) q[0];
sx q[0];
rz(-2.6207663) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.080950254) q[2];
sx q[2];
rz(-1.3737203) q[2];
sx q[2];
rz(1.0287036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0053802) q[1];
sx q[1];
rz(-2.5544205) q[1];
sx q[1];
rz(-3.0239952) q[1];
rz(1.0677797) q[3];
sx q[3];
rz(-1.5340337) q[3];
sx q[3];
rz(-0.063677841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2785953) q[2];
sx q[2];
rz(-1.7509165) q[2];
sx q[2];
rz(2.7421303) q[2];
rz(-2.5755889) q[3];
sx q[3];
rz(-2.1750906) q[3];
sx q[3];
rz(-0.81542265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852916) q[0];
sx q[0];
rz(-1.4194019) q[0];
sx q[0];
rz(2.5592819) q[0];
rz(-2.9583926) q[1];
sx q[1];
rz(-0.87545005) q[1];
sx q[1];
rz(0.80642548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88993011) q[0];
sx q[0];
rz(-1.4668873) q[0];
sx q[0];
rz(0.088168747) q[0];
x q[1];
rz(-0.66566531) q[2];
sx q[2];
rz(-2.7241926) q[2];
sx q[2];
rz(1.1102499) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6567071) q[1];
sx q[1];
rz(-2.2995512) q[1];
sx q[1];
rz(1.2970112) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34727283) q[3];
sx q[3];
rz(-1.4662305) q[3];
sx q[3];
rz(2.4109651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23239423) q[2];
sx q[2];
rz(-0.71467233) q[2];
sx q[2];
rz(0.8052899) q[2];
rz(0.14690873) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.88937) q[0];
sx q[0];
rz(-1.7563553) q[0];
sx q[0];
rz(1.8808543) q[0];
rz(1.5994785) q[1];
sx q[1];
rz(-1.5972932) q[1];
sx q[1];
rz(1.6398026) q[1];
rz(-1.4412104) q[2];
sx q[2];
rz(-0.72740462) q[2];
sx q[2];
rz(-1.8047756) q[2];
rz(-2.4965548) q[3];
sx q[3];
rz(-1.2923664) q[3];
sx q[3];
rz(-0.016135767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
