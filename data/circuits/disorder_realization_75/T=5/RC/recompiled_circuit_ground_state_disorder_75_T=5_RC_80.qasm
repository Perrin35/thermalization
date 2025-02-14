OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45159856) q[0];
sx q[0];
rz(-0.30878433) q[0];
sx q[0];
rz(-0.2398332) q[0];
rz(-3.7026703) q[1];
sx q[1];
rz(2.3993888) q[1];
sx q[1];
rz(11.053434) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6282677) q[0];
sx q[0];
rz(-1.4428992) q[0];
sx q[0];
rz(1.2012175) q[0];
rz(-2.0432908) q[2];
sx q[2];
rz(-1.3368946) q[2];
sx q[2];
rz(-0.59532524) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1354243) q[1];
sx q[1];
rz(-1.0841771) q[1];
sx q[1];
rz(-1.4044589) q[1];
rz(-pi) q[2];
rz(1.2414819) q[3];
sx q[3];
rz(-0.17767492) q[3];
sx q[3];
rz(0.8575646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98770398) q[2];
sx q[2];
rz(-0.99606267) q[2];
sx q[2];
rz(1.055701) q[2];
rz(-2.740247) q[3];
sx q[3];
rz(-1.6258806) q[3];
sx q[3];
rz(-1.5041941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6317247) q[0];
sx q[0];
rz(-2.8724176) q[0];
sx q[0];
rz(-1.4058231) q[0];
rz(0.74613219) q[1];
sx q[1];
rz(-1.965062) q[1];
sx q[1];
rz(3.0768118) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5279605) q[0];
sx q[0];
rz(-2.0191231) q[0];
sx q[0];
rz(-1.1290068) q[0];
x q[1];
rz(-2.9322036) q[2];
sx q[2];
rz(-0.69062606) q[2];
sx q[2];
rz(-0.90589452) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3105615) q[1];
sx q[1];
rz(-1.4906724) q[1];
sx q[1];
rz(1.4270272) q[1];
rz(-1.4218036) q[3];
sx q[3];
rz(-0.21315609) q[3];
sx q[3];
rz(-3.0557291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.259321) q[2];
sx q[2];
rz(-2.6210625) q[2];
sx q[2];
rz(-3.1003013) q[2];
rz(-1.1193554) q[3];
sx q[3];
rz(-1.7740403) q[3];
sx q[3];
rz(2.0004499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3248046) q[0];
sx q[0];
rz(-0.4275221) q[0];
sx q[0];
rz(2.4422755) q[0];
rz(2.518867) q[1];
sx q[1];
rz(-2.3284349) q[1];
sx q[1];
rz(2.8831388) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39614933) q[0];
sx q[0];
rz(-0.93702261) q[0];
sx q[0];
rz(-1.5556015) q[0];
x q[1];
rz(1.7202366) q[2];
sx q[2];
rz(-0.65776134) q[2];
sx q[2];
rz(0.8873111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.099951) q[1];
sx q[1];
rz(-2.0196242) q[1];
sx q[1];
rz(-1.5931409) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0007802) q[3];
sx q[3];
rz(-1.6039492) q[3];
sx q[3];
rz(1.6044817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2059325) q[2];
sx q[2];
rz(-1.460133) q[2];
sx q[2];
rz(0.70651954) q[2];
rz(-2.5126854) q[3];
sx q[3];
rz(-2.3157412) q[3];
sx q[3];
rz(-2.0188873) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0141456) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(2.1674147) q[0];
rz(1.4840508) q[1];
sx q[1];
rz(-2.0482792) q[1];
sx q[1];
rz(-1.0008224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054996252) q[0];
sx q[0];
rz(-0.97584134) q[0];
sx q[0];
rz(0.18878285) q[0];
x q[1];
rz(-0.77032178) q[2];
sx q[2];
rz(-0.97409596) q[2];
sx q[2];
rz(1.53231) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82881935) q[1];
sx q[1];
rz(-2.3760894) q[1];
sx q[1];
rz(0.70154066) q[1];
rz(-2.9934512) q[3];
sx q[3];
rz(-0.32399789) q[3];
sx q[3];
rz(-1.2202386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4546844) q[2];
sx q[2];
rz(-1.800622) q[2];
sx q[2];
rz(2.3243813) q[2];
rz(0.59598437) q[3];
sx q[3];
rz(-1.7242566) q[3];
sx q[3];
rz(-1.7433085) q[3];
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
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3572094) q[0];
sx q[0];
rz(-1.7359808) q[0];
sx q[0];
rz(1.0719517) q[0];
rz(-2.0042073) q[1];
sx q[1];
rz(-1.6268566) q[1];
sx q[1];
rz(-0.17328182) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6023077) q[0];
sx q[0];
rz(-2.0208911) q[0];
sx q[0];
rz(-2.0099239) q[0];
rz(-pi) q[1];
rz(2.7543147) q[2];
sx q[2];
rz(-1.1114632) q[2];
sx q[2];
rz(-1.1371374) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3751602) q[1];
sx q[1];
rz(-0.81534895) q[1];
sx q[1];
rz(-0.098618193) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4777484) q[3];
sx q[3];
rz(-1.4085839) q[3];
sx q[3];
rz(-1.4402267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5002354) q[2];
sx q[2];
rz(-0.45318979) q[2];
sx q[2];
rz(0.87453169) q[2];
rz(-2.4747961) q[3];
sx q[3];
rz(-1.4614481) q[3];
sx q[3];
rz(-2.3165406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2905529) q[0];
sx q[0];
rz(-2.9930826) q[0];
sx q[0];
rz(0.17459757) q[0];
rz(0.54706508) q[1];
sx q[1];
rz(-2.6262296) q[1];
sx q[1];
rz(-1.863716) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9491315) q[0];
sx q[0];
rz(-1.2843848) q[0];
sx q[0];
rz(-3.0770296) q[0];
rz(0.49680423) q[2];
sx q[2];
rz(-2.3251109) q[2];
sx q[2];
rz(1.7330488) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4365579) q[1];
sx q[1];
rz(-2.3113657) q[1];
sx q[1];
rz(-2.6565246) q[1];
rz(-pi) q[2];
rz(-0.61056633) q[3];
sx q[3];
rz(-0.99733099) q[3];
sx q[3];
rz(1.8449044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7586907) q[2];
sx q[2];
rz(-0.26472696) q[2];
sx q[2];
rz(-0.36941377) q[2];
rz(2.3139125) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(0.15414342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.5675885) q[0];
sx q[0];
rz(-2.2249157) q[0];
sx q[0];
rz(-2.9045203) q[0];
rz(-1.392662) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(-1.0038092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8477551) q[0];
sx q[0];
rz(-1.4042353) q[0];
sx q[0];
rz(0.44035797) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5215724) q[2];
sx q[2];
rz(-0.78465377) q[2];
sx q[2];
rz(-1.705738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8449515) q[1];
sx q[1];
rz(-1.2403204) q[1];
sx q[1];
rz(-0.56092324) q[1];
rz(0.5587033) q[3];
sx q[3];
rz(-1.1538343) q[3];
sx q[3];
rz(0.95181634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6134593) q[2];
sx q[2];
rz(-1.5626835) q[2];
sx q[2];
rz(0.51631874) q[2];
rz(2.3033219) q[3];
sx q[3];
rz(-1.6603989) q[3];
sx q[3];
rz(-2.9576438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2876005) q[0];
sx q[0];
rz(-2.8599399) q[0];
sx q[0];
rz(2.2273492) q[0];
rz(-2.5573348) q[1];
sx q[1];
rz(-1.7353568) q[1];
sx q[1];
rz(-0.90726888) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6833008) q[0];
sx q[0];
rz(-2.2193546) q[0];
sx q[0];
rz(3.1187727) q[0];
x q[1];
rz(-1.3970988) q[2];
sx q[2];
rz(-0.96250421) q[2];
sx q[2];
rz(-2.4243958) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3104035) q[1];
sx q[1];
rz(-2.0210279) q[1];
sx q[1];
rz(-2.4163209) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26772883) q[3];
sx q[3];
rz(-0.98891034) q[3];
sx q[3];
rz(-1.8393593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82627901) q[2];
sx q[2];
rz(-1.7118914) q[2];
sx q[2];
rz(0.86177525) q[2];
rz(-1.2365384) q[3];
sx q[3];
rz(-3.0573513) q[3];
sx q[3];
rz(-1.2110565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009907) q[0];
sx q[0];
rz(-2.3589098) q[0];
sx q[0];
rz(-2.1114517) q[0];
rz(2.8252699) q[1];
sx q[1];
rz(-1.373469) q[1];
sx q[1];
rz(2.8533459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.41053) q[0];
sx q[0];
rz(-1.9707435) q[0];
sx q[0];
rz(1.4131141) q[0];
rz(-pi) q[1];
rz(0.64916237) q[2];
sx q[2];
rz(-1.1466951) q[2];
sx q[2];
rz(2.8555388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94142524) q[1];
sx q[1];
rz(-1.9899448) q[1];
sx q[1];
rz(0.9113542) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8602198) q[3];
sx q[3];
rz(-2.4503539) q[3];
sx q[3];
rz(2.6760677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2890702) q[2];
sx q[2];
rz(-1.4129637) q[2];
sx q[2];
rz(2.9150325) q[2];
rz(-1.6864927) q[3];
sx q[3];
rz(-0.89613599) q[3];
sx q[3];
rz(2.4494825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15014547) q[0];
sx q[0];
rz(-1.2757855) q[0];
sx q[0];
rz(2.5164497) q[0];
rz(-1.9236247) q[1];
sx q[1];
rz(-0.70911276) q[1];
sx q[1];
rz(-1.2120754) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9218413) q[0];
sx q[0];
rz(-0.80972396) q[0];
sx q[0];
rz(-0.68897665) q[0];
x q[1];
rz(1.3122968) q[2];
sx q[2];
rz(-2.2257559) q[2];
sx q[2];
rz(-1.0884681) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.068262488) q[1];
sx q[1];
rz(-2.0113328) q[1];
sx q[1];
rz(3.0777626) q[1];
rz(-pi) q[2];
rz(-0.76833581) q[3];
sx q[3];
rz(-2.7339122) q[3];
sx q[3];
rz(2.3967495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7071699) q[2];
sx q[2];
rz(-1.7662798) q[2];
sx q[2];
rz(2.0909069) q[2];
rz(-2.5896942) q[3];
sx q[3];
rz(-2.4514908) q[3];
sx q[3];
rz(-1.2845854) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62591775) q[0];
sx q[0];
rz(-0.82763012) q[0];
sx q[0];
rz(-0.81830842) q[0];
rz(-0.7863518) q[1];
sx q[1];
rz(-0.66023371) q[1];
sx q[1];
rz(0.22088851) q[1];
rz(-3.1037504) q[2];
sx q[2];
rz(-1.9339682) q[2];
sx q[2];
rz(0.54724271) q[2];
rz(0.10765392) q[3];
sx q[3];
rz(-0.68895491) q[3];
sx q[3];
rz(0.72910492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
