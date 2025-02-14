OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1516079) q[0];
sx q[0];
rz(-0.94010544) q[0];
sx q[0];
rz(0.54036933) q[0];
rz(-2.6481533) q[1];
sx q[1];
rz(-0.72055888) q[1];
sx q[1];
rz(0.61385733) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2342398) q[0];
sx q[0];
rz(-2.1932903) q[0];
sx q[0];
rz(-1.7982152) q[0];
rz(-1.5510606) q[2];
sx q[2];
rz(-2.8497549) q[2];
sx q[2];
rz(0.29668754) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7481374) q[1];
sx q[1];
rz(-0.84987133) q[1];
sx q[1];
rz(-1.1925405) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5491539) q[3];
sx q[3];
rz(-1.8179245) q[3];
sx q[3];
rz(-0.46268845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2962239) q[2];
sx q[2];
rz(-1.8989398) q[2];
sx q[2];
rz(-2.5073012) q[2];
rz(0.93572179) q[3];
sx q[3];
rz(-0.24565419) q[3];
sx q[3];
rz(-0.74265695) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82035404) q[0];
sx q[0];
rz(-1.5411493) q[0];
sx q[0];
rz(-0.4441922) q[0];
rz(-2.3815637) q[1];
sx q[1];
rz(-2.0030231) q[1];
sx q[1];
rz(-2.1601423) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8151756) q[0];
sx q[0];
rz(-1.8627286) q[0];
sx q[0];
rz(2.6089704) q[0];
rz(-pi) q[1];
rz(-1.74238) q[2];
sx q[2];
rz(-1.3918456) q[2];
sx q[2];
rz(1.6323665) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.54320286) q[1];
sx q[1];
rz(-1.1335422) q[1];
sx q[1];
rz(-1.6679126) q[1];
rz(-pi) q[2];
rz(2.7090453) q[3];
sx q[3];
rz(-1.3290231) q[3];
sx q[3];
rz(1.6530619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11357073) q[2];
sx q[2];
rz(-2.5271723) q[2];
sx q[2];
rz(-1.2051955) q[2];
rz(-0.53660721) q[3];
sx q[3];
rz(-1.9034932) q[3];
sx q[3];
rz(-1.4046148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2236915) q[0];
sx q[0];
rz(-1.8603928) q[0];
sx q[0];
rz(2.3028288) q[0];
rz(-0.63610786) q[1];
sx q[1];
rz(-1.6517703) q[1];
sx q[1];
rz(0.020523358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18552173) q[0];
sx q[0];
rz(-0.98617109) q[0];
sx q[0];
rz(-0.67618518) q[0];
rz(-pi) q[1];
rz(0.6095992) q[2];
sx q[2];
rz(-2.4188899) q[2];
sx q[2];
rz(-2.056207) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.62333362) q[1];
sx q[1];
rz(-1.9837399) q[1];
sx q[1];
rz(-1.5790126) q[1];
rz(-pi) q[2];
rz(-2.0666615) q[3];
sx q[3];
rz(-0.32204667) q[3];
sx q[3];
rz(-2.5888458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3402349) q[2];
sx q[2];
rz(-1.6131718) q[2];
sx q[2];
rz(-2.197263) q[2];
rz(2.1186192) q[3];
sx q[3];
rz(-1.2228271) q[3];
sx q[3];
rz(2.003722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.870938) q[0];
sx q[0];
rz(-1.9671257) q[0];
sx q[0];
rz(1.9842072) q[0];
rz(-1.4986787) q[1];
sx q[1];
rz(-1.4721556) q[1];
sx q[1];
rz(-1.9532983) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76283264) q[0];
sx q[0];
rz(-2.4597557) q[0];
sx q[0];
rz(-0.76343151) q[0];
x q[1];
rz(0.54862587) q[2];
sx q[2];
rz(-1.6092192) q[2];
sx q[2];
rz(0.98603934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21056023) q[1];
sx q[1];
rz(-1.0068486) q[1];
sx q[1];
rz(2.2497247) q[1];
rz(1.2318939) q[3];
sx q[3];
rz(-1.7247685) q[3];
sx q[3];
rz(3.0018501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0231861) q[2];
sx q[2];
rz(-1.1772757) q[2];
sx q[2];
rz(-2.4510621) q[2];
rz(-2.0265419) q[3];
sx q[3];
rz(-0.77459049) q[3];
sx q[3];
rz(1.640865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0951776) q[0];
sx q[0];
rz(-1.4354118) q[0];
sx q[0];
rz(2.114356) q[0];
rz(-1.8401624) q[1];
sx q[1];
rz(-2.4899028) q[1];
sx q[1];
rz(-0.32381907) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3394748) q[0];
sx q[0];
rz(-2.3314752) q[0];
sx q[0];
rz(-1.0056061) q[0];
x q[1];
rz(1.5478163) q[2];
sx q[2];
rz(-2.5886664) q[2];
sx q[2];
rz(-1.9000017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65460881) q[1];
sx q[1];
rz(-1.2506079) q[1];
sx q[1];
rz(2.5694808) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15669723) q[3];
sx q[3];
rz(-2.0771871) q[3];
sx q[3];
rz(-1.8739669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7258437) q[2];
sx q[2];
rz(-0.21012935) q[2];
sx q[2];
rz(0.48247639) q[2];
rz(1.9735362) q[3];
sx q[3];
rz(-1.9246293) q[3];
sx q[3];
rz(0.67679685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93550682) q[0];
sx q[0];
rz(-0.85010234) q[0];
sx q[0];
rz(0.8859984) q[0];
rz(2.50792) q[1];
sx q[1];
rz(-0.88992563) q[1];
sx q[1];
rz(2.450313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6768178) q[0];
sx q[0];
rz(-0.64660836) q[0];
sx q[0];
rz(3.0809513) q[0];
x q[1];
rz(-0.076926546) q[2];
sx q[2];
rz(-2.4483213) q[2];
sx q[2];
rz(-1.6037387) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4659652) q[1];
sx q[1];
rz(-1.6576997) q[1];
sx q[1];
rz(-0.81710941) q[1];
rz(2.8915358) q[3];
sx q[3];
rz(-0.80509201) q[3];
sx q[3];
rz(2.6893733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1401691) q[2];
sx q[2];
rz(-2.3052577) q[2];
sx q[2];
rz(-2.8720065) q[2];
rz(-0.94830281) q[3];
sx q[3];
rz(-1.6092665) q[3];
sx q[3];
rz(1.74291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36169323) q[0];
sx q[0];
rz(-0.43710709) q[0];
sx q[0];
rz(-2.1345188) q[0];
rz(-0.40564793) q[1];
sx q[1];
rz(-2.5463153) q[1];
sx q[1];
rz(-2.5880623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8853698) q[0];
sx q[0];
rz(-1.5134583) q[0];
sx q[0];
rz(1.3159412) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78041665) q[2];
sx q[2];
rz(-2.4274106) q[2];
sx q[2];
rz(1.2250021) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2731253) q[1];
sx q[1];
rz(-0.9573862) q[1];
sx q[1];
rz(2.1347743) q[1];
x q[2];
rz(2.7615879) q[3];
sx q[3];
rz(-1.5940574) q[3];
sx q[3];
rz(1.301736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8290528) q[2];
sx q[2];
rz(-2.0487787) q[2];
sx q[2];
rz(1.2139758) q[2];
rz(-0.78643262) q[3];
sx q[3];
rz(-1.395547) q[3];
sx q[3];
rz(-3.1184375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19602747) q[0];
sx q[0];
rz(-2.8493311) q[0];
sx q[0];
rz(1.3319525) q[0];
rz(2.5281483) q[1];
sx q[1];
rz(-2.1211801) q[1];
sx q[1];
rz(2.6920998) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17380781) q[0];
sx q[0];
rz(-1.0840084) q[0];
sx q[0];
rz(0.094102458) q[0];
x q[1];
rz(-1.9067326) q[2];
sx q[2];
rz(-1.6523696) q[2];
sx q[2];
rz(2.5334266) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.044607) q[1];
sx q[1];
rz(-1.7076483) q[1];
sx q[1];
rz(-0.78949331) q[1];
rz(-pi) q[2];
rz(-1.6291796) q[3];
sx q[3];
rz(-1.3135785) q[3];
sx q[3];
rz(1.2793878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.31676644) q[2];
sx q[2];
rz(-0.68244857) q[2];
sx q[2];
rz(0.60834926) q[2];
rz(-1.9781205) q[3];
sx q[3];
rz(-1.7721662) q[3];
sx q[3];
rz(-0.56330645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.0332396) q[0];
sx q[0];
rz(-0.9110564) q[0];
sx q[0];
rz(2.5392927) q[0];
rz(0.96744084) q[1];
sx q[1];
rz(-0.93092218) q[1];
sx q[1];
rz(-1.1036576) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7129015) q[0];
sx q[0];
rz(-1.0279623) q[0];
sx q[0];
rz(0.45962861) q[0];
x q[1];
rz(3.121762) q[2];
sx q[2];
rz(-2.5049372) q[2];
sx q[2];
rz(3.0002468) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5491268) q[1];
sx q[1];
rz(-1.571972) q[1];
sx q[1];
rz(-1.1196613) q[1];
rz(2.1488701) q[3];
sx q[3];
rz(-1.5141271) q[3];
sx q[3];
rz(0.68777675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2063107) q[2];
sx q[2];
rz(-1.3056359) q[2];
sx q[2];
rz(0.3375816) q[2];
rz(2.857699) q[3];
sx q[3];
rz(-2.4604535) q[3];
sx q[3];
rz(0.18686992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5481446) q[0];
sx q[0];
rz(-1.5080867) q[0];
sx q[0];
rz(3.0737851) q[0];
rz(1.1495122) q[1];
sx q[1];
rz(-1.5905453) q[1];
sx q[1];
rz(2.6670719) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15856811) q[0];
sx q[0];
rz(-1.8095542) q[0];
sx q[0];
rz(-0.013834133) q[0];
rz(-pi) q[1];
rz(0.24093012) q[2];
sx q[2];
rz(-1.7027156) q[2];
sx q[2];
rz(-1.2031738) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.671512) q[1];
sx q[1];
rz(-1.5776411) q[1];
sx q[1];
rz(1.6575302) q[1];
x q[2];
rz(2.4429604) q[3];
sx q[3];
rz(-2.1967874) q[3];
sx q[3];
rz(-2.0076942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48156753) q[2];
sx q[2];
rz(-0.93067545) q[2];
sx q[2];
rz(1.0732667) q[2];
rz(-0.16452161) q[3];
sx q[3];
rz(-1.3273032) q[3];
sx q[3];
rz(-2.8209414) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58707033) q[0];
sx q[0];
rz(-1.5656492) q[0];
sx q[0];
rz(1.5026305) q[0];
rz(1.5578237) q[1];
sx q[1];
rz(-2.0638034) q[1];
sx q[1];
rz(2.6149909) q[1];
rz(1.2100958) q[2];
sx q[2];
rz(-2.0059235) q[2];
sx q[2];
rz(-0.37034482) q[2];
rz(2.6216636) q[3];
sx q[3];
rz(-0.16362301) q[3];
sx q[3];
rz(2.5273821) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
