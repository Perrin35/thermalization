OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.23586805) q[0];
sx q[0];
rz(1.9960825) q[0];
sx q[0];
rz(9.3762015) q[0];
rz(-1.6254758) q[1];
sx q[1];
rz(-1.0928417) q[1];
sx q[1];
rz(0.71576524) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1224533) q[0];
sx q[0];
rz(-1.4000443) q[0];
sx q[0];
rz(3.1117361) q[0];
rz(2.9346132) q[2];
sx q[2];
rz(-1.2964777) q[2];
sx q[2];
rz(-1.5766608) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29808334) q[1];
sx q[1];
rz(-0.88551003) q[1];
sx q[1];
rz(-2.09649) q[1];
rz(-2.4972102) q[3];
sx q[3];
rz(-0.61840671) q[3];
sx q[3];
rz(-1.1856841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48892659) q[2];
sx q[2];
rz(-0.55997866) q[2];
sx q[2];
rz(-2.0860591) q[2];
rz(2.7355898) q[3];
sx q[3];
rz(-1.3136274) q[3];
sx q[3];
rz(2.3199911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15776289) q[0];
sx q[0];
rz(-1.0394179) q[0];
sx q[0];
rz(2.5854172) q[0];
rz(1.3293386) q[1];
sx q[1];
rz(-0.61873299) q[1];
sx q[1];
rz(3.0583196) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0389685) q[0];
sx q[0];
rz(-2.4100465) q[0];
sx q[0];
rz(2.6464173) q[0];
rz(-pi) q[1];
rz(2.3131392) q[2];
sx q[2];
rz(-0.76702416) q[2];
sx q[2];
rz(-1.3924862) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1412388) q[1];
sx q[1];
rz(-1.3508995) q[1];
sx q[1];
rz(2.6882437) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2733685) q[3];
sx q[3];
rz(-0.69528841) q[3];
sx q[3];
rz(-0.27588683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1978153) q[2];
sx q[2];
rz(-1.752172) q[2];
sx q[2];
rz(-2.6488292) q[2];
rz(-1.0265776) q[3];
sx q[3];
rz(-2.8387098) q[3];
sx q[3];
rz(-1.4125642) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3781085) q[0];
sx q[0];
rz(-0.56270993) q[0];
sx q[0];
rz(-0.43877959) q[0];
rz(-1.6890866) q[1];
sx q[1];
rz(-0.32183281) q[1];
sx q[1];
rz(2.7486393) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753112) q[0];
sx q[0];
rz(-1.5987168) q[0];
sx q[0];
rz(-1.0795388) q[0];
rz(-1.6268756) q[2];
sx q[2];
rz(-2.6125557) q[2];
sx q[2];
rz(-2.7521486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4247228) q[1];
sx q[1];
rz(-0.42415127) q[1];
sx q[1];
rz(1.7440026) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64928253) q[3];
sx q[3];
rz(-1.5074132) q[3];
sx q[3];
rz(-1.5017623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4517333) q[2];
sx q[2];
rz(-2.2971575) q[2];
sx q[2];
rz(0.27352697) q[2];
rz(-2.5290329) q[3];
sx q[3];
rz(-1.7275683) q[3];
sx q[3];
rz(0.85649049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81371152) q[0];
sx q[0];
rz(-0.59309816) q[0];
sx q[0];
rz(2.3408422) q[0];
rz(2.4294991) q[1];
sx q[1];
rz(-2.021603) q[1];
sx q[1];
rz(-0.48316479) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2321076) q[0];
sx q[0];
rz(-1.5714025) q[0];
sx q[0];
rz(2.3988947) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4078385) q[2];
sx q[2];
rz(-0.68063191) q[2];
sx q[2];
rz(-1.5524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9312716) q[1];
sx q[1];
rz(-2.4948947) q[1];
sx q[1];
rz(0.047860459) q[1];
rz(1.2340868) q[3];
sx q[3];
rz(-0.90970618) q[3];
sx q[3];
rz(-2.8291836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.157865) q[2];
sx q[2];
rz(-0.71061504) q[2];
sx q[2];
rz(-1.2845117) q[2];
rz(-0.48786783) q[3];
sx q[3];
rz(-1.7043461) q[3];
sx q[3];
rz(1.454486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81216413) q[0];
sx q[0];
rz(-2.282113) q[0];
sx q[0];
rz(-0.46464768) q[0];
rz(-2.1931785) q[1];
sx q[1];
rz(-0.83810884) q[1];
sx q[1];
rz(-1.4533739) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8888055) q[0];
sx q[0];
rz(-1.545419) q[0];
sx q[0];
rz(-1.5359098) q[0];
rz(-pi) q[1];
rz(-0.74044944) q[2];
sx q[2];
rz(-1.3794823) q[2];
sx q[2];
rz(-2.180661) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0665474) q[1];
sx q[1];
rz(-2.5270224) q[1];
sx q[1];
rz(-0.11654186) q[1];
rz(0.57027633) q[3];
sx q[3];
rz(-0.81109427) q[3];
sx q[3];
rz(2.9549884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1571618) q[2];
sx q[2];
rz(-1.0474019) q[2];
sx q[2];
rz(-2.7471527) q[2];
rz(-1.1557584) q[3];
sx q[3];
rz(-0.93005669) q[3];
sx q[3];
rz(-1.3084779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30094576) q[0];
sx q[0];
rz(-0.41579682) q[0];
sx q[0];
rz(-2.0796602) q[0];
rz(-2.1901219) q[1];
sx q[1];
rz(-2.2679592) q[1];
sx q[1];
rz(-1.5207312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7436696) q[0];
sx q[0];
rz(-0.28215223) q[0];
sx q[0];
rz(1.762853) q[0];
x q[1];
rz(-2.4967147) q[2];
sx q[2];
rz(-2.0487924) q[2];
sx q[2];
rz(-1.0005815) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1870059) q[1];
sx q[1];
rz(-2.329064) q[1];
sx q[1];
rz(1.345326) q[1];
rz(-0.39395515) q[3];
sx q[3];
rz(-1.7749014) q[3];
sx q[3];
rz(2.1597854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1157397) q[2];
sx q[2];
rz(-1.226475) q[2];
sx q[2];
rz(-0.24570492) q[2];
rz(1.4541516) q[3];
sx q[3];
rz(-1.5116296) q[3];
sx q[3];
rz(-0.83034849) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65569735) q[0];
sx q[0];
rz(-2.4035154) q[0];
sx q[0];
rz(-0.6947211) q[0];
rz(-3.0409536) q[1];
sx q[1];
rz(-1.0357608) q[1];
sx q[1];
rz(-0.86520854) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9964906) q[0];
sx q[0];
rz(-1.9887513) q[0];
sx q[0];
rz(-1.5254024) q[0];
x q[1];
rz(0.91036441) q[2];
sx q[2];
rz(-1.5558793) q[2];
sx q[2];
rz(-0.93176022) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76825095) q[1];
sx q[1];
rz(-0.96186559) q[1];
sx q[1];
rz(-0.18280289) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1221625) q[3];
sx q[3];
rz(-2.1688281) q[3];
sx q[3];
rz(-0.74265146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9575017) q[2];
sx q[2];
rz(-2.3881193) q[2];
sx q[2];
rz(0.91020477) q[2];
rz(-1.4522067) q[3];
sx q[3];
rz(-1.2319177) q[3];
sx q[3];
rz(-0.46428251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8978187) q[0];
sx q[0];
rz(-1.047736) q[0];
sx q[0];
rz(2.5183103) q[0];
rz(-0.15577623) q[1];
sx q[1];
rz(-2.364295) q[1];
sx q[1];
rz(2.8782841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1415188) q[0];
sx q[0];
rz(-0.21528582) q[0];
sx q[0];
rz(1.6208642) q[0];
x q[1];
rz(-1.510235) q[2];
sx q[2];
rz(-1.6471631) q[2];
sx q[2];
rz(-1.5983943) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0075078) q[1];
sx q[1];
rz(-2.1459487) q[1];
sx q[1];
rz(-1.6624727) q[1];
x q[2];
rz(1.0271163) q[3];
sx q[3];
rz(-0.86302084) q[3];
sx q[3];
rz(0.88524997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.89173633) q[2];
sx q[2];
rz(-2.8920434) q[2];
sx q[2];
rz(2.2535394) q[2];
rz(2.213721) q[3];
sx q[3];
rz(-2.0136621) q[3];
sx q[3];
rz(2.1685062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71858281) q[0];
sx q[0];
rz(-2.6319478) q[0];
sx q[0];
rz(2.6408559) q[0];
rz(2.0210733) q[1];
sx q[1];
rz(-2.3605533) q[1];
sx q[1];
rz(-0.80148554) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716176) q[0];
sx q[0];
rz(-1.5688251) q[0];
sx q[0];
rz(3.0706186) q[0];
x q[1];
rz(-2.5762465) q[2];
sx q[2];
rz(-1.1156811) q[2];
sx q[2];
rz(-2.3836977) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75569713) q[1];
sx q[1];
rz(-0.79078005) q[1];
sx q[1];
rz(0.62259953) q[1];
rz(-1.3128223) q[3];
sx q[3];
rz(-1.3822864) q[3];
sx q[3];
rz(0.21125716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4215309) q[2];
sx q[2];
rz(-1.820182) q[2];
sx q[2];
rz(2.9113972) q[2];
rz(0.0097533334) q[3];
sx q[3];
rz(-1.6246656) q[3];
sx q[3];
rz(-2.9856288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4474051) q[0];
sx q[0];
rz(-1.367584) q[0];
sx q[0];
rz(0.76960027) q[0];
rz(-2.1790478) q[1];
sx q[1];
rz(-1.3178408) q[1];
sx q[1];
rz(-2.1544971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55961032) q[0];
sx q[0];
rz(-1.501069) q[0];
sx q[0];
rz(2.0980673) q[0];
rz(-pi) q[1];
rz(1.0133621) q[2];
sx q[2];
rz(-2.8575572) q[2];
sx q[2];
rz(1.3576053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1137352) q[1];
sx q[1];
rz(-0.052885508) q[1];
sx q[1];
rz(-0.89370339) q[1];
x q[2];
rz(-2.3651667) q[3];
sx q[3];
rz(-1.2033495) q[3];
sx q[3];
rz(2.020379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4276245) q[2];
sx q[2];
rz(-2.2043113) q[2];
sx q[2];
rz(3.0734708) q[2];
rz(2.6918329) q[3];
sx q[3];
rz(-2.7263548) q[3];
sx q[3];
rz(1.4382039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.149067) q[0];
sx q[0];
rz(-1.608792) q[0];
sx q[0];
rz(1.7444862) q[0];
rz(0.29009157) q[1];
sx q[1];
rz(-0.42872226) q[1];
sx q[1];
rz(0.52679481) q[1];
rz(-2.6424825) q[2];
sx q[2];
rz(-1.8815831) q[2];
sx q[2];
rz(-0.56167233) q[2];
rz(-2.7431106) q[3];
sx q[3];
rz(-1.5881817) q[3];
sx q[3];
rz(2.6276799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
