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
rz(-2.9057246) q[0];
sx q[0];
rz(-1.9960825) q[0];
sx q[0];
rz(3.0930162) q[0];
rz(1.5161169) q[1];
sx q[1];
rz(4.2344344) q[1];
sx q[1];
rz(8.7090127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5950105) q[0];
sx q[0];
rz(-1.6002185) q[0];
sx q[0];
rz(1.741623) q[0];
rz(-pi) q[1];
rz(-0.20697941) q[2];
sx q[2];
rz(-1.845115) q[2];
sx q[2];
rz(-1.5649318) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6246031) q[1];
sx q[1];
rz(-1.9698242) q[1];
sx q[1];
rz(2.3844403) q[1];
rz(-1.1668901) q[3];
sx q[3];
rz(-2.0527186) q[3];
sx q[3];
rz(-1.9306077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48892659) q[2];
sx q[2];
rz(-2.581614) q[2];
sx q[2];
rz(1.0555335) q[2];
rz(2.7355898) q[3];
sx q[3];
rz(-1.3136274) q[3];
sx q[3];
rz(-0.82160151) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838298) q[0];
sx q[0];
rz(-1.0394179) q[0];
sx q[0];
rz(0.55617547) q[0];
rz(-1.8122541) q[1];
sx q[1];
rz(-0.61873299) q[1];
sx q[1];
rz(-0.083273085) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2912968) q[0];
sx q[0];
rz(-1.2477739) q[0];
sx q[0];
rz(-0.6685386) q[0];
rz(-pi) q[1];
rz(-2.5640798) q[2];
sx q[2];
rz(-1.0339875) q[2];
sx q[2];
rz(2.2978738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.992112) q[1];
sx q[1];
rz(-2.6410823) q[1];
sx q[1];
rz(-0.47187279) q[1];
rz(2.2733685) q[3];
sx q[3];
rz(-2.4463042) q[3];
sx q[3];
rz(0.27588683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.1978153) q[2];
sx q[2];
rz(-1.3894206) q[2];
sx q[2];
rz(0.49276349) q[2];
rz(1.0265776) q[3];
sx q[3];
rz(-0.30288282) q[3];
sx q[3];
rz(1.7290285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76348412) q[0];
sx q[0];
rz(-2.5788827) q[0];
sx q[0];
rz(-0.43877959) q[0];
rz(1.4525061) q[1];
sx q[1];
rz(-2.8197598) q[1];
sx q[1];
rz(-2.7486393) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2344491) q[0];
sx q[0];
rz(-0.4919855) q[0];
sx q[0];
rz(1.6299295) q[0];
rz(-1.5147171) q[2];
sx q[2];
rz(-0.529037) q[2];
sx q[2];
rz(0.38944405) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.90653268) q[1];
sx q[1];
rz(-1.1533914) q[1];
sx q[1];
rz(0.077666186) q[1];
rz(-pi) q[2];
rz(-0.64928253) q[3];
sx q[3];
rz(-1.5074132) q[3];
sx q[3];
rz(1.6398304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68985933) q[2];
sx q[2];
rz(-0.84443513) q[2];
sx q[2];
rz(-2.8680657) q[2];
rz(0.61255974) q[3];
sx q[3];
rz(-1.4140244) q[3];
sx q[3];
rz(2.2851022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81371152) q[0];
sx q[0];
rz(-2.5484945) q[0];
sx q[0];
rz(-2.3408422) q[0];
rz(2.4294991) q[1];
sx q[1];
rz(-2.021603) q[1];
sx q[1];
rz(-0.48316479) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8034604) q[0];
sx q[0];
rz(-0.82809859) q[0];
sx q[0];
rz(-1.5699734) q[0];
rz(-pi) q[1];
rz(-2.4078385) q[2];
sx q[2];
rz(-0.68063191) q[2];
sx q[2];
rz(-1.5524) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9312716) q[1];
sx q[1];
rz(-0.64669791) q[1];
sx q[1];
rz(-3.0937322) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2340868) q[3];
sx q[3];
rz(-0.90970618) q[3];
sx q[3];
rz(-0.31240901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.157865) q[2];
sx q[2];
rz(-2.4309776) q[2];
sx q[2];
rz(1.2845117) q[2];
rz(-2.6537248) q[3];
sx q[3];
rz(-1.4372466) q[3];
sx q[3];
rz(1.454486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.81216413) q[0];
sx q[0];
rz(-0.85947961) q[0];
sx q[0];
rz(-2.676945) q[0];
rz(-0.94841415) q[1];
sx q[1];
rz(-0.83810884) q[1];
sx q[1];
rz(-1.6882187) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4522355) q[0];
sx q[0];
rz(-0.043137155) q[0];
sx q[0];
rz(-0.94170155) q[0];
rz(-0.27957966) q[2];
sx q[2];
rz(-2.381392) q[2];
sx q[2];
rz(-2.3265944) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.067277597) q[1];
sx q[1];
rz(-2.1805873) q[1];
sx q[1];
rz(1.6526755) q[1];
rz(1.0539758) q[3];
sx q[3];
rz(-2.2272416) q[3];
sx q[3];
rz(0.93641989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1571618) q[2];
sx q[2];
rz(-1.0474019) q[2];
sx q[2];
rz(0.39443991) q[2];
rz(-1.9858342) q[3];
sx q[3];
rz(-0.93005669) q[3];
sx q[3];
rz(-1.8331147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30094576) q[0];
sx q[0];
rz(-0.41579682) q[0];
sx q[0];
rz(-1.0619324) q[0];
rz(-0.95147079) q[1];
sx q[1];
rz(-2.2679592) q[1];
sx q[1];
rz(-1.6208614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5439107) q[0];
sx q[0];
rz(-1.8476227) q[0];
sx q[0];
rz(-3.0863161) q[0];
x q[1];
rz(0.64487793) q[2];
sx q[2];
rz(-1.0928003) q[2];
sx q[2];
rz(-2.1410112) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86506063) q[1];
sx q[1];
rz(-0.78462696) q[1];
sx q[1];
rz(2.9097981) q[1];
rz(-0.49455182) q[3];
sx q[3];
rz(-0.44122094) q[3];
sx q[3];
rz(2.0989204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1157397) q[2];
sx q[2];
rz(-1.9151177) q[2];
sx q[2];
rz(0.24570492) q[2];
rz(1.687441) q[3];
sx q[3];
rz(-1.629963) q[3];
sx q[3];
rz(-0.83034849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.65569735) q[0];
sx q[0];
rz(-2.4035154) q[0];
sx q[0];
rz(-0.6947211) q[0];
rz(-0.10063902) q[1];
sx q[1];
rz(-2.1058319) q[1];
sx q[1];
rz(-0.86520854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14510205) q[0];
sx q[0];
rz(-1.1528413) q[0];
sx q[0];
rz(-1.6161902) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1227046) q[2];
sx q[2];
rz(-0.91045083) q[2];
sx q[2];
rz(0.65062269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45564902) q[1];
sx q[1];
rz(-2.5091617) q[1];
sx q[1];
rz(1.3157822) q[1];
x q[2];
rz(1.0194302) q[3];
sx q[3];
rz(-0.97276455) q[3];
sx q[3];
rz(0.74265146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18409099) q[2];
sx q[2];
rz(-0.75347334) q[2];
sx q[2];
rz(-2.2313879) q[2];
rz(-1.6893859) q[3];
sx q[3];
rz(-1.9096749) q[3];
sx q[3];
rz(-0.46428251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24377395) q[0];
sx q[0];
rz(-2.0938566) q[0];
sx q[0];
rz(2.5183103) q[0];
rz(0.15577623) q[1];
sx q[1];
rz(-2.364295) q[1];
sx q[1];
rz(-2.8782841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1415188) q[0];
sx q[0];
rz(-0.21528582) q[0];
sx q[0];
rz(1.5207284) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6313577) q[2];
sx q[2];
rz(-1.6471631) q[2];
sx q[2];
rz(-1.5431984) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7548482) q[1];
sx q[1];
rz(-1.6476908) q[1];
sx q[1];
rz(-2.5645178) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[1];
rz(-2.2498563) q[2];
sx q[2];
rz(-2.8920434) q[2];
sx q[2];
rz(-2.2535394) q[2];
rz(-0.92787162) q[3];
sx q[3];
rz(-1.1279305) q[3];
sx q[3];
rz(0.97308648) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4230098) q[0];
sx q[0];
rz(-2.6319478) q[0];
sx q[0];
rz(0.50073671) q[0];
rz(-1.1205193) q[1];
sx q[1];
rz(-0.78103939) q[1];
sx q[1];
rz(0.80148554) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028541738) q[0];
sx q[0];
rz(-3.0705912) q[0];
sx q[0];
rz(0.027790471) q[0];
rz(1.0455442) q[2];
sx q[2];
rz(-2.0727951) q[2];
sx q[2];
rz(-2.0567304) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5513331) q[1];
sx q[1];
rz(-0.9551183) q[1];
sx q[1];
rz(2.1034296) q[1];
x q[2];
rz(1.8287703) q[3];
sx q[3];
rz(-1.3822864) q[3];
sx q[3];
rz(0.21125716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.72006172) q[2];
sx q[2];
rz(-1.3214107) q[2];
sx q[2];
rz(0.23019543) q[2];
rz(-3.1318393) q[3];
sx q[3];
rz(-1.516927) q[3];
sx q[3];
rz(2.9856288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4474051) q[0];
sx q[0];
rz(-1.367584) q[0];
sx q[0];
rz(-0.76960027) q[0];
rz(0.96254483) q[1];
sx q[1];
rz(-1.8237518) q[1];
sx q[1];
rz(-0.98709551) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0113163) q[0];
sx q[0];
rz(-0.53142953) q[0];
sx q[0];
rz(-1.4328782) q[0];
x q[1];
rz(2.9883698) q[2];
sx q[2];
rz(-1.330687) q[2];
sx q[2];
rz(0.78165141) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.35008123) q[1];
sx q[1];
rz(-1.6120076) q[1];
sx q[1];
rz(-3.1084395) q[1];
x q[2];
rz(2.6392699) q[3];
sx q[3];
rz(-2.2993616) q[3];
sx q[3];
rz(0.80020927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7139682) q[2];
sx q[2];
rz(-2.2043113) q[2];
sx q[2];
rz(3.0734708) q[2];
rz(-0.44975975) q[3];
sx q[3];
rz(-2.7263548) q[3];
sx q[3];
rz(-1.7033887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.149067) q[0];
sx q[0];
rz(-1.5328007) q[0];
sx q[0];
rz(-1.3971064) q[0];
rz(-0.29009157) q[1];
sx q[1];
rz(-2.7128704) q[1];
sx q[1];
rz(-2.6147978) q[1];
rz(-1.9214966) q[2];
sx q[2];
rz(-2.0439707) q[2];
sx q[2];
rz(1.1743152) q[2];
rz(-3.0968128) q[3];
sx q[3];
rz(-2.7427518) q[3];
sx q[3];
rz(1.0156143) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
