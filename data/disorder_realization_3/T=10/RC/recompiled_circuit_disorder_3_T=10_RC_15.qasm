OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(5.8689868) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6611377) q[0];
sx q[0];
rz(-1.9248795) q[0];
sx q[0];
rz(-2.7829091) q[0];
x q[1];
rz(1.4235052) q[2];
sx q[2];
rz(-2.5669332) q[2];
sx q[2];
rz(-1.6342083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.011885) q[1];
sx q[1];
rz(-1.1457448) q[1];
sx q[1];
rz(0.070406291) q[1];
rz(-pi) q[2];
rz(2.8768086) q[3];
sx q[3];
rz(-1.6695108) q[3];
sx q[3];
rz(-2.7717154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2177314) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(0.031575354) q[2];
rz(-1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(-0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8388222) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(-0.1698499) q[0];
rz(2.4376712) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(0.53952113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9233421) q[0];
sx q[0];
rz(-1.5241511) q[0];
sx q[0];
rz(-2.0204861) q[0];
x q[1];
rz(2.4194854) q[2];
sx q[2];
rz(-1.8849843) q[2];
sx q[2];
rz(-0.22382852) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9425548) q[1];
sx q[1];
rz(-2.6032762) q[1];
sx q[1];
rz(-1.946279) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2803454) q[3];
sx q[3];
rz(-1.1303139) q[3];
sx q[3];
rz(-1.7064106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66723055) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(-0.66611755) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(2.6932122) q[0];
rz(-1.7547296) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(2.8853436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019779531) q[0];
sx q[0];
rz(-0.67987961) q[0];
sx q[0];
rz(-0.42827423) q[0];
rz(2.0492378) q[2];
sx q[2];
rz(-1.97176) q[2];
sx q[2];
rz(1.608633) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8032854) q[1];
sx q[1];
rz(-2.3803664) q[1];
sx q[1];
rz(2.0558946) q[1];
rz(-2.933421) q[3];
sx q[3];
rz(-0.18067193) q[3];
sx q[3];
rz(-2.4908623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3391352) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(-0.57470542) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280076) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(-0.73192275) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5607802) q[0];
sx q[0];
rz(-0.73045759) q[0];
sx q[0];
rz(2.3985582) q[0];
rz(0.71728431) q[2];
sx q[2];
rz(-1.6589763) q[2];
sx q[2];
rz(-0.35536534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6432453) q[1];
sx q[1];
rz(-0.32532641) q[1];
sx q[1];
rz(2.494032) q[1];
rz(1.6963523) q[3];
sx q[3];
rz(-1.351965) q[3];
sx q[3];
rz(-0.20727508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(2.990492) q[2];
rz(-2.5949196) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(1.2623825) q[0];
rz(-1.6732015) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(2.343822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.939643) q[0];
sx q[0];
rz(-1.573283) q[0];
sx q[0];
rz(1.6909088) q[0];
x q[1];
rz(1.7902137) q[2];
sx q[2];
rz(-1.8674388) q[2];
sx q[2];
rz(2.9079633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5835727) q[1];
sx q[1];
rz(-0.94545525) q[1];
sx q[1];
rz(-3.0261092) q[1];
rz(-1.9005152) q[3];
sx q[3];
rz(-1.9922678) q[3];
sx q[3];
rz(-0.7082522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34565869) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(0.30203715) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(-0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40925947) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(-1.0844768) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(0.070080431) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2698343) q[0];
sx q[0];
rz(-1.4176798) q[0];
sx q[0];
rz(1.397875) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5885133) q[2];
sx q[2];
rz(-0.86420176) q[2];
sx q[2];
rz(-1.7318219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6447883) q[1];
sx q[1];
rz(-1.7323238) q[1];
sx q[1];
rz(-2.5700388) q[1];
x q[2];
rz(-1.4962247) q[3];
sx q[3];
rz(-1.6061022) q[3];
sx q[3];
rz(-1.1535742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30248102) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(2.5202259) q[2];
rz(1.7012043) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(-0.85987464) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(-3.133657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36695776) q[0];
sx q[0];
rz(-1.8263706) q[0];
sx q[0];
rz(-0.5704244) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43086149) q[2];
sx q[2];
rz(-0.92509809) q[2];
sx q[2];
rz(0.25260392) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.36564) q[1];
sx q[1];
rz(-2.0740293) q[1];
sx q[1];
rz(0.788049) q[1];
x q[2];
rz(-1.1778529) q[3];
sx q[3];
rz(-1.8624536) q[3];
sx q[3];
rz(1.1100811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(-2.725214) q[2];
rz(1.3683866) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(-2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(-2.7767048) q[0];
rz(-2.2015613) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(-1.5023124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4125) q[0];
sx q[0];
rz(-1.8209429) q[0];
sx q[0];
rz(-1.581122) q[0];
x q[1];
rz(-2.084923) q[2];
sx q[2];
rz(-2.3377315) q[2];
sx q[2];
rz(-2.4664997) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.9099721) q[1];
sx q[1];
rz(-2.2044704) q[1];
sx q[1];
rz(2.2909067) q[1];
rz(-pi) q[2];
rz(1.7765462) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(-2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49729785) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(2.0765182) q[2];
rz(0.30125695) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97312462) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(2.705943) q[0];
rz(1.3849974) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(2.7246144) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57396736) q[0];
sx q[0];
rz(-1.5331393) q[0];
sx q[0];
rz(0.91659878) q[0];
rz(-2.4776393) q[2];
sx q[2];
rz(-1.10154) q[2];
sx q[2];
rz(-2.9479153) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6777842) q[1];
sx q[1];
rz(-0.52075547) q[1];
sx q[1];
rz(1.8522472) q[1];
rz(-pi) q[2];
rz(1.4407518) q[3];
sx q[3];
rz(-1.4308883) q[3];
sx q[3];
rz(-0.035294447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0372662) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(-2.9157675) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(-0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.726783) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.6171932) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(-1.3226002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.246829) q[0];
sx q[0];
rz(-1.6920977) q[0];
sx q[0];
rz(1.9508719) q[0];
rz(-pi) q[1];
rz(-1.1216713) q[2];
sx q[2];
rz(-1.3438864) q[2];
sx q[2];
rz(-2.801148) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1111856) q[1];
sx q[1];
rz(-2.2135995) q[1];
sx q[1];
rz(-0.91099693) q[1];
x q[2];
rz(-1.3721458) q[3];
sx q[3];
rz(-0.66505265) q[3];
sx q[3];
rz(2.8952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7913197) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(2.1255169) q[2];
rz(-1.919205) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14324698) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(-0.017756391) q[2];
sx q[2];
rz(-0.50928558) q[2];
sx q[2];
rz(-1.7094564) q[2];
rz(-1.3397459) q[3];
sx q[3];
rz(-1.6587202) q[3];
sx q[3];
rz(-0.29826577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
