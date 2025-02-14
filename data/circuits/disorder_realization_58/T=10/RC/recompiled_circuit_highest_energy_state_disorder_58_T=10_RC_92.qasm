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
rz(3.0296221) q[0];
sx q[0];
rz(-1.0267216) q[0];
sx q[0];
rz(-1.8486899) q[0];
rz(-1.0678043) q[1];
sx q[1];
rz(-0.81967241) q[1];
sx q[1];
rz(1.1991062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0584058) q[0];
sx q[0];
rz(-1.633344) q[0];
sx q[0];
rz(3.0330171) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3875628) q[2];
sx q[2];
rz(-2.2361045) q[2];
sx q[2];
rz(2.5380425) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.407302) q[1];
sx q[1];
rz(-1.3774187) q[1];
sx q[1];
rz(2.1466682) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85924826) q[3];
sx q[3];
rz(-1.0785127) q[3];
sx q[3];
rz(-1.2788206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5045515) q[2];
sx q[2];
rz(-1.7402288) q[2];
sx q[2];
rz(-1.6999647) q[2];
rz(-2.1291034) q[3];
sx q[3];
rz(-1.1463405) q[3];
sx q[3];
rz(0.47154021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(0.83990324) q[0];
sx q[0];
rz(-2.58044) q[0];
sx q[0];
rz(1.0294234) q[0];
rz(-1.7816211) q[1];
sx q[1];
rz(-1.1861035) q[1];
sx q[1];
rz(-0.50672466) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8520607) q[0];
sx q[0];
rz(-1.4665876) q[0];
sx q[0];
rz(-3.122217) q[0];
x q[1];
rz(-2.0015745) q[2];
sx q[2];
rz(-0.74249858) q[2];
sx q[2];
rz(3.1378821) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9190127) q[1];
sx q[1];
rz(-1.6363385) q[1];
sx q[1];
rz(2.1497141) q[1];
rz(-0.84552879) q[3];
sx q[3];
rz(-1.7983899) q[3];
sx q[3];
rz(1.9967772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8257251) q[2];
sx q[2];
rz(-0.59810144) q[2];
sx q[2];
rz(-0.21710795) q[2];
rz(1.0202967) q[3];
sx q[3];
rz(-1.812457) q[3];
sx q[3];
rz(-0.98062688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26647907) q[0];
sx q[0];
rz(-1.7784235) q[0];
sx q[0];
rz(-2.3734221) q[0];
rz(-2.8504596) q[1];
sx q[1];
rz(-1.7765287) q[1];
sx q[1];
rz(1.1870144) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6662749) q[0];
sx q[0];
rz(-0.22724477) q[0];
sx q[0];
rz(-0.3792309) q[0];
rz(-1.4974383) q[2];
sx q[2];
rz(-1.5268451) q[2];
sx q[2];
rz(0.32273656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2666013) q[1];
sx q[1];
rz(-0.35540798) q[1];
sx q[1];
rz(0.0796109) q[1];
rz(2.6701982) q[3];
sx q[3];
rz(-2.1389353) q[3];
sx q[3];
rz(-2.527984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27176303) q[2];
sx q[2];
rz(-2.2972079) q[2];
sx q[2];
rz(-0.96770206) q[2];
rz(-2.9076231) q[3];
sx q[3];
rz(-2.440019) q[3];
sx q[3];
rz(-0.35681891) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02803239) q[0];
sx q[0];
rz(-1.5505294) q[0];
sx q[0];
rz(-2.0618942) q[0];
rz(-2.8375541) q[1];
sx q[1];
rz(-1.1540776) q[1];
sx q[1];
rz(1.8255723) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42342351) q[0];
sx q[0];
rz(-1.6945413) q[0];
sx q[0];
rz(2.607671) q[0];
x q[1];
rz(-1.448668) q[2];
sx q[2];
rz(-0.59905648) q[2];
sx q[2];
rz(-1.0775104) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7217654) q[1];
sx q[1];
rz(-2.160795) q[1];
sx q[1];
rz(-0.71728398) q[1];
rz(0.06240978) q[3];
sx q[3];
rz(-1.9963328) q[3];
sx q[3];
rz(2.0151013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5217846) q[2];
sx q[2];
rz(-1.9079756) q[2];
sx q[2];
rz(-3.0002777) q[2];
rz(-2.5610793) q[3];
sx q[3];
rz(-1.6655191) q[3];
sx q[3];
rz(0.22835246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.88437951) q[0];
sx q[0];
rz(-2.3452106) q[0];
sx q[0];
rz(1.4759395) q[0];
rz(-2.2085704) q[1];
sx q[1];
rz(-0.7111744) q[1];
sx q[1];
rz(1.3915541) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0071268) q[0];
sx q[0];
rz(-2.3023805) q[0];
sx q[0];
rz(-0.15739077) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7180842) q[2];
sx q[2];
rz(-1.711634) q[2];
sx q[2];
rz(-0.25885669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7793624) q[1];
sx q[1];
rz(-2.0915151) q[1];
sx q[1];
rz(0.034265072) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0972033) q[3];
sx q[3];
rz(-0.59848753) q[3];
sx q[3];
rz(2.583556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89620227) q[2];
sx q[2];
rz(-0.37084493) q[2];
sx q[2];
rz(-2.89213) q[2];
rz(1.6438515) q[3];
sx q[3];
rz(-1.7353053) q[3];
sx q[3];
rz(0.41516414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6849218) q[0];
sx q[0];
rz(-2.6455854) q[0];
sx q[0];
rz(1.7556835) q[0];
rz(-1.2617525) q[1];
sx q[1];
rz(-1.933814) q[1];
sx q[1];
rz(0.52880803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.38976) q[0];
sx q[0];
rz(-0.44527894) q[0];
sx q[0];
rz(-1.5543429) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98189378) q[2];
sx q[2];
rz(-2.1236651) q[2];
sx q[2];
rz(2.5585563) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7976378) q[1];
sx q[1];
rz(-0.64227637) q[1];
sx q[1];
rz(-1.6828868) q[1];
x q[2];
rz(2.4566023) q[3];
sx q[3];
rz(-1.4695952) q[3];
sx q[3];
rz(1.2736959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1888107) q[2];
sx q[2];
rz(-2.6691801) q[2];
sx q[2];
rz(1.7702276) q[2];
rz(2.93907) q[3];
sx q[3];
rz(-1.6731508) q[3];
sx q[3];
rz(-1.9523581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8093569) q[0];
sx q[0];
rz(-1.3487331) q[0];
sx q[0];
rz(0.15383823) q[0];
rz(2.4065252) q[1];
sx q[1];
rz(-1.091833) q[1];
sx q[1];
rz(2.99248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99054333) q[0];
sx q[0];
rz(-1.4874389) q[0];
sx q[0];
rz(0.91611422) q[0];
rz(-pi) q[1];
rz(2.9682069) q[2];
sx q[2];
rz(-1.2116222) q[2];
sx q[2];
rz(-0.92196354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9184409) q[1];
sx q[1];
rz(-2.0479124) q[1];
sx q[1];
rz(2.3972307) q[1];
rz(-0.090773067) q[3];
sx q[3];
rz(-2.0730264) q[3];
sx q[3];
rz(1.2365562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6814802) q[2];
sx q[2];
rz(-1.6296547) q[2];
sx q[2];
rz(-0.45664772) q[2];
rz(-1.7234507) q[3];
sx q[3];
rz(-2.2599615) q[3];
sx q[3];
rz(-2.452623) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6316471) q[0];
sx q[0];
rz(-2.8896285) q[0];
sx q[0];
rz(3.0889567) q[0];
rz(1.3098199) q[1];
sx q[1];
rz(-2.8927264) q[1];
sx q[1];
rz(-1.2059258) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0154738) q[0];
sx q[0];
rz(-2.0778225) q[0];
sx q[0];
rz(0.4457723) q[0];
x q[1];
rz(-0.3605901) q[2];
sx q[2];
rz(-2.5506488) q[2];
sx q[2];
rz(2.6481368) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5284023) q[1];
sx q[1];
rz(-2.2397016) q[1];
sx q[1];
rz(-1.3701128) q[1];
rz(-pi) q[2];
rz(2.1931529) q[3];
sx q[3];
rz(-1.3788169) q[3];
sx q[3];
rz(0.49302671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2974818) q[2];
sx q[2];
rz(-1.7606807) q[2];
sx q[2];
rz(2.5816176) q[2];
rz(2.0848134) q[3];
sx q[3];
rz(-0.70928514) q[3];
sx q[3];
rz(-2.3287676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63447222) q[0];
sx q[0];
rz(-1.3993323) q[0];
sx q[0];
rz(-1.5647474) q[0];
rz(-1.5722081) q[1];
sx q[1];
rz(-1.9805311) q[1];
sx q[1];
rz(-2.3326468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29306889) q[0];
sx q[0];
rz(-1.5331755) q[0];
sx q[0];
rz(2.6462808) q[0];
rz(-pi) q[1];
rz(-1.9845029) q[2];
sx q[2];
rz(-0.80360694) q[2];
sx q[2];
rz(1.4860019) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8508243) q[1];
sx q[1];
rz(-1.5577258) q[1];
sx q[1];
rz(-2.1816861) q[1];
rz(-pi) q[2];
rz(-1.393717) q[3];
sx q[3];
rz(-1.5054323) q[3];
sx q[3];
rz(-3.1013427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.481015) q[2];
sx q[2];
rz(-0.95260859) q[2];
sx q[2];
rz(-0.048952254) q[2];
rz(1.8598716) q[3];
sx q[3];
rz(-2.6661524) q[3];
sx q[3];
rz(-1.6767282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82185164) q[0];
sx q[0];
rz(-0.26092437) q[0];
sx q[0];
rz(1.9191746) q[0];
rz(1.6659196) q[1];
sx q[1];
rz(-1.2011352) q[1];
sx q[1];
rz(-2.6840721) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0442866) q[0];
sx q[0];
rz(-2.4032018) q[0];
sx q[0];
rz(0.13412688) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0779764) q[2];
sx q[2];
rz(-1.9263144) q[2];
sx q[2];
rz(1.7001482) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19472518) q[1];
sx q[1];
rz(-1.3683142) q[1];
sx q[1];
rz(0.35670403) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4549834) q[3];
sx q[3];
rz(-1.2401738) q[3];
sx q[3];
rz(2.7329684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7484625) q[2];
sx q[2];
rz(-0.25105432) q[2];
sx q[2];
rz(-1.040323) q[2];
rz(0.48804247) q[3];
sx q[3];
rz(-1.4405684) q[3];
sx q[3];
rz(-2.603781) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3671065) q[0];
sx q[0];
rz(-1.0716866) q[0];
sx q[0];
rz(-2.5197784) q[0];
rz(1.294301) q[1];
sx q[1];
rz(-2.558567) q[1];
sx q[1];
rz(-0.8898215) q[1];
rz(-2.1637259) q[2];
sx q[2];
rz(-1.821283) q[2];
sx q[2];
rz(-1.4878426) q[2];
rz(-2.4518012) q[3];
sx q[3];
rz(-1.2648911) q[3];
sx q[3];
rz(-2.779724) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
