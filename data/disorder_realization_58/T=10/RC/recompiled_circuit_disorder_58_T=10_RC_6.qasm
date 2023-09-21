OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(-2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93847371) q[0];
sx q[0];
rz(-1.9460558) q[0];
sx q[0];
rz(-1.9012326) q[0];
rz(-0.73859282) q[2];
sx q[2];
rz(-1.2054218) q[2];
sx q[2];
rz(1.5473168) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9484529) q[1];
sx q[1];
rz(-1.2558736) q[1];
sx q[1];
rz(0.83955168) q[1];
x q[2];
rz(-2.41483) q[3];
sx q[3];
rz(-0.57075497) q[3];
sx q[3];
rz(-1.4445514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7044907) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(0.1201771) q[2];
rz(-1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(-2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-2.8149014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17900285) q[0];
sx q[0];
rz(-1.0674745) q[0];
sx q[0];
rz(-2.3110564) q[0];
x q[1];
rz(1.1268483) q[2];
sx q[2];
rz(-1.4269281) q[2];
sx q[2];
rz(-0.680188) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.28194004) q[1];
sx q[1];
rz(-1.9743866) q[1];
sx q[1];
rz(-0.90359009) q[1];
rz(3.0558415) q[3];
sx q[3];
rz(-0.88480703) q[3];
sx q[3];
rz(-0.12967295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5102753) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(-3.0316947) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(-1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54863769) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(2.864481) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(2.0842016) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43305106) q[0];
sx q[0];
rz(-1.2305224) q[0];
sx q[0];
rz(0.021854594) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8003642) q[2];
sx q[2];
rz(-1.7592906) q[2];
sx q[2];
rz(0.54268062) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1765269) q[1];
sx q[1];
rz(-0.9170734) q[1];
sx q[1];
rz(2.7214126) q[1];
x q[2];
rz(-1.7626761) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(2.2172745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(-1.7791629) q[2];
rz(2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3574922) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(-2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-2.9761956) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20277682) q[0];
sx q[0];
rz(-1.6949777) q[0];
sx q[0];
rz(1.3846272) q[0];
x q[1];
rz(2.9781614) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(-0.93726678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0796697) q[1];
sx q[1];
rz(-1.8975782) q[1];
sx q[1];
rz(-0.57481874) q[1];
rz(-pi) q[2];
rz(-0.27407077) q[3];
sx q[3];
rz(-2.1087397) q[3];
sx q[3];
rz(-0.9243954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7455204) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(0.20544927) q[2];
rz(-1.127634) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(-1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(-2.7868295) q[0];
rz(-1.9873437) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-0.23194557) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97868279) q[0];
sx q[0];
rz(-1.9918348) q[0];
sx q[0];
rz(1.0324423) q[0];
x q[1];
rz(-2.9734128) q[2];
sx q[2];
rz(-0.7025223) q[2];
sx q[2];
rz(-1.0934193) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9310589) q[1];
sx q[1];
rz(-1.7909044) q[1];
sx q[1];
rz(1.7935497) q[1];
rz(-pi) q[2];
rz(2.7730745) q[3];
sx q[3];
rz(-1.5526062) q[3];
sx q[3];
rz(0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.140124) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(1.9011964) q[2];
rz(-2.5455348) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.8381455) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(2.1441377) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31842445) q[0];
sx q[0];
rz(-2.423375) q[0];
sx q[0];
rz(1.7757925) q[0];
rz(2.4100254) q[2];
sx q[2];
rz(-1.288207) q[2];
sx q[2];
rz(0.76991316) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3012078) q[1];
sx q[1];
rz(-1.9921229) q[1];
sx q[1];
rz(-2.8788484) q[1];
x q[2];
rz(-1.0783844) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(-1.9099964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9138907) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(0.4513936) q[2];
rz(-0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(-1.2020948) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.30615) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(0.62414449) q[0];
rz(1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(0.61378941) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1873338) q[0];
sx q[0];
rz(-1.7566534) q[0];
sx q[0];
rz(-2.3128187) q[0];
x q[1];
rz(2.8579312) q[2];
sx q[2];
rz(-1.5093056) q[2];
sx q[2];
rz(-0.21360699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9871414) q[1];
sx q[1];
rz(-0.11949355) q[1];
sx q[1];
rz(0.51388545) q[1];
x q[2];
rz(-3.0244163) q[3];
sx q[3];
rz(-2.0332608) q[3];
sx q[3];
rz(-1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7081786) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(-1.9667352) q[2];
rz(2.5332149) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(1.7181989) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90010086) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(0.4883782) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(2.1571295) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2135293) q[0];
sx q[0];
rz(-0.54740471) q[0];
sx q[0];
rz(-0.071896032) q[0];
rz(-pi) q[1];
rz(-1.1659053) q[2];
sx q[2];
rz(-2.3003909) q[2];
sx q[2];
rz(-1.0798432) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.24482803) q[1];
sx q[1];
rz(-2.7762189) q[1];
sx q[1];
rz(0.016896292) q[1];
x q[2];
rz(-0.47655388) q[3];
sx q[3];
rz(-1.8503975) q[3];
sx q[3];
rz(-2.9001146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-0.53517503) q[2];
rz(-2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64131367) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(-2.4556659) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(-2.2156782) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6322964) q[0];
sx q[0];
rz(-0.25521454) q[0];
sx q[0];
rz(2.1584828) q[0];
rz(0.8530059) q[2];
sx q[2];
rz(-0.91663137) q[2];
sx q[2];
rz(2.5329563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33854252) q[1];
sx q[1];
rz(-2.0085137) q[1];
sx q[1];
rz(2.9796757) q[1];
rz(-pi) q[2];
rz(0.20415281) q[3];
sx q[3];
rz(-1.4545868) q[3];
sx q[3];
rz(-1.0128563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(-1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(1.013247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17908827) q[0];
sx q[0];
rz(-2.8694186) q[0];
sx q[0];
rz(-1.0286691) q[0];
x q[1];
rz(-0.35372325) q[2];
sx q[2];
rz(-1.612066) q[2];
sx q[2];
rz(1.8281787) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3328247) q[1];
sx q[1];
rz(-1.9908394) q[1];
sx q[1];
rz(1.4648449) q[1];
rz(-pi) q[2];
rz(-0.69443955) q[3];
sx q[3];
rz(-2.246292) q[3];
sx q[3];
rz(1.465786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(0.75941336) q[2];
rz(1.7761207) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(-2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(-2.5683174) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(3.0647031) q[2];
sx q[2];
rz(-1.1270317) q[2];
sx q[2];
rz(1.8539853) q[2];
rz(0.063354062) q[3];
sx q[3];
rz(-0.8739211) q[3];
sx q[3];
rz(-2.7779761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];