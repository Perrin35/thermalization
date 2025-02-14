OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5819017) q[0];
sx q[0];
rz(-3.0540967) q[0];
sx q[0];
rz(3.0425332) q[0];
rz(2.8506408) q[1];
sx q[1];
rz(-1.7506316) q[1];
sx q[1];
rz(1.6265534) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6884424) q[0];
sx q[0];
rz(-0.52061284) q[0];
sx q[0];
rz(-2.0474035) q[0];
rz(-0.85364437) q[2];
sx q[2];
rz(-2.5042297) q[2];
sx q[2];
rz(-1.846922) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5196621) q[1];
sx q[1];
rz(-0.77072137) q[1];
sx q[1];
rz(-1.8008485) q[1];
rz(-pi) q[2];
x q[2];
rz(1.280925) q[3];
sx q[3];
rz(-1.8767281) q[3];
sx q[3];
rz(-3.0973697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0809975) q[2];
sx q[2];
rz(-0.011757714) q[2];
sx q[2];
rz(0.11059977) q[2];
rz(-0.009875385) q[3];
sx q[3];
rz(-0.014009352) q[3];
sx q[3];
rz(-0.88147718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.1572384) q[0];
sx q[0];
rz(-0.088033661) q[0];
sx q[0];
rz(1.7617759) q[0];
rz(0.012160483) q[1];
sx q[1];
rz(-0.98406839) q[1];
sx q[1];
rz(-1.5252569) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75689689) q[0];
sx q[0];
rz(-1.8840356) q[0];
sx q[0];
rz(1.5114853) q[0];
rz(0.41864606) q[2];
sx q[2];
rz(-1.3867298) q[2];
sx q[2];
rz(-0.41773645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5489651) q[1];
sx q[1];
rz(-1.4969259) q[1];
sx q[1];
rz(0.17358257) q[1];
x q[2];
rz(0.2539999) q[3];
sx q[3];
rz(-0.45886654) q[3];
sx q[3];
rz(2.7953058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81185594) q[2];
sx q[2];
rz(-3.1168823) q[2];
sx q[2];
rz(-0.032026637) q[2];
rz(-2.5305667) q[3];
sx q[3];
rz(-1.9910087) q[3];
sx q[3];
rz(1.376763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5585612) q[0];
sx q[0];
rz(-0.91201454) q[0];
sx q[0];
rz(-1.4645905) q[0];
rz(-1.5088082) q[1];
sx q[1];
rz(-1.2631515) q[1];
sx q[1];
rz(0.060922932) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39696844) q[0];
sx q[0];
rz(-1.3875657) q[0];
sx q[0];
rz(-2.7559962) q[0];
rz(-pi) q[1];
rz(0.82229422) q[2];
sx q[2];
rz(-0.15317433) q[2];
sx q[2];
rz(2.8399452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12184252) q[1];
sx q[1];
rz(-1.3977244) q[1];
sx q[1];
rz(-2.9815841) q[1];
x q[2];
rz(0.64843602) q[3];
sx q[3];
rz(-0.47659527) q[3];
sx q[3];
rz(-0.70515448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9691201) q[2];
sx q[2];
rz(-0.0027593297) q[2];
sx q[2];
rz(-0.74392444) q[2];
rz(2.7562691) q[3];
sx q[3];
rz(-0.0016366882) q[3];
sx q[3];
rz(-2.2374432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0124403) q[0];
sx q[0];
rz(-0.04239447) q[0];
sx q[0];
rz(2.2989035) q[0];
rz(2.1878302) q[1];
sx q[1];
rz(-1.6396435) q[1];
sx q[1];
rz(-1.5927673) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020352954) q[0];
sx q[0];
rz(-1.5505393) q[0];
sx q[0];
rz(0.11718226) q[0];
x q[1];
rz(1.1723123) q[2];
sx q[2];
rz(-1.7165802) q[2];
sx q[2];
rz(-0.1383986) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5332414) q[1];
sx q[1];
rz(-0.63801879) q[1];
sx q[1];
rz(1.649943) q[1];
rz(-pi) q[2];
rz(1.6706659) q[3];
sx q[3];
rz(-2.348987) q[3];
sx q[3];
rz(0.85505337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.10161764) q[2];
sx q[2];
rz(-0.016923252) q[2];
sx q[2];
rz(-2.7162111) q[2];
rz(-1.3500805) q[3];
sx q[3];
rz(-1.5451508) q[3];
sx q[3];
rz(-2.8369956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199583) q[0];
sx q[0];
rz(-0.0036792734) q[0];
sx q[0];
rz(-0.71596181) q[0];
rz(1.6250027) q[1];
sx q[1];
rz(-0.45418987) q[1];
sx q[1];
rz(-2.9862278) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5276521) q[0];
sx q[0];
rz(-3.0342073) q[0];
sx q[0];
rz(-1.7201109) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22414872) q[2];
sx q[2];
rz(-2.756065) q[2];
sx q[2];
rz(1.7412501) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6581056) q[1];
sx q[1];
rz(-2.8466821) q[1];
sx q[1];
rz(-0.43842478) q[1];
x q[2];
rz(-0.89837667) q[3];
sx q[3];
rz(-1.2247318) q[3];
sx q[3];
rz(-2.4806674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4510437) q[2];
sx q[2];
rz(-2.0205708) q[2];
sx q[2];
rz(-1.5345908) q[2];
rz(2.0524041) q[3];
sx q[3];
rz(-0.023622731) q[3];
sx q[3];
rz(-1.7781809) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4264193) q[0];
sx q[0];
rz(-0.024113163) q[0];
sx q[0];
rz(0.56060785) q[0];
rz(0.86588612) q[1];
sx q[1];
rz(-0.0090323369) q[1];
sx q[1];
rz(0.83585709) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7880121) q[0];
sx q[0];
rz(-2.8451974) q[0];
sx q[0];
rz(-0.061876492) q[0];
x q[1];
rz(3.1070527) q[2];
sx q[2];
rz(-1.1636834) q[2];
sx q[2];
rz(-1.7742334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2627203) q[1];
sx q[1];
rz(-1.6769028) q[1];
sx q[1];
rz(-0.2596007) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3243598) q[3];
sx q[3];
rz(-1.2643688) q[3];
sx q[3];
rz(-2.5401194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5777638) q[2];
sx q[2];
rz(-1.3059629) q[2];
sx q[2];
rz(-1.5888265) q[2];
rz(1.0519489) q[3];
sx q[3];
rz(-2.6237539) q[3];
sx q[3];
rz(2.6036105) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2688667) q[0];
sx q[0];
rz(-2.5476542) q[0];
sx q[0];
rz(-1.3237704) q[0];
rz(-0.84953228) q[1];
sx q[1];
rz(-0.001566611) q[1];
sx q[1];
rz(2.3193287) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9233775) q[0];
sx q[0];
rz(-1.4750622) q[0];
sx q[0];
rz(-1.3908007) q[0];
x q[1];
rz(0.15531628) q[2];
sx q[2];
rz(-0.34201038) q[2];
sx q[2];
rz(1.6928455) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2660127) q[1];
sx q[1];
rz(-3.071738) q[1];
sx q[1];
rz(1.6641947) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5093171) q[3];
sx q[3];
rz(-2.816449) q[3];
sx q[3];
rz(-2.7141311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5620586) q[2];
sx q[2];
rz(-1.2349393) q[2];
sx q[2];
rz(0.84260064) q[2];
rz(0.54251999) q[3];
sx q[3];
rz(-0.018535651) q[3];
sx q[3];
rz(2.5491469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056996718) q[0];
sx q[0];
rz(-1.5472941) q[0];
sx q[0];
rz(2.090276) q[0];
rz(1.7683138) q[1];
sx q[1];
rz(-3.1128502) q[1];
sx q[1];
rz(-1.9017259) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1418041) q[0];
sx q[0];
rz(-1.8134789) q[0];
sx q[0];
rz(1.61116) q[0];
x q[1];
rz(1.6592024) q[2];
sx q[2];
rz(-2.7475221) q[2];
sx q[2];
rz(1.5752156) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.925773) q[1];
sx q[1];
rz(-2.5160393) q[1];
sx q[1];
rz(1.9240153) q[1];
rz(0.9628125) q[3];
sx q[3];
rz(-0.32139489) q[3];
sx q[3];
rz(-2.8803745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5956868) q[2];
sx q[2];
rz(-3.1133856) q[2];
sx q[2];
rz(2.9941881) q[2];
rz(1.6029415) q[3];
sx q[3];
rz(-1.7037062) q[3];
sx q[3];
rz(0.182972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.1786757) q[0];
sx q[0];
rz(-0.33483949) q[0];
sx q[0];
rz(-1.3382925) q[0];
rz(-1.6637404) q[1];
sx q[1];
rz(-1.1344974) q[1];
sx q[1];
rz(-0.17325625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4497183) q[0];
sx q[0];
rz(-1.0074179) q[0];
sx q[0];
rz(1.059235) q[0];
x q[1];
rz(2.1384754) q[2];
sx q[2];
rz(-0.67174339) q[2];
sx q[2];
rz(-2.9152108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0782287) q[1];
sx q[1];
rz(-1.3241819) q[1];
sx q[1];
rz(1.4117086) q[1];
rz(-pi) q[2];
rz(-1.9789437) q[3];
sx q[3];
rz(-0.014685304) q[3];
sx q[3];
rz(2.5774308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0080537) q[2];
sx q[2];
rz(-0.011913813) q[2];
sx q[2];
rz(-1.0201721) q[2];
rz(-0.080848761) q[3];
sx q[3];
rz(-3.1404218) q[3];
sx q[3];
rz(-1.1297869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6218277) q[0];
sx q[0];
rz(-2.8563359) q[0];
sx q[0];
rz(1.4813625) q[0];
rz(-0.18352428) q[1];
sx q[1];
rz(-0.32569519) q[1];
sx q[1];
rz(3.0537925) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0073951) q[0];
sx q[0];
rz(-2.0124276) q[0];
sx q[0];
rz(-1.6820119) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8351458) q[2];
sx q[2];
rz(-1.6493091) q[2];
sx q[2];
rz(1.7959838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1488249) q[1];
sx q[1];
rz(-1.449905) q[1];
sx q[1];
rz(0.24538998) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7664639) q[3];
sx q[3];
rz(-2.4371456) q[3];
sx q[3];
rz(2.7337027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8537019) q[2];
sx q[2];
rz(-0.012306865) q[2];
sx q[2];
rz(0.88407174) q[2];
rz(-0.70841241) q[3];
sx q[3];
rz(-3.1413779) q[3];
sx q[3];
rz(-2.8475672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6207599) q[0];
sx q[0];
rz(-1.5694869) q[0];
sx q[0];
rz(-1.6318305) q[0];
rz(0.028342551) q[1];
sx q[1];
rz(-0.62348532) q[1];
sx q[1];
rz(-3.0709406) q[1];
rz(0.6762517) q[2];
sx q[2];
rz(-2.5898503) q[2];
sx q[2];
rz(2.5838625) q[2];
rz(0.36793637) q[3];
sx q[3];
rz(-0.76303861) q[3];
sx q[3];
rz(-1.7581024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
