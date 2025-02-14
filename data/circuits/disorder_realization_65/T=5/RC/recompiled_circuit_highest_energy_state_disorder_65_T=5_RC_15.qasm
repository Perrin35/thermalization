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
rz(1.7841568) q[0];
sx q[0];
rz(-0.62499243) q[0];
sx q[0];
rz(-1.619119) q[0];
rz(-1.9982665) q[1];
sx q[1];
rz(-0.63653094) q[1];
sx q[1];
rz(-1.7826537) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5769437) q[0];
sx q[0];
rz(-1.2822106) q[0];
sx q[0];
rz(2.1493069) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.236249) q[2];
sx q[2];
rz(-2.6839089) q[2];
sx q[2];
rz(-2.8913218) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.57558838) q[1];
sx q[1];
rz(-1.8270565) q[1];
sx q[1];
rz(2.9017067) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0380958) q[3];
sx q[3];
rz(-2.0465133) q[3];
sx q[3];
rz(0.58799839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4437359) q[2];
sx q[2];
rz(-0.82252684) q[2];
sx q[2];
rz(-2.6653384) q[2];
rz(-3.0016628) q[3];
sx q[3];
rz(-1.5014476) q[3];
sx q[3];
rz(0.073277624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78854617) q[0];
sx q[0];
rz(-1.6206425) q[0];
sx q[0];
rz(0.36343685) q[0];
rz(2.4711171) q[1];
sx q[1];
rz(-1.8164219) q[1];
sx q[1];
rz(-1.0796116) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9417861) q[0];
sx q[0];
rz(-0.58871709) q[0];
sx q[0];
rz(-2.0090282) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7893546) q[2];
sx q[2];
rz(-2.5982476) q[2];
sx q[2];
rz(-1.8107514) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.55077165) q[1];
sx q[1];
rz(-1.6147922) q[1];
sx q[1];
rz(1.3293056) q[1];
rz(-pi) q[2];
rz(1.8325975) q[3];
sx q[3];
rz(-2.1356694) q[3];
sx q[3];
rz(1.3606646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67538992) q[2];
sx q[2];
rz(-2.4066996) q[2];
sx q[2];
rz(0.43851635) q[2];
rz(-2.130326) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(-1.6605759) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587875) q[0];
sx q[0];
rz(-2.6130982) q[0];
sx q[0];
rz(-0.82861376) q[0];
rz(1.8939182) q[1];
sx q[1];
rz(-1.9745461) q[1];
sx q[1];
rz(-2.8819328) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22406507) q[0];
sx q[0];
rz(-2.4026786) q[0];
sx q[0];
rz(-1.9784443) q[0];
x q[1];
rz(-2.3724062) q[2];
sx q[2];
rz(-2.63317) q[2];
sx q[2];
rz(-1.891013) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.34777113) q[1];
sx q[1];
rz(-2.035537) q[1];
sx q[1];
rz(-0.79885599) q[1];
x q[2];
rz(-2.8494495) q[3];
sx q[3];
rz(-2.202987) q[3];
sx q[3];
rz(1.2911673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9344249) q[2];
sx q[2];
rz(-1.577689) q[2];
sx q[2];
rz(0.820532) q[2];
rz(-2.3151243) q[3];
sx q[3];
rz(-0.99244899) q[3];
sx q[3];
rz(1.3124527) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8823223) q[0];
sx q[0];
rz(-2.3311908) q[0];
sx q[0];
rz(0.93233863) q[0];
rz(-1.3177634) q[1];
sx q[1];
rz(-1.980314) q[1];
sx q[1];
rz(-0.12277776) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6500191) q[0];
sx q[0];
rz(-0.54560018) q[0];
sx q[0];
rz(-2.9789717) q[0];
rz(-0.81715314) q[2];
sx q[2];
rz(-1.7981525) q[2];
sx q[2];
rz(2.0968693) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.40298324) q[1];
sx q[1];
rz(-1.5718286) q[1];
sx q[1];
rz(-2.3139075) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6969321) q[3];
sx q[3];
rz(-2.529699) q[3];
sx q[3];
rz(-1.9652308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54470283) q[2];
sx q[2];
rz(-2.2201316) q[2];
sx q[2];
rz(0.60453647) q[2];
rz(2.4321411) q[3];
sx q[3];
rz(-0.76990288) q[3];
sx q[3];
rz(-0.82367045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17968793) q[0];
sx q[0];
rz(-0.12374319) q[0];
sx q[0];
rz(1.3667579) q[0];
rz(2.1429515) q[1];
sx q[1];
rz(-0.81592453) q[1];
sx q[1];
rz(-0.5407812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3546412) q[0];
sx q[0];
rz(-1.4313414) q[0];
sx q[0];
rz(-2.9324173) q[0];
rz(-pi) q[1];
rz(2.2774463) q[2];
sx q[2];
rz(-2.1808778) q[2];
sx q[2];
rz(0.44877258) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65366983) q[1];
sx q[1];
rz(-1.8663632) q[1];
sx q[1];
rz(-1.771677) q[1];
x q[2];
rz(-0.56638797) q[3];
sx q[3];
rz(-1.4175804) q[3];
sx q[3];
rz(2.9569616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46405408) q[2];
sx q[2];
rz(-2.3276734) q[2];
sx q[2];
rz(2.3133254) q[2];
rz(-1.6895435) q[3];
sx q[3];
rz(-1.4938846) q[3];
sx q[3];
rz(2.233708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0238817) q[0];
sx q[0];
rz(-1.1957059) q[0];
sx q[0];
rz(1.1718132) q[0];
rz(-0.68012971) q[1];
sx q[1];
rz(-1.9155733) q[1];
sx q[1];
rz(2.5853058) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9754575) q[0];
sx q[0];
rz(-1.8248471) q[0];
sx q[0];
rz(-1.182876) q[0];
x q[1];
rz(-0.0027782253) q[2];
sx q[2];
rz(-1.7173212) q[2];
sx q[2];
rz(-2.3511982) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.43832512) q[1];
sx q[1];
rz(-2.5653953) q[1];
sx q[1];
rz(-0.11771113) q[1];
rz(-pi) q[2];
rz(2.5056434) q[3];
sx q[3];
rz(-1.1750286) q[3];
sx q[3];
rz(-1.1102939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3351626) q[2];
sx q[2];
rz(-1.282629) q[2];
sx q[2];
rz(-0.20509091) q[2];
rz(0.80275503) q[3];
sx q[3];
rz(-1.1057248) q[3];
sx q[3];
rz(0.55454379) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71451181) q[0];
sx q[0];
rz(-1.0733805) q[0];
sx q[0];
rz(2.9141973) q[0];
rz(0.79044509) q[1];
sx q[1];
rz(-2.3643654) q[1];
sx q[1];
rz(0.8367742) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1632501) q[0];
sx q[0];
rz(-2.0866835) q[0];
sx q[0];
rz(1.8229026) q[0];
rz(-1.2143986) q[2];
sx q[2];
rz(-2.1292392) q[2];
sx q[2];
rz(-1.3791549) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62722271) q[1];
sx q[1];
rz(-0.52241507) q[1];
sx q[1];
rz(-2.475513) q[1];
rz(-pi) q[2];
rz(2.3804139) q[3];
sx q[3];
rz(-1.2952779) q[3];
sx q[3];
rz(-0.31106424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5493912) q[2];
sx q[2];
rz(-1.4791919) q[2];
sx q[2];
rz(0.57360348) q[2];
rz(-2.8163689) q[3];
sx q[3];
rz(-2.2461788) q[3];
sx q[3];
rz(0.4044683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0322872) q[0];
sx q[0];
rz(-0.17806299) q[0];
sx q[0];
rz(0.67894116) q[0];
rz(0.13488723) q[1];
sx q[1];
rz(-2.0811681) q[1];
sx q[1];
rz(-1.7178242) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8979864) q[0];
sx q[0];
rz(-0.73315128) q[0];
sx q[0];
rz(-0.85444684) q[0];
rz(-pi) q[1];
rz(-0.55439722) q[2];
sx q[2];
rz(-0.94292484) q[2];
sx q[2];
rz(-2.6822821) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.4383761) q[1];
sx q[1];
rz(-1.5408311) q[1];
sx q[1];
rz(-2.334146) q[1];
rz(-pi) q[2];
rz(0.28935953) q[3];
sx q[3];
rz(-1.123424) q[3];
sx q[3];
rz(-0.20848955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5941102) q[2];
sx q[2];
rz(-0.91071931) q[2];
sx q[2];
rz(-2.9098848) q[2];
rz(2.1314651) q[3];
sx q[3];
rz(-1.2332656) q[3];
sx q[3];
rz(-0.89099425) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6758839) q[0];
sx q[0];
rz(-0.75513419) q[0];
sx q[0];
rz(-2.6278507) q[0];
rz(-1.4328009) q[1];
sx q[1];
rz(-1.865973) q[1];
sx q[1];
rz(-1.5076216) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9378273) q[0];
sx q[0];
rz(-1.4368334) q[0];
sx q[0];
rz(-1.521827) q[0];
rz(-0.25357004) q[2];
sx q[2];
rz(-1.7268983) q[2];
sx q[2];
rz(-1.1461604) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6027939) q[1];
sx q[1];
rz(-1.4414296) q[1];
sx q[1];
rz(-1.6954697) q[1];
x q[2];
rz(1.8225372) q[3];
sx q[3];
rz(-1.5949524) q[3];
sx q[3];
rz(-2.5509953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0214009) q[2];
sx q[2];
rz(-2.2970565) q[2];
sx q[2];
rz(-0.25944844) q[2];
rz(-1.9868959) q[3];
sx q[3];
rz(-0.76609937) q[3];
sx q[3];
rz(1.8416789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65086377) q[0];
sx q[0];
rz(-1.6614953) q[0];
sx q[0];
rz(1.65253) q[0];
rz(2.5015855) q[1];
sx q[1];
rz(-1.0970486) q[1];
sx q[1];
rz(-2.1515813) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2442349) q[0];
sx q[0];
rz(-1.4264297) q[0];
sx q[0];
rz(-1.6598527) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1009388) q[2];
sx q[2];
rz(-0.85950101) q[2];
sx q[2];
rz(2.2180722) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9089936) q[1];
sx q[1];
rz(-0.4505583) q[1];
sx q[1];
rz(1.9480991) q[1];
rz(-1.7933937) q[3];
sx q[3];
rz(-1.5464467) q[3];
sx q[3];
rz(-2.3241732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9330357) q[2];
sx q[2];
rz(-1.2988337) q[2];
sx q[2];
rz(-0.2253069) q[2];
rz(-0.92132583) q[3];
sx q[3];
rz(-0.39042979) q[3];
sx q[3];
rz(-0.050617378) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13023547) q[0];
sx q[0];
rz(-0.36778944) q[0];
sx q[0];
rz(0.854048) q[0];
rz(-1.8232952) q[1];
sx q[1];
rz(-1.5250991) q[1];
sx q[1];
rz(-0.93223882) q[1];
rz(0.097250289) q[2];
sx q[2];
rz(-1.8293867) q[2];
sx q[2];
rz(0.23501227) q[2];
rz(2.1193567) q[3];
sx q[3];
rz(-1.4132186) q[3];
sx q[3];
rz(1.3506387) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
