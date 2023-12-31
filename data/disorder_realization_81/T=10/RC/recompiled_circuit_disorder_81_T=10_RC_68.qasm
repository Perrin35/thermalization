OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(1.3786432) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(4.6255914) q[1];
sx q[1];
rz(8.9738823) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3987797) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(0.16422693) q[0];
rz(-pi) q[1];
rz(0.11159201) q[2];
sx q[2];
rz(-1.0942232) q[2];
sx q[2];
rz(3.1389719) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7887468) q[1];
sx q[1];
rz(-0.55410085) q[1];
sx q[1];
rz(-0.43177859) q[1];
x q[2];
rz(1.6879184) q[3];
sx q[3];
rz(-2.4544567) q[3];
sx q[3];
rz(0.62108921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1575872) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(-2.297304) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(2.5355693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5490897) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(0.26309183) q[0];
rz(0.94353765) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(1.1862322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77957905) q[0];
sx q[0];
rz(-1.5148666) q[0];
sx q[0];
rz(-1.5520142) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5086753) q[2];
sx q[2];
rz(-1.3717692) q[2];
sx q[2];
rz(-1.8387427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1520878) q[1];
sx q[1];
rz(-2.4651335) q[1];
sx q[1];
rz(-1.3701887) q[1];
x q[2];
rz(1.3783185) q[3];
sx q[3];
rz(-0.88497439) q[3];
sx q[3];
rz(1.6845077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(1.1594695) q[2];
rz(-0.37108478) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3804669) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(-2.3348715) q[0];
rz(-2.9280248) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0778753) q[0];
sx q[0];
rz(-2.2532007) q[0];
sx q[0];
rz(2.9966485) q[0];
rz(-pi) q[1];
rz(-2.9224612) q[2];
sx q[2];
rz(-1.0872772) q[2];
sx q[2];
rz(-0.52106524) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.080938235) q[1];
sx q[1];
rz(-1.9229691) q[1];
sx q[1];
rz(-2.0209795) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0355755) q[3];
sx q[3];
rz(-1.6764063) q[3];
sx q[3];
rz(-1.6184023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.31072581) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-0.93079981) q[2];
rz(-0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(-2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(0.43831929) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(-1.320425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60162773) q[0];
sx q[0];
rz(-0.40973445) q[0];
sx q[0];
rz(-1.5967303) q[0];
rz(-2.690372) q[2];
sx q[2];
rz(-1.3845452) q[2];
sx q[2];
rz(-2.1860683) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1530694) q[1];
sx q[1];
rz(-1.3077902) q[1];
sx q[1];
rz(2.4371229) q[1];
rz(0.78762357) q[3];
sx q[3];
rz(-2.0260603) q[3];
sx q[3];
rz(0.61592197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-0.92869174) q[2];
sx q[2];
rz(2.7992115) q[2];
rz(0.17677447) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(-2.0006196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(1.226549) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(-1.3006166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9601701) q[0];
sx q[0];
rz(-1.9248065) q[0];
sx q[0];
rz(-1.5976853) q[0];
x q[1];
rz(1.498921) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(2.9938811) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3805441) q[1];
sx q[1];
rz(-2.492978) q[1];
sx q[1];
rz(-2.5785239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6404068) q[3];
sx q[3];
rz(-0.29705829) q[3];
sx q[3];
rz(0.81851573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(2.664393) q[2];
rz(0.19208433) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(-0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79648298) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(-3.1298424) q[0];
rz(-2.5911962) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(1.5884429) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5396744) q[0];
sx q[0];
rz(-1.9214905) q[0];
sx q[0];
rz(0.35777103) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9366829) q[2];
sx q[2];
rz(-1.5503746) q[2];
sx q[2];
rz(0.4991971) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8223871) q[1];
sx q[1];
rz(-0.86383312) q[1];
sx q[1];
rz(1.6824526) q[1];
x q[2];
rz(1.4828959) q[3];
sx q[3];
rz(-0.70471901) q[3];
sx q[3];
rz(0.3375012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.50756303) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(-1.8590415) q[2];
rz(-1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(-1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(-0.67725956) q[0];
rz(0.15180763) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(-2.1645434) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0316232) q[0];
sx q[0];
rz(-0.97951802) q[0];
sx q[0];
rz(1.6615608) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5708959) q[2];
sx q[2];
rz(-1.7028371) q[2];
sx q[2];
rz(-0.11167234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.2077142) q[1];
sx q[1];
rz(-1.9095451) q[1];
sx q[1];
rz(-3.0599041) q[1];
rz(-1.8924176) q[3];
sx q[3];
rz(-1.8971895) q[3];
sx q[3];
rz(1.4332989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7523505) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(2.1288669) q[2];
rz(1.1879454) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1241207) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(-2.0195122) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(1.2493856) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69985336) q[0];
sx q[0];
rz(-2.925736) q[0];
sx q[0];
rz(2.9237843) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0935358) q[2];
sx q[2];
rz(-1.153423) q[2];
sx q[2];
rz(-1.0440895) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5640806) q[1];
sx q[1];
rz(-3.1187594) q[1];
sx q[1];
rz(-0.24502416) q[1];
rz(-pi) q[2];
rz(0.28835339) q[3];
sx q[3];
rz(-0.95803146) q[3];
sx q[3];
rz(-1.8307277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32968783) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(-1.1784941) q[2];
rz(1.684749) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6417398) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(-1.8485803) q[0];
rz(1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(-0.59757772) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7062089) q[0];
sx q[0];
rz(-1.7058813) q[0];
sx q[0];
rz(0.038890966) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75002807) q[2];
sx q[2];
rz(-2.4366597) q[2];
sx q[2];
rz(-2.7761369) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5850726) q[1];
sx q[1];
rz(-0.44610281) q[1];
sx q[1];
rz(0.82533605) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4937917) q[3];
sx q[3];
rz(-2.1553851) q[3];
sx q[3];
rz(1.0494174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.22275816) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(1.9082327) q[2];
rz(-2.2402066) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(-1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047886588) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(3.1179324) q[0];
rz(-0.95611447) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(0.67217174) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8925079) q[0];
sx q[0];
rz(-1.5605643) q[0];
sx q[0];
rz(-3.1363048) q[0];
rz(-pi) q[1];
rz(1.5800843) q[2];
sx q[2];
rz(-2.0822968) q[2];
sx q[2];
rz(-1.4052504) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8476665) q[1];
sx q[1];
rz(-0.77984174) q[1];
sx q[1];
rz(-0.49077175) q[1];
rz(-1.724733) q[3];
sx q[3];
rz(-1.142475) q[3];
sx q[3];
rz(-2.0314914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6293634) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(0.36995861) q[2];
rz(-1.6379179) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(-1.2009719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(-0.72369408) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(-2.9818515) q[2];
sx q[2];
rz(-1.1764871) q[2];
sx q[2];
rz(2.7574678) q[2];
rz(-1.9772114) q[3];
sx q[3];
rz(-1.2953399) q[3];
sx q[3];
rz(-3.0084707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
