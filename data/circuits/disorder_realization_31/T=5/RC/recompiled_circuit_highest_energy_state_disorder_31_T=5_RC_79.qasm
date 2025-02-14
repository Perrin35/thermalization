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
rz(1.2827058) q[0];
sx q[0];
rz(1.016322) q[0];
sx q[0];
rz(10.680351) q[0];
rz(-1.3610871) q[1];
sx q[1];
rz(-0.95870107) q[1];
sx q[1];
rz(0.60671848) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17598303) q[0];
sx q[0];
rz(-2.2874766) q[0];
sx q[0];
rz(0.59881439) q[0];
rz(0.032564596) q[2];
sx q[2];
rz(-1.3059907) q[2];
sx q[2];
rz(2.3015568) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7156034) q[1];
sx q[1];
rz(-1.2066325) q[1];
sx q[1];
rz(-3.0028572) q[1];
rz(1.6912724) q[3];
sx q[3];
rz(-1.7306149) q[3];
sx q[3];
rz(-0.32630703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2029734) q[2];
sx q[2];
rz(-1.1569269) q[2];
sx q[2];
rz(-0.17091664) q[2];
rz(-2.0140698) q[3];
sx q[3];
rz(-0.5564965) q[3];
sx q[3];
rz(3.0881622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.663986) q[0];
sx q[0];
rz(-2.8150616) q[0];
sx q[0];
rz(2.4531181) q[0];
rz(1.7636048) q[1];
sx q[1];
rz(-2.1207899) q[1];
sx q[1];
rz(-0.14150208) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6571638) q[0];
sx q[0];
rz(-1.9701029) q[0];
sx q[0];
rz(2.5998373) q[0];
x q[1];
rz(-0.21464084) q[2];
sx q[2];
rz(-0.43834375) q[2];
sx q[2];
rz(0.12210309) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0558171) q[1];
sx q[1];
rz(-1.2165099) q[1];
sx q[1];
rz(1.8721168) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0612437) q[3];
sx q[3];
rz(-1.1871205) q[3];
sx q[3];
rz(1.2017045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3839533) q[2];
sx q[2];
rz(-2.6944104) q[2];
sx q[2];
rz(-0.028707061) q[2];
rz(2.7590397) q[3];
sx q[3];
rz(-1.5111204) q[3];
sx q[3];
rz(2.9401275) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460227) q[0];
sx q[0];
rz(-2.0992794) q[0];
sx q[0];
rz(0.37164715) q[0];
rz(2.9314575) q[1];
sx q[1];
rz(-2.7919283) q[1];
sx q[1];
rz(1.9336611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2164566) q[0];
sx q[0];
rz(-1.0724004) q[0];
sx q[0];
rz(-0.67277661) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5919331) q[2];
sx q[2];
rz(-1.653228) q[2];
sx q[2];
rz(1.3940514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.74151285) q[1];
sx q[1];
rz(-0.99435213) q[1];
sx q[1];
rz(2.4680572) q[1];
rz(-3.0284027) q[3];
sx q[3];
rz(-2.7940911) q[3];
sx q[3];
rz(0.45873911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5189884) q[2];
sx q[2];
rz(-3.0486139) q[2];
sx q[2];
rz(2.5034215) q[2];
rz(-2.7291164) q[3];
sx q[3];
rz(-1.6715489) q[3];
sx q[3];
rz(2.0392058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0647122) q[0];
sx q[0];
rz(-0.91016722) q[0];
sx q[0];
rz(1.3077211) q[0];
rz(0.67963302) q[1];
sx q[1];
rz(-2.4145587) q[1];
sx q[1];
rz(1.1839428) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1935496) q[0];
sx q[0];
rz(-1.2468766) q[0];
sx q[0];
rz(0.69635038) q[0];
x q[1];
rz(1.4550736) q[2];
sx q[2];
rz(-1.3365796) q[2];
sx q[2];
rz(0.17143347) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.19540994) q[1];
sx q[1];
rz(-1.2210746) q[1];
sx q[1];
rz(-1.6433952) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3700962) q[3];
sx q[3];
rz(-2.960254) q[3];
sx q[3];
rz(2.4490631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.022698553) q[2];
sx q[2];
rz(-0.1476295) q[2];
sx q[2];
rz(-2.0856608) q[2];
rz(-0.28327709) q[3];
sx q[3];
rz(-1.2669468) q[3];
sx q[3];
rz(1.7578846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088575514) q[0];
sx q[0];
rz(-2.179189) q[0];
sx q[0];
rz(-1.4085294) q[0];
rz(-0.22964302) q[1];
sx q[1];
rz(-0.85669986) q[1];
sx q[1];
rz(-1.1913258) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4588683) q[0];
sx q[0];
rz(-2.7153569) q[0];
sx q[0];
rz(-1.2221673) q[0];
rz(-pi) q[1];
rz(1.1896594) q[2];
sx q[2];
rz(-0.32763619) q[2];
sx q[2];
rz(-2.0608847) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3670259) q[1];
sx q[1];
rz(-1.8205678) q[1];
sx q[1];
rz(2.7109409) q[1];
rz(-pi) q[2];
rz(2.5727083) q[3];
sx q[3];
rz(-0.99939049) q[3];
sx q[3];
rz(2.6517154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64629897) q[2];
sx q[2];
rz(-2.7806492) q[2];
sx q[2];
rz(1.4781282) q[2];
rz(-2.9804001) q[3];
sx q[3];
rz(-1.6094004) q[3];
sx q[3];
rz(0.78981367) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6394871) q[0];
sx q[0];
rz(-2.7171071) q[0];
sx q[0];
rz(2.2232527) q[0];
rz(0.25686747) q[1];
sx q[1];
rz(-0.99595064) q[1];
sx q[1];
rz(1.1161944) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74280069) q[0];
sx q[0];
rz(-1.0519807) q[0];
sx q[0];
rz(-0.98720819) q[0];
rz(-pi) q[1];
rz(0.29735469) q[2];
sx q[2];
rz(-2.8562244) q[2];
sx q[2];
rz(-0.076208027) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65016547) q[1];
sx q[1];
rz(-2.8313447) q[1];
sx q[1];
rz(2.4568471) q[1];
rz(-0.61577101) q[3];
sx q[3];
rz(-0.29524657) q[3];
sx q[3];
rz(-2.8930882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18672289) q[2];
sx q[2];
rz(-0.8064417) q[2];
sx q[2];
rz(2.7346129) q[2];
rz(-0.61378971) q[3];
sx q[3];
rz(-1.5301306) q[3];
sx q[3];
rz(-0.49155244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3864022) q[0];
sx q[0];
rz(-2.3004005) q[0];
sx q[0];
rz(-2.2156773) q[0];
rz(0.47261247) q[1];
sx q[1];
rz(-1.0541397) q[1];
sx q[1];
rz(-0.47017631) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42379728) q[0];
sx q[0];
rz(-1.4873624) q[0];
sx q[0];
rz(0.83441894) q[0];
rz(-pi) q[1];
rz(1.0711925) q[2];
sx q[2];
rz(-2.4322369) q[2];
sx q[2];
rz(1.5200523) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2963841) q[1];
sx q[1];
rz(-1.607158) q[1];
sx q[1];
rz(-0.78486376) q[1];
rz(-2.3005465) q[3];
sx q[3];
rz(-1.5435092) q[3];
sx q[3];
rz(0.59852615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4818695) q[2];
sx q[2];
rz(-2.0707776) q[2];
sx q[2];
rz(-0.83079633) q[2];
rz(-0.080549084) q[3];
sx q[3];
rz(-1.0815257) q[3];
sx q[3];
rz(-0.95675937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.047121) q[0];
sx q[0];
rz(-1.6166649) q[0];
sx q[0];
rz(-2.9834874) q[0];
rz(0.72599167) q[1];
sx q[1];
rz(-0.43015614) q[1];
sx q[1];
rz(-2.7714738) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4649432) q[0];
sx q[0];
rz(-1.9029362) q[0];
sx q[0];
rz(1.667883) q[0];
rz(-1.0962531) q[2];
sx q[2];
rz(-2.2703746) q[2];
sx q[2];
rz(-0.84700218) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8565048) q[1];
sx q[1];
rz(-2.5743432) q[1];
sx q[1];
rz(-2.3598266) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1922969) q[3];
sx q[3];
rz(-0.64084508) q[3];
sx q[3];
rz(-2.8261938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.98081723) q[2];
sx q[2];
rz(-1.2951916) q[2];
sx q[2];
rz(-0.13713169) q[2];
rz(-1.49617) q[3];
sx q[3];
rz(-0.6936332) q[3];
sx q[3];
rz(-0.48367286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54336035) q[0];
sx q[0];
rz(-2.1393263) q[0];
sx q[0];
rz(3.015236) q[0];
rz(2.5670746) q[1];
sx q[1];
rz(-0.9205598) q[1];
sx q[1];
rz(0.7472907) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68796038) q[0];
sx q[0];
rz(-1.9575141) q[0];
sx q[0];
rz(2.9403015) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1529601) q[2];
sx q[2];
rz(-0.95541856) q[2];
sx q[2];
rz(-2.5948713) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.14366985) q[1];
sx q[1];
rz(-1.2005245) q[1];
sx q[1];
rz(2.0475494) q[1];
rz(0.46183698) q[3];
sx q[3];
rz(-2.850252) q[3];
sx q[3];
rz(0.1801462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1025461) q[2];
sx q[2];
rz(-2.0206082) q[2];
sx q[2];
rz(0.54753629) q[2];
rz(1.840379) q[3];
sx q[3];
rz(-2.0042714) q[3];
sx q[3];
rz(-3.014452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(1.8319156) q[0];
sx q[0];
rz(-1.3818106) q[0];
sx q[0];
rz(0.58050138) q[0];
rz(0.25513395) q[1];
sx q[1];
rz(-0.83108035) q[1];
sx q[1];
rz(1.1855804) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0806633) q[0];
sx q[0];
rz(-2.1619045) q[0];
sx q[0];
rz(-0.76158686) q[0];
x q[1];
rz(-1.2078778) q[2];
sx q[2];
rz(-1.5722154) q[2];
sx q[2];
rz(-3.0848173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45764582) q[1];
sx q[1];
rz(-2.6630029) q[1];
sx q[1];
rz(-0.19400383) q[1];
rz(-0.2106726) q[3];
sx q[3];
rz(-0.86189358) q[3];
sx q[3];
rz(1.3900932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4127976) q[2];
sx q[2];
rz(-2.2788861) q[2];
sx q[2];
rz(-1.0069138) q[2];
rz(1.6179196) q[3];
sx q[3];
rz(-0.98098743) q[3];
sx q[3];
rz(1.1716918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11378743) q[0];
sx q[0];
rz(-1.8814977) q[0];
sx q[0];
rz(0.30497288) q[0];
rz(1.8995359) q[1];
sx q[1];
rz(-1.8546974) q[1];
sx q[1];
rz(0.7484662) q[1];
rz(0.17340809) q[2];
sx q[2];
rz(-1.1447403) q[2];
sx q[2];
rz(-0.8131342) q[2];
rz(-2.5438944) q[3];
sx q[3];
rz(-2.2168474) q[3];
sx q[3];
rz(-0.26099152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
