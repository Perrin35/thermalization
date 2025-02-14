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
rz(1.500904) q[0];
sx q[0];
rz(-0.7839497) q[0];
sx q[0];
rz(6.1824829) q[0];
rz(1.7164224) q[1];
sx q[1];
rz(-2.4727688) q[1];
sx q[1];
rz(1.5727023) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4948499) q[0];
sx q[0];
rz(-1.6892489) q[0];
sx q[0];
rz(-0.54991566) q[0];
rz(0.35694795) q[2];
sx q[2];
rz(-0.28951807) q[2];
sx q[2];
rz(-3.0384921) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0364337) q[1];
sx q[1];
rz(-1.729125) q[1];
sx q[1];
rz(1.4033844) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80683462) q[3];
sx q[3];
rz(-2.5111622) q[3];
sx q[3];
rz(2.2866319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0147741) q[2];
sx q[2];
rz(-2.7744881) q[2];
sx q[2];
rz(-1.611562) q[2];
rz(-1.4989467) q[3];
sx q[3];
rz(-0.25224125) q[3];
sx q[3];
rz(2.4422395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67795578) q[0];
sx q[0];
rz(-1.3709443) q[0];
sx q[0];
rz(-0.26588765) q[0];
rz(1.6193341) q[1];
sx q[1];
rz(-1.3109173) q[1];
sx q[1];
rz(1.486385) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35039577) q[0];
sx q[0];
rz(-2.3662218) q[0];
sx q[0];
rz(-3.1309137) q[0];
rz(0.3703385) q[2];
sx q[2];
rz(-2.5458434) q[2];
sx q[2];
rz(-2.0496313) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3656908) q[1];
sx q[1];
rz(-1.2925102) q[1];
sx q[1];
rz(2.6551618) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9363181) q[3];
sx q[3];
rz(-2.6418002) q[3];
sx q[3];
rz(-2.5427853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8365525) q[2];
sx q[2];
rz(-1.0779447) q[2];
sx q[2];
rz(-3.0744699) q[2];
rz(-0.51131311) q[3];
sx q[3];
rz(-1.3498243) q[3];
sx q[3];
rz(-0.93222031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.2387282) q[0];
sx q[0];
rz(-1.4816062) q[0];
sx q[0];
rz(0.1486775) q[0];
rz(-2.6909434) q[1];
sx q[1];
rz(-1.5842178) q[1];
sx q[1];
rz(0.91948909) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4212358) q[0];
sx q[0];
rz(-1.8381869) q[0];
sx q[0];
rz(1.8084099) q[0];
rz(0.8834768) q[2];
sx q[2];
rz(-1.7224489) q[2];
sx q[2];
rz(2.6014181) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0120079) q[1];
sx q[1];
rz(-1.1430446) q[1];
sx q[1];
rz(2.5431551) q[1];
rz(-0.33691562) q[3];
sx q[3];
rz(-1.7506699) q[3];
sx q[3];
rz(2.3726187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86842361) q[2];
sx q[2];
rz(-1.0852852) q[2];
sx q[2];
rz(-0.85624179) q[2];
rz(0.87598261) q[3];
sx q[3];
rz(-0.3951422) q[3];
sx q[3];
rz(-1.1538848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3530537) q[0];
sx q[0];
rz(-1.5667229) q[0];
sx q[0];
rz(0.62623155) q[0];
rz(-1.7010472) q[1];
sx q[1];
rz(-2.3150621) q[1];
sx q[1];
rz(0.70762077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2224049) q[0];
sx q[0];
rz(-2.4116431) q[0];
sx q[0];
rz(-1.4428424) q[0];
rz(-pi) q[1];
rz(-1.2265251) q[2];
sx q[2];
rz(-0.81721899) q[2];
sx q[2];
rz(2.5622649) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6615078) q[1];
sx q[1];
rz(-0.28107496) q[1];
sx q[1];
rz(0.29647372) q[1];
rz(-pi) q[2];
rz(-1.6314916) q[3];
sx q[3];
rz(-0.48572054) q[3];
sx q[3];
rz(2.751089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6669199) q[2];
sx q[2];
rz(-1.2771353) q[2];
sx q[2];
rz(-0.39371583) q[2];
rz(0.79931021) q[3];
sx q[3];
rz(-0.36472133) q[3];
sx q[3];
rz(-1.6126311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2213152) q[0];
sx q[0];
rz(-2.3096313) q[0];
sx q[0];
rz(0.078201683) q[0];
rz(-2.4529264) q[1];
sx q[1];
rz(-0.93910256) q[1];
sx q[1];
rz(0.92855531) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79828145) q[0];
sx q[0];
rz(-2.4791299) q[0];
sx q[0];
rz(-2.6210253) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8357115) q[2];
sx q[2];
rz(-1.2291996) q[2];
sx q[2];
rz(-0.46457967) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21338687) q[1];
sx q[1];
rz(-1.7891856) q[1];
sx q[1];
rz(2.7534927) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5220248) q[3];
sx q[3];
rz(-1.2396401) q[3];
sx q[3];
rz(-1.1201657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3277305) q[2];
sx q[2];
rz(-2.6321415) q[2];
sx q[2];
rz(0.65208411) q[2];
rz(-0.70513519) q[3];
sx q[3];
rz(-1.6465681) q[3];
sx q[3];
rz(2.7130073) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.240165) q[0];
sx q[0];
rz(-2.9225898) q[0];
sx q[0];
rz(-1.3448311) q[0];
rz(2.6118458) q[1];
sx q[1];
rz(-1.2836645) q[1];
sx q[1];
rz(1.8240428) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.844125) q[0];
sx q[0];
rz(-0.58124761) q[0];
sx q[0];
rz(-1.7183003) q[0];
x q[1];
rz(2.8805783) q[2];
sx q[2];
rz(-1.3512638) q[2];
sx q[2];
rz(-2.9308776) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3403794) q[1];
sx q[1];
rz(-1.6928343) q[1];
sx q[1];
rz(-1.9122047) q[1];
x q[2];
rz(-1.9002135) q[3];
sx q[3];
rz(-2.6086411) q[3];
sx q[3];
rz(-0.50607133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.162447) q[2];
sx q[2];
rz(-1.2787168) q[2];
sx q[2];
rz(-1.0865967) q[2];
rz(-3.047591) q[3];
sx q[3];
rz(-2.0407929) q[3];
sx q[3];
rz(0.47959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78733665) q[0];
sx q[0];
rz(-0.84699637) q[0];
sx q[0];
rz(-3.060044) q[0];
rz(1.6185919) q[1];
sx q[1];
rz(-1.0146419) q[1];
sx q[1];
rz(0.60752216) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475246) q[0];
sx q[0];
rz(-0.77558625) q[0];
sx q[0];
rz(-1.8344384) q[0];
rz(-pi) q[1];
rz(2.4685988) q[2];
sx q[2];
rz(-0.57306) q[2];
sx q[2];
rz(2.429643) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41740356) q[1];
sx q[1];
rz(-0.10358394) q[1];
sx q[1];
rz(-1.4404907) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19316407) q[3];
sx q[3];
rz(-0.37988362) q[3];
sx q[3];
rz(1.3182321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0532694) q[2];
sx q[2];
rz(-0.20603267) q[2];
sx q[2];
rz(2.6982809) q[2];
rz(-3.0961224) q[3];
sx q[3];
rz(-1.9732405) q[3];
sx q[3];
rz(2.5209691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6670253) q[0];
sx q[0];
rz(-2.1587631) q[0];
sx q[0];
rz(-0.56812754) q[0];
rz(-1.8727632) q[1];
sx q[1];
rz(-1.9867691) q[1];
sx q[1];
rz(2.5297129) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4418714) q[0];
sx q[0];
rz(-1.6590902) q[0];
sx q[0];
rz(-0.37847145) q[0];
x q[1];
rz(2.1558495) q[2];
sx q[2];
rz(-0.32801706) q[2];
sx q[2];
rz(-2.8190921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.001732262) q[1];
sx q[1];
rz(-1.0852768) q[1];
sx q[1];
rz(-0.27487572) q[1];
rz(-1.3723721) q[3];
sx q[3];
rz(-1.4003318) q[3];
sx q[3];
rz(2.9476019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2511217) q[2];
sx q[2];
rz(-1.8939563) q[2];
sx q[2];
rz(-2.9337511) q[2];
rz(-3.1244997) q[3];
sx q[3];
rz(-1.0216252) q[3];
sx q[3];
rz(1.0926854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.51181483) q[0];
sx q[0];
rz(-2.9371174) q[0];
sx q[0];
rz(1.7275607) q[0];
rz(-0.066970197) q[1];
sx q[1];
rz(-1.8194865) q[1];
sx q[1];
rz(-0.13967839) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.139555) q[0];
sx q[0];
rz(-1.25601) q[0];
sx q[0];
rz(-0.33483968) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5508444) q[2];
sx q[2];
rz(-0.70066888) q[2];
sx q[2];
rz(0.76393647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2181704) q[1];
sx q[1];
rz(-1.7814444) q[1];
sx q[1];
rz(2.9107374) q[1];
x q[2];
rz(0.051316694) q[3];
sx q[3];
rz(-1.7398549) q[3];
sx q[3];
rz(-3.0640107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4329873) q[2];
sx q[2];
rz(-1.8726417) q[2];
sx q[2];
rz(-0.73680669) q[2];
rz(-1.6144729) q[3];
sx q[3];
rz(-1.0090642) q[3];
sx q[3];
rz(-1.8710776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0438185) q[0];
sx q[0];
rz(-0.96776861) q[0];
sx q[0];
rz(-1.1678102) q[0];
rz(-0.66520005) q[1];
sx q[1];
rz(-1.2780814) q[1];
sx q[1];
rz(-2.1591689) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8609527) q[0];
sx q[0];
rz(-0.5361983) q[0];
sx q[0];
rz(1.1711189) q[0];
rz(0.63913156) q[2];
sx q[2];
rz(-1.7878727) q[2];
sx q[2];
rz(2.9242196) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1878882) q[1];
sx q[1];
rz(-0.70548234) q[1];
sx q[1];
rz(3.0445552) q[1];
rz(-pi) q[2];
rz(1.1725964) q[3];
sx q[3];
rz(-1.6329582) q[3];
sx q[3];
rz(-0.27649319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3490225) q[2];
sx q[2];
rz(-2.4945365) q[2];
sx q[2];
rz(-1.6046074) q[2];
rz(2.1215306) q[3];
sx q[3];
rz(-1.6198817) q[3];
sx q[3];
rz(-0.1608688) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0237324) q[0];
sx q[0];
rz(-1.5845789) q[0];
sx q[0];
rz(-1.3316863) q[0];
rz(-0.82904077) q[1];
sx q[1];
rz(-1.6571028) q[1];
sx q[1];
rz(2.170457) q[1];
rz(-1.2185417) q[2];
sx q[2];
rz(-0.92377077) q[2];
sx q[2];
rz(0.9818264) q[2];
rz(0.42468023) q[3];
sx q[3];
rz(-1.1241354) q[3];
sx q[3];
rz(-0.507171) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
