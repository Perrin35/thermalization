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
rz(-3.0814085) q[0];
sx q[0];
rz(-1.0589851) q[0];
sx q[0];
rz(-2.1314148) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(-0.29616907) q[1];
sx q[1];
rz(0.30997601) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.298807) q[0];
sx q[0];
rz(-1.1514542) q[0];
sx q[0];
rz(-1.692665) q[0];
rz(3.0344738) q[2];
sx q[2];
rz(-2.9351882) q[2];
sx q[2];
rz(-0.53910461) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1811318) q[1];
sx q[1];
rz(-2.7560661) q[1];
sx q[1];
rz(-1.0268282) q[1];
x q[2];
rz(-2.3100389) q[3];
sx q[3];
rz(-1.420701) q[3];
sx q[3];
rz(-2.7287366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4787204) q[2];
sx q[2];
rz(-0.24820776) q[2];
sx q[2];
rz(-0.09566801) q[2];
rz(-0.11224789) q[3];
sx q[3];
rz(-0.87533689) q[3];
sx q[3];
rz(1.7101425) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2410759) q[0];
sx q[0];
rz(-2.2523585) q[0];
sx q[0];
rz(-0.61479968) q[0];
rz(-2.2531033) q[1];
sx q[1];
rz(-1.5526086) q[1];
sx q[1];
rz(-0.49957553) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7730656) q[0];
sx q[0];
rz(-1.245226) q[0];
sx q[0];
rz(-2.9286309) q[0];
rz(-pi) q[1];
rz(2.8728054) q[2];
sx q[2];
rz(-2.6013881) q[2];
sx q[2];
rz(-2.7160783) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0608637) q[1];
sx q[1];
rz(-1.398096) q[1];
sx q[1];
rz(0.91367803) q[1];
rz(1.0811483) q[3];
sx q[3];
rz(-2.5003365) q[3];
sx q[3];
rz(-1.7288417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2483612) q[2];
sx q[2];
rz(-1.2373368) q[2];
sx q[2];
rz(3.1207747) q[2];
rz(1.9013532) q[3];
sx q[3];
rz(-1.4695243) q[3];
sx q[3];
rz(-1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6388539) q[0];
sx q[0];
rz(-0.26832142) q[0];
sx q[0];
rz(-2.8079206) q[0];
rz(0.73792136) q[1];
sx q[1];
rz(-1.2260022) q[1];
sx q[1];
rz(3.0121682) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0282818) q[0];
sx q[0];
rz(-2.2419562) q[0];
sx q[0];
rz(-2.9379803) q[0];
x q[1];
rz(-0.67141165) q[2];
sx q[2];
rz(-1.08403) q[2];
sx q[2];
rz(-1.1306131) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94679615) q[1];
sx q[1];
rz(-1.692564) q[1];
sx q[1];
rz(0.090784723) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3593216) q[3];
sx q[3];
rz(-1.2201628) q[3];
sx q[3];
rz(0.80814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1239329) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(1.370149) q[2];
rz(2.395199) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(-0.082775041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.133404) q[0];
sx q[0];
rz(-1.9318102) q[0];
sx q[0];
rz(-2.5764537) q[0];
rz(0.42916974) q[1];
sx q[1];
rz(-2.4950835) q[1];
sx q[1];
rz(2.3133004) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9681184) q[0];
sx q[0];
rz(-2.3096363) q[0];
sx q[0];
rz(1.4627187) q[0];
x q[1];
rz(-2.0798551) q[2];
sx q[2];
rz(-2.5974496) q[2];
sx q[2];
rz(0.5296112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8737303) q[1];
sx q[1];
rz(-2.5791188) q[1];
sx q[1];
rz(-2.6270694) q[1];
rz(-pi) q[2];
rz(1.1513814) q[3];
sx q[3];
rz(-2.5852301) q[3];
sx q[3];
rz(-2.0030947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4565178) q[2];
sx q[2];
rz(-2.0802616) q[2];
sx q[2];
rz(1.7100517) q[2];
rz(-2.2020014) q[3];
sx q[3];
rz(-2.6288433) q[3];
sx q[3];
rz(-0.00069869839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13570304) q[0];
sx q[0];
rz(-2.8764184) q[0];
sx q[0];
rz(-1.3294719) q[0];
rz(-1.1520518) q[1];
sx q[1];
rz(-1.0877129) q[1];
sx q[1];
rz(0.00024814127) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0572994) q[0];
sx q[0];
rz(-1.3868185) q[0];
sx q[0];
rz(-1.3217682) q[0];
x q[1];
rz(1.423944) q[2];
sx q[2];
rz(-0.53895437) q[2];
sx q[2];
rz(-2.1951064) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.50284144) q[1];
sx q[1];
rz(-1.0783245) q[1];
sx q[1];
rz(0.3943464) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4100634) q[3];
sx q[3];
rz(-1.6247182) q[3];
sx q[3];
rz(-1.1889072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0033215) q[2];
sx q[2];
rz(-1.891581) q[2];
sx q[2];
rz(0.0058343466) q[2];
rz(0.88998574) q[3];
sx q[3];
rz(-2.4599288) q[3];
sx q[3];
rz(1.1267004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24790813) q[0];
sx q[0];
rz(-1.4434781) q[0];
sx q[0];
rz(1.4877315) q[0];
rz(-3.1112025) q[1];
sx q[1];
rz(-1.1266212) q[1];
sx q[1];
rz(-0.94246513) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2104032) q[0];
sx q[0];
rz(-1.382834) q[0];
sx q[0];
rz(0.24954777) q[0];
x q[1];
rz(-1.312227) q[2];
sx q[2];
rz(-1.8601918) q[2];
sx q[2];
rz(-1.0602151) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16557517) q[1];
sx q[1];
rz(-1.8569131) q[1];
sx q[1];
rz(2.0598434) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0464155) q[3];
sx q[3];
rz(-1.3381357) q[3];
sx q[3];
rz(-1.3568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59948644) q[2];
sx q[2];
rz(-2.3607871) q[2];
sx q[2];
rz(-1.3055118) q[2];
rz(-0.54715884) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(-0.80593306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18043537) q[0];
sx q[0];
rz(-1.0338217) q[0];
sx q[0];
rz(-0.10051522) q[0];
rz(-2.6590977) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(2.2844792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53062886) q[0];
sx q[0];
rz(-2.2016494) q[0];
sx q[0];
rz(-1.1683589) q[0];
rz(-pi) q[1];
rz(-2.4947462) q[2];
sx q[2];
rz(-1.610272) q[2];
sx q[2];
rz(-2.0496862) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.313011) q[1];
sx q[1];
rz(-1.9691113) q[1];
sx q[1];
rz(-1.7969153) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1825652) q[3];
sx q[3];
rz(-2.534158) q[3];
sx q[3];
rz(1.5608112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7414005) q[2];
sx q[2];
rz(-2.1705748) q[2];
sx q[2];
rz(-0.3024438) q[2];
rz(0.97964573) q[3];
sx q[3];
rz(-2.1392348) q[3];
sx q[3];
rz(-1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.361146) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(3.1106023) q[0];
rz(-2.567645) q[1];
sx q[1];
rz(-1.6684883) q[1];
sx q[1];
rz(1.2773638) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.663186) q[0];
sx q[0];
rz(-2.660523) q[0];
sx q[0];
rz(-0.78669725) q[0];
x q[1];
rz(-1.1436815) q[2];
sx q[2];
rz(-0.53887212) q[2];
sx q[2];
rz(0.32921916) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3801703) q[1];
sx q[1];
rz(-1.4974277) q[1];
sx q[1];
rz(1.0835034) q[1];
rz(-3.0798036) q[3];
sx q[3];
rz(-1.8553858) q[3];
sx q[3];
rz(-2.9348843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0370827) q[2];
sx q[2];
rz(-0.13585486) q[2];
sx q[2];
rz(2.6541397) q[2];
rz(-0.69665748) q[3];
sx q[3];
rz(-2.2127377) q[3];
sx q[3];
rz(-3.0776183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697407) q[0];
sx q[0];
rz(-2.6672279) q[0];
sx q[0];
rz(2.790614) q[0];
rz(-2.9674496) q[1];
sx q[1];
rz(-1.6330279) q[1];
sx q[1];
rz(1.7399656) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1202209) q[0];
sx q[0];
rz(-1.8184156) q[0];
sx q[0];
rz(1.6970474) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3910037) q[2];
sx q[2];
rz(-1.2341414) q[2];
sx q[2];
rz(-0.85699334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9822787) q[1];
sx q[1];
rz(-1.3873006) q[1];
sx q[1];
rz(3.112041) q[1];
rz(-0.63302682) q[3];
sx q[3];
rz(-0.84869519) q[3];
sx q[3];
rz(-2.4304469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.031124) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(2.5049211) q[2];
rz(-2.0311671) q[3];
sx q[3];
rz(-1.5444376) q[3];
sx q[3];
rz(2.6231664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3670032) q[0];
sx q[0];
rz(-2.84802) q[0];
sx q[0];
rz(1.0585744) q[0];
rz(3.1029347) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(1.0640594) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2702613) q[0];
sx q[0];
rz(-3.1360021) q[0];
sx q[0];
rz(0.63951512) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27711192) q[2];
sx q[2];
rz(-1.849035) q[2];
sx q[2];
rz(0.81467512) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41760379) q[1];
sx q[1];
rz(-0.072712459) q[1];
sx q[1];
rz(0.59007187) q[1];
rz(-pi) q[2];
rz(-2.6820716) q[3];
sx q[3];
rz(-1.6020613) q[3];
sx q[3];
rz(0.36324901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71558636) q[2];
sx q[2];
rz(-0.49804372) q[2];
sx q[2];
rz(-0.074020298) q[2];
rz(-2.7600539) q[3];
sx q[3];
rz(-1.3149202) q[3];
sx q[3];
rz(0.35416245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.030180177) q[0];
sx q[0];
rz(-1.8288061) q[0];
sx q[0];
rz(-0.70963138) q[0];
rz(-2.5297655) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(0.50449087) q[2];
sx q[2];
rz(-1.2332543) q[2];
sx q[2];
rz(2.6674657) q[2];
rz(2.6388219) q[3];
sx q[3];
rz(-2.3270861) q[3];
sx q[3];
rz(2.6181639) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
