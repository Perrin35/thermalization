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
rz(-2.8895145) q[0];
sx q[0];
rz(-0.57555389) q[0];
sx q[0];
rz(-3.0192896) q[0];
rz(1.1072493) q[1];
sx q[1];
rz(5.3218359) q[1];
sx q[1];
rz(9.3864592) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092889812) q[0];
sx q[0];
rz(-0.99703353) q[0];
sx q[0];
rz(0.4452197) q[0];
x q[1];
rz(-1.6898481) q[2];
sx q[2];
rz(-1.6887293) q[2];
sx q[2];
rz(-2.8182507) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9895175) q[1];
sx q[1];
rz(-0.65836009) q[1];
sx q[1];
rz(-2.1170298) q[1];
rz(-1.5502717) q[3];
sx q[3];
rz(-1.3234659) q[3];
sx q[3];
rz(1.0350682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64500874) q[2];
sx q[2];
rz(-1.714548) q[2];
sx q[2];
rz(-1.4487779) q[2];
rz(2.951176) q[3];
sx q[3];
rz(-1.0551635) q[3];
sx q[3];
rz(-0.14934389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377624) q[0];
sx q[0];
rz(-1.2685403) q[0];
sx q[0];
rz(2.8835836) q[0];
rz(1.5915271) q[1];
sx q[1];
rz(-1.240088) q[1];
sx q[1];
rz(-2.9749427) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7014871) q[0];
sx q[0];
rz(-2.1512554) q[0];
sx q[0];
rz(1.6850182) q[0];
x q[1];
rz(-0.3853674) q[2];
sx q[2];
rz(-1.8358942) q[2];
sx q[2];
rz(-2.6814658) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.38063654) q[1];
sx q[1];
rz(-1.6957307) q[1];
sx q[1];
rz(2.4255468) q[1];
x q[2];
rz(0.54661481) q[3];
sx q[3];
rz(-2.0677462) q[3];
sx q[3];
rz(0.95595804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0717281) q[2];
sx q[2];
rz(-1.122415) q[2];
sx q[2];
rz(1.2321164) q[2];
rz(-0.92464906) q[3];
sx q[3];
rz(-2.2053714) q[3];
sx q[3];
rz(1.1054976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.449618) q[0];
sx q[0];
rz(-1.6716577) q[0];
sx q[0];
rz(-0.0019419226) q[0];
rz(-3.107403) q[1];
sx q[1];
rz(-1.9296153) q[1];
sx q[1];
rz(-1.5431822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6337834) q[0];
sx q[0];
rz(-0.52716053) q[0];
sx q[0];
rz(2.6109004) q[0];
rz(1.9667186) q[2];
sx q[2];
rz(-2.2668419) q[2];
sx q[2];
rz(-0.093122236) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0405827) q[1];
sx q[1];
rz(-1.9242745) q[1];
sx q[1];
rz(1.6124875) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8887599) q[3];
sx q[3];
rz(-2.468716) q[3];
sx q[3];
rz(1.8102136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2446642) q[2];
sx q[2];
rz(-2.7074773) q[2];
sx q[2];
rz(-1.0373235) q[2];
rz(2.0364929) q[3];
sx q[3];
rz(-1.5367855) q[3];
sx q[3];
rz(1.1387811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26745519) q[0];
sx q[0];
rz(-0.48478165) q[0];
sx q[0];
rz(-1.4111891) q[0];
rz(0.63938582) q[1];
sx q[1];
rz(-1.4566028) q[1];
sx q[1];
rz(-0.04714084) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1829202) q[0];
sx q[0];
rz(-2.483568) q[0];
sx q[0];
rz(0.47112314) q[0];
x q[1];
rz(1.3762952) q[2];
sx q[2];
rz(-1.10656) q[2];
sx q[2];
rz(0.08438202) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2915708) q[1];
sx q[1];
rz(-2.154989) q[1];
sx q[1];
rz(2.7136346) q[1];
rz(-1.9023667) q[3];
sx q[3];
rz(-0.75115381) q[3];
sx q[3];
rz(-2.0560665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3186657) q[2];
sx q[2];
rz(-1.1456127) q[2];
sx q[2];
rz(1.9226496) q[2];
rz(0.31560358) q[3];
sx q[3];
rz(-0.054840755) q[3];
sx q[3];
rz(2.992673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(1.3217992) q[0];
sx q[0];
rz(-2.5892374) q[0];
sx q[0];
rz(-2.7929982) q[0];
rz(1.2527342) q[1];
sx q[1];
rz(-1.9786973) q[1];
sx q[1];
rz(-1.3016275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0785261) q[0];
sx q[0];
rz(-1.6498897) q[0];
sx q[0];
rz(2.3344759) q[0];
rz(-1.9867861) q[2];
sx q[2];
rz(-2.092053) q[2];
sx q[2];
rz(-0.74300569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0730473) q[1];
sx q[1];
rz(-1.7918192) q[1];
sx q[1];
rz(2.9566492) q[1];
rz(-pi) q[2];
rz(2.1215277) q[3];
sx q[3];
rz(-1.8740843) q[3];
sx q[3];
rz(2.1316949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18315135) q[2];
sx q[2];
rz(-2.8330467) q[2];
sx q[2];
rz(1.6661673) q[2];
rz(-0.685855) q[3];
sx q[3];
rz(-2.1619022) q[3];
sx q[3];
rz(0.53387749) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1118065) q[0];
sx q[0];
rz(-2.989558) q[0];
sx q[0];
rz(1.6492122) q[0];
rz(-0.97914186) q[1];
sx q[1];
rz(-1.5359595) q[1];
sx q[1];
rz(-1.2243366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92517744) q[0];
sx q[0];
rz(-1.6013711) q[0];
sx q[0];
rz(0.19801099) q[0];
x q[1];
rz(0.54927214) q[2];
sx q[2];
rz(-1.7584929) q[2];
sx q[2];
rz(-0.21085462) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8414014) q[1];
sx q[1];
rz(-2.527664) q[1];
sx q[1];
rz(1.5924686) q[1];
rz(-pi) q[2];
rz(2.200279) q[3];
sx q[3];
rz(-1.8555897) q[3];
sx q[3];
rz(-0.79508699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7225723) q[2];
sx q[2];
rz(-1.6442862) q[2];
sx q[2];
rz(2.2255619) q[2];
rz(1.7570868) q[3];
sx q[3];
rz(-1.8839096) q[3];
sx q[3];
rz(-0.70703435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(3.1322121) q[0];
sx q[0];
rz(-3.0896602) q[0];
sx q[0];
rz(2.1873271) q[0];
rz(-0.43680278) q[1];
sx q[1];
rz(-1.5780459) q[1];
sx q[1];
rz(-0.11016914) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4797929) q[0];
sx q[0];
rz(-0.40783007) q[0];
sx q[0];
rz(-2.5026425) q[0];
rz(-pi) q[1];
rz(-1.0344347) q[2];
sx q[2];
rz(-2.1948994) q[2];
sx q[2];
rz(1.0076866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84444204) q[1];
sx q[1];
rz(-1.5664829) q[1];
sx q[1];
rz(1.9164852) q[1];
rz(-pi) q[2];
rz(2.0721335) q[3];
sx q[3];
rz(-1.7946464) q[3];
sx q[3];
rz(1.7373794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43180141) q[2];
sx q[2];
rz(-3.1352037) q[2];
sx q[2];
rz(-2.1978281) q[2];
rz(0.56378311) q[3];
sx q[3];
rz(-1.0958593) q[3];
sx q[3];
rz(1.110466) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875882) q[0];
sx q[0];
rz(-2.8546951) q[0];
sx q[0];
rz(2.7960844) q[0];
rz(2.0461138) q[1];
sx q[1];
rz(-1.6637207) q[1];
sx q[1];
rz(-1.5798205) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51673543) q[0];
sx q[0];
rz(-0.78291946) q[0];
sx q[0];
rz(1.5262414) q[0];
x q[1];
rz(1.6510294) q[2];
sx q[2];
rz(-0.65872619) q[2];
sx q[2];
rz(1.4874489) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53198481) q[1];
sx q[1];
rz(-0.72276211) q[1];
sx q[1];
rz(-2.1400129) q[1];
rz(1.5388266) q[3];
sx q[3];
rz(-1.073146) q[3];
sx q[3];
rz(-2.9835193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3788508) q[2];
sx q[2];
rz(-0.26143917) q[2];
sx q[2];
rz(1.4307107) q[2];
rz(-2.6070969) q[3];
sx q[3];
rz(-1.4662687) q[3];
sx q[3];
rz(-0.9526332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6119824) q[0];
sx q[0];
rz(-3.072325) q[0];
sx q[0];
rz(-1.8310504) q[0];
rz(0.88690859) q[1];
sx q[1];
rz(-2.3014018) q[1];
sx q[1];
rz(-1.8720253) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16989141) q[0];
sx q[0];
rz(-1.9215688) q[0];
sx q[0];
rz(3.0049075) q[0];
x q[1];
rz(-2.9863556) q[2];
sx q[2];
rz(-1.9433776) q[2];
sx q[2];
rz(-2.9205204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9880105) q[1];
sx q[1];
rz(-1.8021823) q[1];
sx q[1];
rz(-1.2910976) q[1];
rz(-pi) q[2];
rz(2.9501602) q[3];
sx q[3];
rz(-0.20729724) q[3];
sx q[3];
rz(1.8199004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5294007) q[2];
sx q[2];
rz(-1.3357013) q[2];
sx q[2];
rz(0.23294918) q[2];
rz(-1.285078) q[3];
sx q[3];
rz(-0.67684567) q[3];
sx q[3];
rz(2.7555833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1583629) q[0];
sx q[0];
rz(-0.60929275) q[0];
sx q[0];
rz(0.097271517) q[0];
rz(-2.3376047) q[1];
sx q[1];
rz(-2.4318047) q[1];
sx q[1];
rz(2.9248617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0618503) q[0];
sx q[0];
rz(-1.5155795) q[0];
sx q[0];
rz(-1.7415857) q[0];
rz(-pi) q[1];
rz(-1.5243657) q[2];
sx q[2];
rz(-1.7123607) q[2];
sx q[2];
rz(-2.9449449) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35516675) q[1];
sx q[1];
rz(-0.92235074) q[1];
sx q[1];
rz(-0.72079682) q[1];
rz(-pi) q[2];
rz(-0.43301591) q[3];
sx q[3];
rz(-0.70954126) q[3];
sx q[3];
rz(-2.7482928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8486166) q[2];
sx q[2];
rz(-2.8605707) q[2];
sx q[2];
rz(2.6463553) q[2];
rz(-1.5589335) q[3];
sx q[3];
rz(-2.0571183) q[3];
sx q[3];
rz(-2.9787279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0611298) q[0];
sx q[0];
rz(-1.6292138) q[0];
sx q[0];
rz(-0.52393352) q[0];
rz(-1.0864661) q[1];
sx q[1];
rz(-2.6230984) q[1];
sx q[1];
rz(-2.7883504) q[1];
rz(-1.5389961) q[2];
sx q[2];
rz(-1.563579) q[2];
sx q[2];
rz(-2.210571) q[2];
rz(-1.466352) q[3];
sx q[3];
rz(-0.71487311) q[3];
sx q[3];
rz(-2.5345595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
