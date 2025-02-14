OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7393957) q[0];
sx q[0];
rz(-2.2576809) q[0];
sx q[0];
rz(-0.85564268) q[0];
rz(-0.17172509) q[1];
sx q[1];
rz(-0.11556927) q[1];
sx q[1];
rz(2.5791383) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2320975) q[0];
sx q[0];
rz(-1.2065071) q[0];
sx q[0];
rz(-3.0845736) q[0];
x q[1];
rz(2.0777656) q[2];
sx q[2];
rz(-1.7071144) q[2];
sx q[2];
rz(2.9796114) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3823279) q[1];
sx q[1];
rz(-2.538343) q[1];
sx q[1];
rz(-0.77253491) q[1];
rz(-pi) q[2];
rz(-1.9770369) q[3];
sx q[3];
rz(-1.2066168) q[3];
sx q[3];
rz(1.6842357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90613753) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(2.3423024) q[2];
rz(2.6702787) q[3];
sx q[3];
rz(-0.91336942) q[3];
sx q[3];
rz(-0.93588626) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039463194) q[0];
sx q[0];
rz(-1.3251745) q[0];
sx q[0];
rz(-1.7934196) q[0];
rz(-1.7680291) q[1];
sx q[1];
rz(-1.9919688) q[1];
sx q[1];
rz(-1.9893533) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0420899) q[0];
sx q[0];
rz(-1.2782017) q[0];
sx q[0];
rz(-0.83935229) q[0];
rz(-pi) q[1];
rz(0.67518465) q[2];
sx q[2];
rz(-1.0567046) q[2];
sx q[2];
rz(-2.6509283) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3481231) q[1];
sx q[1];
rz(-1.4033699) q[1];
sx q[1];
rz(0.69116418) q[1];
x q[2];
rz(2.3186734) q[3];
sx q[3];
rz(-2.4429818) q[3];
sx q[3];
rz(0.60958344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79933244) q[2];
sx q[2];
rz(-2.7228184) q[2];
sx q[2];
rz(0.93117923) q[2];
rz(3.0350507) q[3];
sx q[3];
rz(-1.1139161) q[3];
sx q[3];
rz(-0.60025269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0837285) q[0];
sx q[0];
rz(-1.3902384) q[0];
sx q[0];
rz(-2.6237543) q[0];
rz(-0.88874108) q[1];
sx q[1];
rz(-2.4338212) q[1];
sx q[1];
rz(0.4471561) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6155375) q[0];
sx q[0];
rz(-1.3096022) q[0];
sx q[0];
rz(1.5047856) q[0];
rz(-2.4141623) q[2];
sx q[2];
rz(-1.1377678) q[2];
sx q[2];
rz(-0.43011452) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39147705) q[1];
sx q[1];
rz(-2.1575054) q[1];
sx q[1];
rz(0.16280414) q[1];
x q[2];
rz(-2.6833862) q[3];
sx q[3];
rz(-0.57469207) q[3];
sx q[3];
rz(-1.1632869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37457028) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(-0.53263295) q[2];
rz(2.8940708) q[3];
sx q[3];
rz(-2.4041924) q[3];
sx q[3];
rz(1.0386764) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6240876) q[0];
sx q[0];
rz(-0.085260304) q[0];
sx q[0];
rz(3.0349773) q[0];
rz(2.7952349) q[1];
sx q[1];
rz(-2.296591) q[1];
sx q[1];
rz(-1.9812298) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896987) q[0];
sx q[0];
rz(-1.3972613) q[0];
sx q[0];
rz(0.24994295) q[0];
rz(-pi) q[1];
rz(0.8125272) q[2];
sx q[2];
rz(-0.90355325) q[2];
sx q[2];
rz(3.0548801) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9306158) q[1];
sx q[1];
rz(-1.2158356) q[1];
sx q[1];
rz(-0.32101722) q[1];
rz(-pi) q[2];
rz(-2.2665146) q[3];
sx q[3];
rz(-2.4949673) q[3];
sx q[3];
rz(2.2602606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0059263) q[2];
sx q[2];
rz(-2.8432507) q[2];
sx q[2];
rz(-3.0653817) q[2];
rz(-2.5557319) q[3];
sx q[3];
rz(-1.9867089) q[3];
sx q[3];
rz(1.7990254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0214486) q[0];
sx q[0];
rz(-1.3055389) q[0];
sx q[0];
rz(2.7914877) q[0];
rz(0.95589751) q[1];
sx q[1];
rz(-1.8616385) q[1];
sx q[1];
rz(1.2299445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6087225) q[0];
sx q[0];
rz(-1.9367095) q[0];
sx q[0];
rz(-2.8844112) q[0];
x q[1];
rz(1.7361264) q[2];
sx q[2];
rz(-1.606719) q[2];
sx q[2];
rz(1.3816116) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38795162) q[1];
sx q[1];
rz(-1.3820096) q[1];
sx q[1];
rz(-1.2354047) q[1];
rz(-pi) q[2];
rz(-0.1996207) q[3];
sx q[3];
rz(-0.99310447) q[3];
sx q[3];
rz(-2.9132089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.89796394) q[2];
sx q[2];
rz(-0.25883365) q[2];
sx q[2];
rz(1.6492856) q[2];
rz(-2.7367075) q[3];
sx q[3];
rz(-2.3758774) q[3];
sx q[3];
rz(-0.83612061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355857) q[0];
sx q[0];
rz(-1.0772935) q[0];
sx q[0];
rz(0.58240044) q[0];
rz(-0.68663418) q[1];
sx q[1];
rz(-1.735894) q[1];
sx q[1];
rz(-1.2581717) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63623896) q[0];
sx q[0];
rz(-2.6864144) q[0];
sx q[0];
rz(-1.8591465) q[0];
rz(-pi) q[1];
rz(-2.7863726) q[2];
sx q[2];
rz(-2.0783011) q[2];
sx q[2];
rz(1.2025646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4542089) q[1];
sx q[1];
rz(-1.6404248) q[1];
sx q[1];
rz(1.5147665) q[1];
x q[2];
rz(0.97719394) q[3];
sx q[3];
rz(-1.753429) q[3];
sx q[3];
rz(0.066150811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76576343) q[2];
sx q[2];
rz(-1.9932237) q[2];
sx q[2];
rz(2.1235535) q[2];
rz(-2.8954519) q[3];
sx q[3];
rz(-1.3956416) q[3];
sx q[3];
rz(-1.6233981) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70596424) q[0];
sx q[0];
rz(-2.3376597) q[0];
sx q[0];
rz(-0.58018082) q[0];
rz(0.14353453) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(2.9339583) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26941368) q[0];
sx q[0];
rz(-2.2996443) q[0];
sx q[0];
rz(0.4704041) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7933664) q[2];
sx q[2];
rz(-2.1219606) q[2];
sx q[2];
rz(-1.4081692) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8676198) q[1];
sx q[1];
rz(-1.8302396) q[1];
sx q[1];
rz(-2.2745489) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24702279) q[3];
sx q[3];
rz(-2.7678856) q[3];
sx q[3];
rz(1.6173687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84270728) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(-2.5353954) q[2];
rz(-0.68743622) q[3];
sx q[3];
rz(-2.0076553) q[3];
sx q[3];
rz(-2.4351951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6523478) q[0];
sx q[0];
rz(-3.1135961) q[0];
sx q[0];
rz(-2.0835173) q[0];
rz(-0.10969133) q[1];
sx q[1];
rz(-2.0249764) q[1];
sx q[1];
rz(1.6995957) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88272754) q[0];
sx q[0];
rz(-1.8804714) q[0];
sx q[0];
rz(0.41559269) q[0];
rz(2.6410854) q[2];
sx q[2];
rz(-1.3277413) q[2];
sx q[2];
rz(-2.1432723) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.31046384) q[1];
sx q[1];
rz(-2.052124) q[1];
sx q[1];
rz(0.83408611) q[1];
rz(1.2992925) q[3];
sx q[3];
rz(-1.7978923) q[3];
sx q[3];
rz(2.5271067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66120061) q[2];
sx q[2];
rz(-0.19442393) q[2];
sx q[2];
rz(1.2083758) q[2];
rz(-2.4783573) q[3];
sx q[3];
rz(-1.6845208) q[3];
sx q[3];
rz(1.2164345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97063589) q[0];
sx q[0];
rz(-2.1747776) q[0];
sx q[0];
rz(-3.048625) q[0];
rz(1.8611106) q[1];
sx q[1];
rz(-0.72851506) q[1];
sx q[1];
rz(3.0063937) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2743535) q[0];
sx q[0];
rz(-1.6238191) q[0];
sx q[0];
rz(1.5252602) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5298631) q[2];
sx q[2];
rz(-0.29680291) q[2];
sx q[2];
rz(-0.57020818) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1982806) q[1];
sx q[1];
rz(-1.3168646) q[1];
sx q[1];
rz(-1.5487681) q[1];
x q[2];
rz(-2.0207094) q[3];
sx q[3];
rz(-2.9723047) q[3];
sx q[3];
rz(1.560488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4325503) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(2.3700628) q[2];
rz(-2.6500474) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(-0.5947203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5563357) q[0];
sx q[0];
rz(-2.519643) q[0];
sx q[0];
rz(-1.1248032) q[0];
rz(-1.8428165) q[1];
sx q[1];
rz(-0.61538428) q[1];
sx q[1];
rz(-0.76464701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1021834) q[0];
sx q[0];
rz(-1.4847857) q[0];
sx q[0];
rz(3.0360704) q[0];
rz(-0.89057335) q[2];
sx q[2];
rz(-1.9720417) q[2];
sx q[2];
rz(-0.97637343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.56988) q[1];
sx q[1];
rz(-0.24001828) q[1];
sx q[1];
rz(-2.1313138) q[1];
rz(-0.57101698) q[3];
sx q[3];
rz(-2.5072376) q[3];
sx q[3];
rz(-0.77658021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7654968) q[2];
sx q[2];
rz(-1.2906047) q[2];
sx q[2];
rz(0.20720227) q[2];
rz(-0.97992212) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(-0.19237147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1205263) q[0];
sx q[0];
rz(-1.6854032) q[0];
sx q[0];
rz(2.2717463) q[0];
rz(-0.5008685) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(-0.7122559) q[2];
sx q[2];
rz(-2.9604572) q[2];
sx q[2];
rz(-1.0618718) q[2];
rz(1.408314) q[3];
sx q[3];
rz(-1.3357031) q[3];
sx q[3];
rz(1.6817844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
