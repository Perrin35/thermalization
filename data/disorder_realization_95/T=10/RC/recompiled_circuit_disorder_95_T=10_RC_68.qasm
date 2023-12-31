OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(-2.2440417) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25047725) q[0];
sx q[0];
rz(-2.5076137) q[0];
sx q[0];
rz(-1.4527713) q[0];
x q[1];
rz(0.35383309) q[2];
sx q[2];
rz(-0.96828038) q[2];
sx q[2];
rz(1.8941855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.861046) q[1];
sx q[1];
rz(-1.4006873) q[1];
sx q[1];
rz(-2.4982846) q[1];
rz(1.3863871) q[3];
sx q[3];
rz(-2.3204436) q[3];
sx q[3];
rz(0.26339312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.477318) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(2.409639) q[2];
rz(2.1814363) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(-1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9663548) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(2.2251341) q[0];
rz(2.6610999) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(-0.8786456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7322313) q[0];
sx q[0];
rz(-1.9525813) q[0];
sx q[0];
rz(0.9071) q[0];
x q[1];
rz(0.61848817) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(2.6478812) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6138184) q[1];
sx q[1];
rz(-1.4973745) q[1];
sx q[1];
rz(0.41927494) q[1];
x q[2];
rz(3.0733172) q[3];
sx q[3];
rz(-1.8797415) q[3];
sx q[3];
rz(-2.6451049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(2.4943165) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(-0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(0.24599427) q[0];
rz(1.4100769) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(-2.0203967) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66378731) q[0];
sx q[0];
rz(-0.87513798) q[0];
sx q[0];
rz(-1.0870766) q[0];
rz(-pi) q[1];
rz(-1.7240702) q[2];
sx q[2];
rz(-0.3970662) q[2];
sx q[2];
rz(0.085445554) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9417291) q[1];
sx q[1];
rz(-0.45048303) q[1];
sx q[1];
rz(-0.73168879) q[1];
x q[2];
rz(-0.49384533) q[3];
sx q[3];
rz(-1.4673127) q[3];
sx q[3];
rz(-2.8007357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(-2.4915063) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(-0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6275416) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(-0.78805584) q[0];
rz(1.0568985) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(3.025211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560293) q[0];
sx q[0];
rz(-0.93755994) q[0];
sx q[0];
rz(-1.602519) q[0];
rz(2.0747244) q[2];
sx q[2];
rz(-1.2328086) q[2];
sx q[2];
rz(2.4525814) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0105646) q[1];
sx q[1];
rz(-0.66337913) q[1];
sx q[1];
rz(2.4478711) q[1];
rz(2.651752) q[3];
sx q[3];
rz(-1.6162989) q[3];
sx q[3];
rz(3.0385466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6115761) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(-0.3616412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5761121) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(0.53260032) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(-1.1486357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2411597) q[0];
sx q[0];
rz(-2.1923089) q[0];
sx q[0];
rz(-0.09089367) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45102851) q[2];
sx q[2];
rz(-2.1941059) q[2];
sx q[2];
rz(0.23652467) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27302882) q[1];
sx q[1];
rz(-1.7595353) q[1];
sx q[1];
rz(-2.7815458) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7025181) q[3];
sx q[3];
rz(-2.3665161) q[3];
sx q[3];
rz(2.4174158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1308412) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(2.1095236) q[2];
rz(-0.71470913) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(-3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59153581) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(-0.39598879) q[0];
rz(1.4453325) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(-0.2063624) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26158953) q[0];
sx q[0];
rz(-2.3149009) q[0];
sx q[0];
rz(-1.4522626) q[0];
rz(-pi) q[1];
rz(-2.4762819) q[2];
sx q[2];
rz(-0.99669257) q[2];
sx q[2];
rz(2.2523508) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.98112647) q[1];
sx q[1];
rz(-1.4741352) q[1];
sx q[1];
rz(-1.4658982) q[1];
rz(2.9456375) q[3];
sx q[3];
rz(-1.8926419) q[3];
sx q[3];
rz(-2.4623507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.53987327) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(-2.8708141) q[2];
rz(2.3932636) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(-1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.2449743) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(-0.57624972) q[0];
rz(0.17310625) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(-1.9304088) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1887218) q[0];
sx q[0];
rz(-1.0591918) q[0];
sx q[0];
rz(2.1923724) q[0];
rz(-1.4326101) q[2];
sx q[2];
rz(-2.0565363) q[2];
sx q[2];
rz(-1.1952343) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3019575) q[1];
sx q[1];
rz(-2.9710899) q[1];
sx q[1];
rz(-1.1595999) q[1];
rz(-pi) q[2];
rz(0.96774775) q[3];
sx q[3];
rz(-1.3021886) q[3];
sx q[3];
rz(2.8318162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(-2.426614) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(-1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054758469) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(-2.1210282) q[0];
rz(-0.15469805) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(-1.221009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2894665) q[0];
sx q[0];
rz(-0.8284491) q[0];
sx q[0];
rz(2.7702987) q[0];
x q[1];
rz(1.4079354) q[2];
sx q[2];
rz(-1.7454299) q[2];
sx q[2];
rz(-0.83144951) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78302661) q[1];
sx q[1];
rz(-1.6811803) q[1];
sx q[1];
rz(3.0049938) q[1];
rz(-2.8052748) q[3];
sx q[3];
rz(-0.85018966) q[3];
sx q[3];
rz(1.1433126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(-0.60338902) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.110638) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(0.13701339) q[0];
rz(-2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(-3.0158214) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.754697) q[0];
sx q[0];
rz(-1.7594975) q[0];
sx q[0];
rz(2.9625091) q[0];
rz(-pi) q[1];
x q[1];
rz(2.618082) q[2];
sx q[2];
rz(-2.7611809) q[2];
sx q[2];
rz(-1.3695804) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3667664) q[1];
sx q[1];
rz(-0.80087304) q[1];
sx q[1];
rz(1.1714539) q[1];
rz(2.845876) q[3];
sx q[3];
rz(-1.3774646) q[3];
sx q[3];
rz(0.28702345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13828364) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(-2.4323145) q[2];
rz(2.5214031) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(1.2003157) q[0];
rz(2.8740846) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(2.0013924) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2518371) q[0];
sx q[0];
rz(-1.4714843) q[0];
sx q[0];
rz(-1.3072144) q[0];
rz(-pi) q[1];
rz(2.0613725) q[2];
sx q[2];
rz(-2.6274101) q[2];
sx q[2];
rz(-2.0928004) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4789341) q[1];
sx q[1];
rz(-1.3824029) q[1];
sx q[1];
rz(1.4807832) q[1];
x q[2];
rz(1.8885918) q[3];
sx q[3];
rz(-2.3386049) q[3];
sx q[3];
rz(-3.0190937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7662979) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(2.2429788) q[2];
rz(3.1344154) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2469149) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-0.49171369) q[2];
sx q[2];
rz(-1.1469054) q[2];
sx q[2];
rz(2.7804874) q[2];
rz(-0.58017147) q[3];
sx q[3];
rz(-2.4709354) q[3];
sx q[3];
rz(1.2448268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
