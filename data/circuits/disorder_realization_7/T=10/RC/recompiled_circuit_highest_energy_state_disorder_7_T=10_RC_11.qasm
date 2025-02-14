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
rz(1.3697019) q[0];
sx q[0];
rz(5.852795) q[0];
sx q[0];
rz(9.2601321) q[0];
rz(-1.6904263) q[1];
sx q[1];
rz(-0.37835205) q[1];
sx q[1];
rz(0.70986706) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71043308) q[0];
sx q[0];
rz(-1.834615) q[0];
sx q[0];
rz(0.13508787) q[0];
x q[1];
rz(1.1426635) q[2];
sx q[2];
rz(-1.1096508) q[2];
sx q[2];
rz(-2.3335339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.134608) q[1];
sx q[1];
rz(-0.35672327) q[1];
sx q[1];
rz(-1.7118042) q[1];
rz(0.48810256) q[3];
sx q[3];
rz(-0.26345601) q[3];
sx q[3];
rz(1.7149705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6115438) q[2];
sx q[2];
rz(-0.350746) q[2];
sx q[2];
rz(-1.9225527) q[2];
rz(2.6576095) q[3];
sx q[3];
rz(-2.2958906) q[3];
sx q[3];
rz(-1.3673966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4265863) q[0];
sx q[0];
rz(-0.33892092) q[0];
sx q[0];
rz(-1.4434848) q[0];
rz(1.8679856) q[1];
sx q[1];
rz(-1.0304291) q[1];
sx q[1];
rz(0.66676203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50567195) q[0];
sx q[0];
rz(-1.2742325) q[0];
sx q[0];
rz(1.3937598) q[0];
rz(-1.6826434) q[2];
sx q[2];
rz(-1.6665272) q[2];
sx q[2];
rz(2.7494631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1135865) q[1];
sx q[1];
rz(-2.0087112) q[1];
sx q[1];
rz(2.6922845) q[1];
rz(-pi) q[2];
rz(-1.1902892) q[3];
sx q[3];
rz(-0.78808053) q[3];
sx q[3];
rz(-0.57315592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.71191177) q[2];
sx q[2];
rz(-1.9671755) q[2];
sx q[2];
rz(-2.5453117) q[2];
rz(2.1176254) q[3];
sx q[3];
rz(-1.7496611) q[3];
sx q[3];
rz(1.1202687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.4100818) q[0];
sx q[0];
rz(-2.5937268) q[0];
sx q[0];
rz(-1.4672853) q[0];
rz(0.038854988) q[1];
sx q[1];
rz(-1.7170186) q[1];
sx q[1];
rz(2.255596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17490921) q[0];
sx q[0];
rz(-1.0499448) q[0];
sx q[0];
rz(2.0684) q[0];
x q[1];
rz(-1.5019234) q[2];
sx q[2];
rz(-2.4703272) q[2];
sx q[2];
rz(1.6948989) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6941144) q[1];
sx q[1];
rz(-1.9285818) q[1];
sx q[1];
rz(-1.612979) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4913849) q[3];
sx q[3];
rz(-1.5089499) q[3];
sx q[3];
rz(2.0275778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2140865) q[2];
sx q[2];
rz(-1.4518041) q[2];
sx q[2];
rz(-0.56373325) q[2];
rz(-2.2659414) q[3];
sx q[3];
rz(-2.0893658) q[3];
sx q[3];
rz(-2.0503069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5617274) q[0];
sx q[0];
rz(-0.848333) q[0];
sx q[0];
rz(-1.0680098) q[0];
rz(2.7470398) q[1];
sx q[1];
rz(-0.5296455) q[1];
sx q[1];
rz(-1.4716757) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36763299) q[0];
sx q[0];
rz(-1.9955705) q[0];
sx q[0];
rz(-1.9494809) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9392233) q[2];
sx q[2];
rz(-0.30469182) q[2];
sx q[2];
rz(-0.42542377) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.79395548) q[1];
sx q[1];
rz(-0.82073602) q[1];
sx q[1];
rz(2.1669037) q[1];
rz(1.6789465) q[3];
sx q[3];
rz(-1.1740854) q[3];
sx q[3];
rz(2.0559529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7431405) q[2];
sx q[2];
rz(-1.3116216) q[2];
sx q[2];
rz(3.0121646) q[2];
rz(-2.5982924) q[3];
sx q[3];
rz(-1.0933417) q[3];
sx q[3];
rz(1.3885385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.038079) q[0];
sx q[0];
rz(-2.1383801) q[0];
sx q[0];
rz(2.6858618) q[0];
rz(-2.575846) q[1];
sx q[1];
rz(-0.63876286) q[1];
sx q[1];
rz(2.9260213) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087061398) q[0];
sx q[0];
rz(-2.2534459) q[0];
sx q[0];
rz(0.21933098) q[0];
x q[1];
rz(-0.59026123) q[2];
sx q[2];
rz(-1.0059085) q[2];
sx q[2];
rz(1.5399726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4734772) q[1];
sx q[1];
rz(-1.501742) q[1];
sx q[1];
rz(1.3051239) q[1];
rz(-pi) q[2];
rz(-1.2108675) q[3];
sx q[3];
rz(-2.8382819) q[3];
sx q[3];
rz(3.123542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.7171219) q[2];
sx q[2];
rz(-2.8346546) q[2];
sx q[2];
rz(-1.4324987) q[2];
rz(1.1411544) q[3];
sx q[3];
rz(-1.9211831) q[3];
sx q[3];
rz(-0.46877638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6345217) q[0];
sx q[0];
rz(-0.85352007) q[0];
sx q[0];
rz(2.5216907) q[0];
rz(1.9339405) q[1];
sx q[1];
rz(-0.67283216) q[1];
sx q[1];
rz(-2.8501453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0196592) q[0];
sx q[0];
rz(-1.013151) q[0];
sx q[0];
rz(-2.9479821) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7038149) q[2];
sx q[2];
rz(-1.8916733) q[2];
sx q[2];
rz(1.4747628) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6537854) q[1];
sx q[1];
rz(-1.8099306) q[1];
sx q[1];
rz(2.5433648) q[1];
rz(-pi) q[2];
rz(3.1278947) q[3];
sx q[3];
rz(-0.65675101) q[3];
sx q[3];
rz(-0.17431549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7181299) q[2];
sx q[2];
rz(-0.66897696) q[2];
sx q[2];
rz(-0.97406975) q[2];
rz(-1.918321) q[3];
sx q[3];
rz(-1.4281102) q[3];
sx q[3];
rz(1.0380925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85422) q[0];
sx q[0];
rz(-1.9725476) q[0];
sx q[0];
rz(-2.8337692) q[0];
rz(0.27990118) q[1];
sx q[1];
rz(-2.8925536) q[1];
sx q[1];
rz(1.261796) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1824337) q[0];
sx q[0];
rz(-1.4366581) q[0];
sx q[0];
rz(3.1073551) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7136015) q[2];
sx q[2];
rz(-2.3624785) q[2];
sx q[2];
rz(0.20881685) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7397464) q[1];
sx q[1];
rz(-1.8658499) q[1];
sx q[1];
rz(-3.1298545) q[1];
rz(-pi) q[2];
rz(3.0715176) q[3];
sx q[3];
rz(-1.106602) q[3];
sx q[3];
rz(-1.7572559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6383003) q[2];
sx q[2];
rz(-1.952848) q[2];
sx q[2];
rz(2.3390181) q[2];
rz(-1.5447626) q[3];
sx q[3];
rz(-2.7464726) q[3];
sx q[3];
rz(1.6667268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.06642) q[0];
sx q[0];
rz(-1.5479227) q[0];
sx q[0];
rz(-2.5286034) q[0];
rz(-1.0523484) q[1];
sx q[1];
rz(-2.3520825) q[1];
sx q[1];
rz(1.8863511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8156453) q[0];
sx q[0];
rz(-1.1425352) q[0];
sx q[0];
rz(-2.0518579) q[0];
rz(2.1034715) q[2];
sx q[2];
rz(-1.7591068) q[2];
sx q[2];
rz(-1.7323158) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8415642) q[1];
sx q[1];
rz(-2.2662913) q[1];
sx q[1];
rz(2.5072875) q[1];
rz(-pi) q[2];
rz(-0.46801059) q[3];
sx q[3];
rz(-2.8492894) q[3];
sx q[3];
rz(1.7287031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.5959979) q[2];
sx q[2];
rz(-0.91384807) q[2];
sx q[2];
rz(2.9929898) q[2];
rz(2.4871067) q[3];
sx q[3];
rz(-1.196967) q[3];
sx q[3];
rz(1.7053846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8504836) q[0];
sx q[0];
rz(-1.5330667) q[0];
sx q[0];
rz(3.1403551) q[0];
rz(-1.2507863) q[1];
sx q[1];
rz(-2.2092399) q[1];
sx q[1];
rz(2.1900182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5244999) q[0];
sx q[0];
rz(-1.6902335) q[0];
sx q[0];
rz(-2.9725084) q[0];
x q[1];
rz(-2.1733093) q[2];
sx q[2];
rz(-0.70455019) q[2];
sx q[2];
rz(2.8271443) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39115649) q[1];
sx q[1];
rz(-1.6200778) q[1];
sx q[1];
rz(1.6154013) q[1];
rz(-2.9195027) q[3];
sx q[3];
rz(-2.7324893) q[3];
sx q[3];
rz(-2.3698185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9963659) q[2];
sx q[2];
rz(-1.9731015) q[2];
sx q[2];
rz(2.1053947) q[2];
rz(0.869831) q[3];
sx q[3];
rz(-1.4709604) q[3];
sx q[3];
rz(-0.74365348) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1096126) q[0];
sx q[0];
rz(-2.5514422) q[0];
sx q[0];
rz(1.1085283) q[0];
rz(-2.176586) q[1];
sx q[1];
rz(-1.3771649) q[1];
sx q[1];
rz(0.56934294) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78380221) q[0];
sx q[0];
rz(-2.0817502) q[0];
sx q[0];
rz(1.3207256) q[0];
rz(1.1064208) q[2];
sx q[2];
rz(-1.9074252) q[2];
sx q[2];
rz(-2.6576633) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9534193) q[1];
sx q[1];
rz(-1.0480289) q[1];
sx q[1];
rz(3.1320577) q[1];
x q[2];
rz(0.7580076) q[3];
sx q[3];
rz(-1.692309) q[3];
sx q[3];
rz(1.9415159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0137332) q[2];
sx q[2];
rz(-2.6467085) q[2];
sx q[2];
rz(-3.0950756) q[2];
rz(-0.30139309) q[3];
sx q[3];
rz(-1.9207585) q[3];
sx q[3];
rz(-0.2218328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2409869) q[0];
sx q[0];
rz(-1.931668) q[0];
sx q[0];
rz(1.3264309) q[0];
rz(1.9164512) q[1];
sx q[1];
rz(-2.6381208) q[1];
sx q[1];
rz(-2.0874964) q[1];
rz(-1.7635119) q[2];
sx q[2];
rz(-0.96619923) q[2];
sx q[2];
rz(-0.33351225) q[2];
rz(-2.3836661) q[3];
sx q[3];
rz(-2.1099542) q[3];
sx q[3];
rz(-1.7655752) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
