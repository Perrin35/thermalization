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
rz(0.47705874) q[0];
sx q[0];
rz(-0.30204371) q[0];
sx q[0];
rz(1.8855236) q[0];
rz(0.37580252) q[1];
sx q[1];
rz(1.8355651) q[1];
sx q[1];
rz(10.553283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9039584) q[0];
sx q[0];
rz(-2.4914143) q[0];
sx q[0];
rz(2.0439434) q[0];
rz(2.520325) q[2];
sx q[2];
rz(-1.318802) q[2];
sx q[2];
rz(-0.14185837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0417718) q[1];
sx q[1];
rz(-1.1879634) q[1];
sx q[1];
rz(-2.113093) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1323186) q[3];
sx q[3];
rz(-0.98377675) q[3];
sx q[3];
rz(-1.7917716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5350554) q[2];
sx q[2];
rz(-2.8474319) q[2];
sx q[2];
rz(1.9425707) q[2];
rz(-2.8626275) q[3];
sx q[3];
rz(-1.4594892) q[3];
sx q[3];
rz(-0.038486686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8798384) q[0];
sx q[0];
rz(-2.6848875) q[0];
sx q[0];
rz(-2.7754011) q[0];
rz(-1.8531894) q[1];
sx q[1];
rz(-2.2346456) q[1];
sx q[1];
rz(-2.8790976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27124559) q[0];
sx q[0];
rz(-2.7923005) q[0];
sx q[0];
rz(0.55089921) q[0];
rz(-pi) q[1];
rz(-3.0832564) q[2];
sx q[2];
rz(-1.8192141) q[2];
sx q[2];
rz(-2.2417541) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9001635) q[1];
sx q[1];
rz(-1.2821322) q[1];
sx q[1];
rz(2.7735387) q[1];
rz(3.1163773) q[3];
sx q[3];
rz(-1.3159015) q[3];
sx q[3];
rz(-1.0043292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3890106) q[2];
sx q[2];
rz(-1.8730619) q[2];
sx q[2];
rz(-0.1788204) q[2];
rz(0.016077476) q[3];
sx q[3];
rz(-0.5355081) q[3];
sx q[3];
rz(-2.2107562) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43612424) q[0];
sx q[0];
rz(-0.12691623) q[0];
sx q[0];
rz(-0.57417589) q[0];
rz(-2.9899959) q[1];
sx q[1];
rz(-2.6972289) q[1];
sx q[1];
rz(1.2281598) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7235665) q[0];
sx q[0];
rz(-2.0834588) q[0];
sx q[0];
rz(2.961757) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0730293) q[2];
sx q[2];
rz(-0.97899635) q[2];
sx q[2];
rz(-2.9983799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2252402) q[1];
sx q[1];
rz(-2.181899) q[1];
sx q[1];
rz(3.0322187) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0841885) q[3];
sx q[3];
rz(-0.54277615) q[3];
sx q[3];
rz(2.0273939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.14153081) q[2];
sx q[2];
rz(-0.63566339) q[2];
sx q[2];
rz(-2.2504811) q[2];
rz(2.7530503) q[3];
sx q[3];
rz(-1.5107692) q[3];
sx q[3];
rz(-0.33963206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6252839) q[0];
sx q[0];
rz(-3.0966274) q[0];
sx q[0];
rz(-2.6594824) q[0];
rz(-0.66857839) q[1];
sx q[1];
rz(-1.5336288) q[1];
sx q[1];
rz(-2.5040928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26698819) q[0];
sx q[0];
rz(-0.98867304) q[0];
sx q[0];
rz(0.28018392) q[0];
rz(-pi) q[1];
rz(0.11580616) q[2];
sx q[2];
rz(-0.90065946) q[2];
sx q[2];
rz(-2.0320803) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60739775) q[1];
sx q[1];
rz(-3.0574905) q[1];
sx q[1];
rz(-2.7264595) q[1];
x q[2];
rz(-1.927078) q[3];
sx q[3];
rz(-1.8064835) q[3];
sx q[3];
rz(-1.5819777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1981226) q[2];
sx q[2];
rz(-1.024647) q[2];
sx q[2];
rz(3.0812954) q[2];
rz(-0.30334011) q[3];
sx q[3];
rz(-3.0637432) q[3];
sx q[3];
rz(-2.8764603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1497134) q[0];
sx q[0];
rz(-2.7880221) q[0];
sx q[0];
rz(2.7283332) q[0];
rz(2.7464271) q[1];
sx q[1];
rz(-1.3576077) q[1];
sx q[1];
rz(-2.1392335) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7089139) q[0];
sx q[0];
rz(-1.029338) q[0];
sx q[0];
rz(-2.0869291) q[0];
rz(-pi) q[1];
rz(-1.4031314) q[2];
sx q[2];
rz(-2.9288026) q[2];
sx q[2];
rz(-2.8820702) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2152948) q[1];
sx q[1];
rz(-1.7855254) q[1];
sx q[1];
rz(0.88820998) q[1];
rz(1.2484364) q[3];
sx q[3];
rz(-1.768422) q[3];
sx q[3];
rz(1.2569497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62785441) q[2];
sx q[2];
rz(-1.0971053) q[2];
sx q[2];
rz(-1.537079) q[2];
rz(1.9419935) q[3];
sx q[3];
rz(-2.0204085) q[3];
sx q[3];
rz(-2.3698923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5143249) q[0];
sx q[0];
rz(-1.5031313) q[0];
sx q[0];
rz(2.7728873) q[0];
rz(-0.74327028) q[1];
sx q[1];
rz(-0.5629881) q[1];
sx q[1];
rz(-0.35400131) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1339851) q[0];
sx q[0];
rz(-0.93080321) q[0];
sx q[0];
rz(-3.0764393) q[0];
x q[1];
rz(2.6145489) q[2];
sx q[2];
rz(-1.8436925) q[2];
sx q[2];
rz(-0.54293406) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1917369) q[1];
sx q[1];
rz(-2.5402148) q[1];
sx q[1];
rz(-1.2527501) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50320585) q[3];
sx q[3];
rz(-0.70069289) q[3];
sx q[3];
rz(-2.1029187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.02504286) q[2];
sx q[2];
rz(-1.2094689) q[2];
sx q[2];
rz(-0.0054736007) q[2];
rz(-0.51760751) q[3];
sx q[3];
rz(-0.36492437) q[3];
sx q[3];
rz(-0.027211729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7356877) q[0];
sx q[0];
rz(-2.5820177) q[0];
sx q[0];
rz(-2.9285808) q[0];
rz(-1.4951911) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(-2.2921553) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0606468) q[0];
sx q[0];
rz(-2.1889926) q[0];
sx q[0];
rz(1.5175876) q[0];
rz(-0.95291887) q[2];
sx q[2];
rz(-0.053013405) q[2];
sx q[2];
rz(1.8840781) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7905459) q[1];
sx q[1];
rz(-0.67706579) q[1];
sx q[1];
rz(0.76404567) q[1];
x q[2];
rz(1.3091503) q[3];
sx q[3];
rz(-2.5980686) q[3];
sx q[3];
rz(1.8147008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21462943) q[2];
sx q[2];
rz(-1.233485) q[2];
sx q[2];
rz(2.5978973) q[2];
rz(0.26116535) q[3];
sx q[3];
rz(-1.5800913) q[3];
sx q[3];
rz(-2.9847667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1866622) q[0];
sx q[0];
rz(-1.0533227) q[0];
sx q[0];
rz(2.9539811) q[0];
rz(-1.1232417) q[1];
sx q[1];
rz(-0.76485991) q[1];
sx q[1];
rz(0.16256464) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32255781) q[0];
sx q[0];
rz(-1.8936689) q[0];
sx q[0];
rz(-2.4486008) q[0];
rz(-pi) q[1];
rz(2.0948903) q[2];
sx q[2];
rz(-1.419029) q[2];
sx q[2];
rz(-2.2693315) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8423262) q[1];
sx q[1];
rz(-0.61533058) q[1];
sx q[1];
rz(-1.6493114) q[1];
rz(-1.4054589) q[3];
sx q[3];
rz(-2.5389034) q[3];
sx q[3];
rz(-1.7053982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2742013) q[2];
sx q[2];
rz(-2.826639) q[2];
sx q[2];
rz(2.9401722) q[2];
rz(-2.6650186) q[3];
sx q[3];
rz(-1.7628935) q[3];
sx q[3];
rz(-0.99982888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1808712) q[0];
sx q[0];
rz(-2.9044386) q[0];
sx q[0];
rz(2.5647105) q[0];
rz(-0.30413973) q[1];
sx q[1];
rz(-1.1241333) q[1];
sx q[1];
rz(2.0374128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21872422) q[0];
sx q[0];
rz(-1.6593242) q[0];
sx q[0];
rz(-1.1909816) q[0];
rz(-pi) q[1];
rz(-0.18851243) q[2];
sx q[2];
rz(-2.5645263) q[2];
sx q[2];
rz(-0.12332502) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34852178) q[1];
sx q[1];
rz(-2.1442815) q[1];
sx q[1];
rz(0.46646935) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.125413) q[3];
sx q[3];
rz(-1.1796629) q[3];
sx q[3];
rz(-1.1333381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1748109) q[2];
sx q[2];
rz(-2.8947783) q[2];
sx q[2];
rz(0.68320572) q[2];
rz(1.5481868) q[3];
sx q[3];
rz(-0.67845172) q[3];
sx q[3];
rz(0.21827503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19112793) q[0];
sx q[0];
rz(-0.8466962) q[0];
sx q[0];
rz(0.51093131) q[0];
rz(-1.1448942) q[1];
sx q[1];
rz(-1.1868008) q[1];
sx q[1];
rz(-1.5085545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7785698) q[0];
sx q[0];
rz(-1.3219044) q[0];
sx q[0];
rz(-2.581434) q[0];
x q[1];
rz(1.0912446) q[2];
sx q[2];
rz(-2.6065718) q[2];
sx q[2];
rz(-1.0045235) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.252987) q[1];
sx q[1];
rz(-1.7481511) q[1];
sx q[1];
rz(-2.6887522) q[1];
rz(-pi) q[2];
rz(-1.8246966) q[3];
sx q[3];
rz(-1.2798556) q[3];
sx q[3];
rz(-1.8559141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7213584) q[2];
sx q[2];
rz(-0.63174641) q[2];
sx q[2];
rz(1.4314123) q[2];
rz(-0.015627705) q[3];
sx q[3];
rz(-1.2651919) q[3];
sx q[3];
rz(0.50610745) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2331727) q[0];
sx q[0];
rz(-1.40042) q[0];
sx q[0];
rz(-0.9216876) q[0];
rz(1.635101) q[1];
sx q[1];
rz(-1.9252621) q[1];
sx q[1];
rz(2.6797934) q[1];
rz(-2.8534081) q[2];
sx q[2];
rz(-1.202604) q[2];
sx q[2];
rz(1.6163608) q[2];
rz(1.9962068) q[3];
sx q[3];
rz(-3.0563995) q[3];
sx q[3];
rz(1.2264768) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
