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
rz(-2.6645339) q[0];
sx q[0];
rz(-2.8395489) q[0];
sx q[0];
rz(-1.8855236) q[0];
rz(0.37580252) q[1];
sx q[1];
rz(-1.3060275) q[1];
sx q[1];
rz(-1.1285055) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94621979) q[0];
sx q[0];
rz(-1.8502619) q[0];
sx q[0];
rz(2.165876) q[0];
x q[1];
rz(-0.62126764) q[2];
sx q[2];
rz(-1.318802) q[2];
sx q[2];
rz(2.9997343) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8339631) q[1];
sx q[1];
rz(-1.071573) q[1];
sx q[1];
rz(-2.7021033) q[1];
rz(2.5735568) q[3];
sx q[3];
rz(-0.71692335) q[3];
sx q[3];
rz(-2.4942945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6065373) q[2];
sx q[2];
rz(-2.8474319) q[2];
sx q[2];
rz(1.199022) q[2];
rz(0.27896518) q[3];
sx q[3];
rz(-1.4594892) q[3];
sx q[3];
rz(-0.038486686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8798384) q[0];
sx q[0];
rz(-2.6848875) q[0];
sx q[0];
rz(2.7754011) q[0];
rz(-1.2884033) q[1];
sx q[1];
rz(-0.90694702) q[1];
sx q[1];
rz(0.2624951) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703471) q[0];
sx q[0];
rz(-0.34929212) q[0];
sx q[0];
rz(0.55089921) q[0];
rz(-1.7967223) q[2];
sx q[2];
rz(-2.886555) q[2];
sx q[2];
rz(-2.0085427) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.69345638) q[1];
sx q[1];
rz(-0.46365204) q[1];
sx q[1];
rz(-0.69000375) q[1];
rz(-1.4743342) q[3];
sx q[3];
rz(-0.25611195) q[3];
sx q[3];
rz(-0.90463582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75258201) q[2];
sx q[2];
rz(-1.8730619) q[2];
sx q[2];
rz(2.9627723) q[2];
rz(0.016077476) q[3];
sx q[3];
rz(-0.5355081) q[3];
sx q[3];
rz(-2.2107562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7054684) q[0];
sx q[0];
rz(-3.0146764) q[0];
sx q[0];
rz(-2.5674168) q[0];
rz(-2.9899959) q[1];
sx q[1];
rz(-2.6972289) q[1];
sx q[1];
rz(-1.9134329) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235665) q[0];
sx q[0];
rz(-2.0834588) q[0];
sx q[0];
rz(0.17983564) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.068563383) q[2];
sx q[2];
rz(-0.97899635) q[2];
sx q[2];
rz(-0.14321274) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0361402) q[1];
sx q[1];
rz(-2.5220118) q[1];
sx q[1];
rz(-1.7253656) q[1];
rz(-2.853573) q[3];
sx q[3];
rz(-1.1041118) q[3];
sx q[3];
rz(-0.53190255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14153081) q[2];
sx q[2];
rz(-2.5059293) q[2];
sx q[2];
rz(2.2504811) q[2];
rz(2.7530503) q[3];
sx q[3];
rz(-1.6308234) q[3];
sx q[3];
rz(-2.8019606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6252839) q[0];
sx q[0];
rz(-0.044965222) q[0];
sx q[0];
rz(2.6594824) q[0];
rz(-2.4730143) q[1];
sx q[1];
rz(-1.6079638) q[1];
sx q[1];
rz(0.63749981) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264049) q[0];
sx q[0];
rz(-2.5026461) q[0];
sx q[0];
rz(-1.9685345) q[0];
x q[1];
rz(-0.89738557) q[2];
sx q[2];
rz(-1.6614794) q[2];
sx q[2];
rz(-0.53340837) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5920221) q[1];
sx q[1];
rz(-1.6046822) q[1];
sx q[1];
rz(3.0646044) q[1];
rz(0.25084514) q[3];
sx q[3];
rz(-1.224784) q[3];
sx q[3];
rz(0.097867997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1981226) q[2];
sx q[2];
rz(-2.1169457) q[2];
sx q[2];
rz(0.060297273) q[2];
rz(-2.8382525) q[3];
sx q[3];
rz(-3.0637432) q[3];
sx q[3];
rz(-0.2651324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1497134) q[0];
sx q[0];
rz(-2.7880221) q[0];
sx q[0];
rz(-2.7283332) q[0];
rz(-0.39516559) q[1];
sx q[1];
rz(-1.7839849) q[1];
sx q[1];
rz(-1.0023592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7089139) q[0];
sx q[0];
rz(-2.1122547) q[0];
sx q[0];
rz(-1.0546636) q[0];
rz(-1.3609028) q[2];
sx q[2];
rz(-1.6060467) q[2];
sx q[2];
rz(1.9942703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.099157028) q[1];
sx q[1];
rz(-0.71036541) q[1];
sx q[1];
rz(-1.9036681) q[1];
rz(-pi) q[2];
rz(-0.20805863) q[3];
sx q[3];
rz(-1.8866619) q[3];
sx q[3];
rz(-2.7622595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5137382) q[2];
sx q[2];
rz(-2.0444874) q[2];
sx q[2];
rz(-1.537079) q[2];
rz(1.1995992) q[3];
sx q[3];
rz(-1.1211841) q[3];
sx q[3];
rz(-2.3698923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5143249) q[0];
sx q[0];
rz(-1.5031313) q[0];
sx q[0];
rz(2.7728873) q[0];
rz(-2.3983224) q[1];
sx q[1];
rz(-2.5786046) q[1];
sx q[1];
rz(2.7875913) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5394595) q[0];
sx q[0];
rz(-1.51855) q[0];
sx q[0];
rz(0.9297855) q[0];
rz(-pi) q[1];
rz(-0.52704372) q[2];
sx q[2];
rz(-1.8436925) q[2];
sx q[2];
rz(-0.54293406) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2555775) q[1];
sx q[1];
rz(-1.7486582) q[1];
sx q[1];
rz(-2.1484003) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1844814) q[3];
sx q[3];
rz(-2.1710178) q[3];
sx q[3];
rz(-1.6627897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.02504286) q[2];
sx q[2];
rz(-1.2094689) q[2];
sx q[2];
rz(-3.1361191) q[2];
rz(0.51760751) q[3];
sx q[3];
rz(-0.36492437) q[3];
sx q[3];
rz(0.027211729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40590498) q[0];
sx q[0];
rz(-2.5820177) q[0];
sx q[0];
rz(-2.9285808) q[0];
rz(-1.4951911) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(-2.2921553) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54100688) q[0];
sx q[0];
rz(-1.6141506) q[0];
sx q[0];
rz(0.61886532) q[0];
rz(3.1108625) q[2];
sx q[2];
rz(-1.6140013) q[2];
sx q[2];
rz(1.876056) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7905459) q[1];
sx q[1];
rz(-2.4645269) q[1];
sx q[1];
rz(0.76404567) q[1];
rz(-1.8324424) q[3];
sx q[3];
rz(-0.54352409) q[3];
sx q[3];
rz(1.3268918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9269632) q[2];
sx q[2];
rz(-1.233485) q[2];
sx q[2];
rz(-2.5978973) q[2];
rz(-0.26116535) q[3];
sx q[3];
rz(-1.5800913) q[3];
sx q[3];
rz(-0.15682596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1866622) q[0];
sx q[0];
rz(-1.0533227) q[0];
sx q[0];
rz(-0.18761158) q[0];
rz(2.018351) q[1];
sx q[1];
rz(-0.76485991) q[1];
sx q[1];
rz(-2.979028) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.150972) q[0];
sx q[0];
rz(-0.91998581) q[0];
sx q[0];
rz(1.9810173) q[0];
rz(0.17485042) q[2];
sx q[2];
rz(-1.0533337) q[2];
sx q[2];
rz(-2.3558921) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9342207) q[1];
sx q[1];
rz(-1.6160864) q[1];
sx q[1];
rz(0.95691935) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7361338) q[3];
sx q[3];
rz(-2.5389034) q[3];
sx q[3];
rz(-1.7053982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8673914) q[2];
sx q[2];
rz(-2.826639) q[2];
sx q[2];
rz(-0.20142041) q[2];
rz(-2.6650186) q[3];
sx q[3];
rz(-1.3786992) q[3];
sx q[3];
rz(0.99982888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9607214) q[0];
sx q[0];
rz(-2.9044386) q[0];
sx q[0];
rz(0.57688212) q[0];
rz(-2.8374529) q[1];
sx q[1];
rz(-1.1241333) q[1];
sx q[1];
rz(-2.0374128) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7542412) q[0];
sx q[0];
rz(-1.1925444) q[0];
sx q[0];
rz(0.095281466) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18851243) q[2];
sx q[2];
rz(-0.57706632) q[2];
sx q[2];
rz(0.12332502) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4889646) q[1];
sx q[1];
rz(-1.9581989) q[1];
sx q[1];
rz(-0.94469686) q[1];
rz(-pi) q[2];
rz(-1.1796168) q[3];
sx q[3];
rz(-1.5857539) q[3];
sx q[3];
rz(0.43128951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96678174) q[2];
sx q[2];
rz(-2.8947783) q[2];
sx q[2];
rz(0.68320572) q[2];
rz(1.5481868) q[3];
sx q[3];
rz(-2.4631409) q[3];
sx q[3];
rz(2.9233176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9504647) q[0];
sx q[0];
rz(-2.2948965) q[0];
sx q[0];
rz(-2.6306613) q[0];
rz(-1.9966985) q[1];
sx q[1];
rz(-1.1868008) q[1];
sx q[1];
rz(1.5085545) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.975113) q[0];
sx q[0];
rz(-0.60752019) q[0];
sx q[0];
rz(-0.44618531) q[0];
rz(-1.0867003) q[2];
sx q[2];
rz(-1.3333313) q[2];
sx q[2];
rz(0.1456085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3737693) q[1];
sx q[1];
rz(-1.1255742) q[1];
sx q[1];
rz(1.7675464) q[1];
rz(-pi) q[2];
rz(2.4436117) q[3];
sx q[3];
rz(-2.7578286) q[3];
sx q[3];
rz(2.021054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42023429) q[2];
sx q[2];
rz(-0.63174641) q[2];
sx q[2];
rz(-1.7101804) q[2];
rz(-0.015627705) q[3];
sx q[3];
rz(-1.2651919) q[3];
sx q[3];
rz(0.50610745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.90841993) q[0];
sx q[0];
rz(-1.7411727) q[0];
sx q[0];
rz(2.2199051) q[0];
rz(-1.5064916) q[1];
sx q[1];
rz(-1.9252621) q[1];
sx q[1];
rz(2.6797934) q[1];
rz(2.2057381) q[2];
sx q[2];
rz(-2.6781185) q[2];
sx q[2];
rz(0.92739633) q[2];
rz(1.1453859) q[3];
sx q[3];
rz(-0.085193188) q[3];
sx q[3];
rz(-1.9151158) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
