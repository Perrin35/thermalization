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
rz(0.25207818) q[0];
sx q[0];
rz(3.7171465) q[0];
sx q[0];
rz(9.3024749) q[0];
rz(-2.0343434) q[1];
sx q[1];
rz(-2.1802433) q[1];
sx q[1];
rz(-3.1032739) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0487028) q[0];
sx q[0];
rz(-2.1445591) q[0];
sx q[0];
rz(2.696373) q[0];
x q[1];
rz(1.6898481) q[2];
sx q[2];
rz(-1.6887293) q[2];
sx q[2];
rz(-0.32334194) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3342383) q[1];
sx q[1];
rz(-2.1209201) q[1];
sx q[1];
rz(-0.38205876) q[1];
rz(1.5913209) q[3];
sx q[3];
rz(-1.8181268) q[3];
sx q[3];
rz(2.1065245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4965839) q[2];
sx q[2];
rz(-1.4270447) q[2];
sx q[2];
rz(1.4487779) q[2];
rz(-2.951176) q[3];
sx q[3];
rz(-1.0551635) q[3];
sx q[3];
rz(-2.9922488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9377624) q[0];
sx q[0];
rz(-1.8730524) q[0];
sx q[0];
rz(2.8835836) q[0];
rz(1.5915271) q[1];
sx q[1];
rz(-1.9015046) q[1];
sx q[1];
rz(2.9749427) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23388966) q[0];
sx q[0];
rz(-2.5512716) q[0];
sx q[0];
rz(-2.9695244) q[0];
x q[1];
rz(-1.2857976) q[2];
sx q[2];
rz(-1.1995595) q[2];
sx q[2];
rz(1.2165537) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0479134) q[1];
sx q[1];
rz(-2.4166345) q[1];
sx q[1];
rz(2.9525501) q[1];
rz(-pi) q[2];
rz(-1.0051425) q[3];
sx q[3];
rz(-1.0962624) q[3];
sx q[3];
rz(2.8090854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0698645) q[2];
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
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.449618) q[0];
sx q[0];
rz(-1.6716577) q[0];
sx q[0];
rz(3.1396507) q[0];
rz(-3.107403) q[1];
sx q[1];
rz(-1.9296153) q[1];
sx q[1];
rz(-1.5431822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6337834) q[0];
sx q[0];
rz(-2.6144321) q[0];
sx q[0];
rz(0.5306923) q[0];
rz(2.7091647) q[2];
sx q[2];
rz(-0.78410599) q[2];
sx q[2];
rz(2.6570005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0405827) q[1];
sx q[1];
rz(-1.9242745) q[1];
sx q[1];
rz(1.5291052) q[1];
x q[2];
rz(0.9228306) q[3];
sx q[3];
rz(-1.3746972) q[3];
sx q[3];
rz(0.01250532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2446642) q[2];
sx q[2];
rz(-2.7074773) q[2];
sx q[2];
rz(-1.0373235) q[2];
rz(-2.0364929) q[3];
sx q[3];
rz(-1.6048071) q[3];
sx q[3];
rz(1.1387811) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8741375) q[0];
sx q[0];
rz(-0.48478165) q[0];
sx q[0];
rz(1.7304035) q[0];
rz(-0.63938582) q[1];
sx q[1];
rz(-1.4566028) q[1];
sx q[1];
rz(0.04714084) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5306803) q[0];
sx q[0];
rz(-2.147104) q[0];
sx q[0];
rz(1.9082101) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7732255) q[2];
sx q[2];
rz(-0.50058156) q[2];
sx q[2];
rz(0.49886242) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85002181) q[1];
sx q[1];
rz(-0.98660368) q[1];
sx q[1];
rz(2.7136346) q[1];
rz(2.8465038) q[3];
sx q[3];
rz(-0.86936823) q[3];
sx q[3];
rz(1.6158582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3186657) q[2];
sx q[2];
rz(-1.9959799) q[2];
sx q[2];
rz(1.9226496) q[2];
rz(-0.31560358) q[3];
sx q[3];
rz(-3.0867519) q[3];
sx q[3];
rz(-0.1489197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8197935) q[0];
sx q[0];
rz(-2.5892374) q[0];
sx q[0];
rz(-2.7929982) q[0];
rz(1.2527342) q[1];
sx q[1];
rz(-1.1628954) q[1];
sx q[1];
rz(-1.8399651) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0785261) q[0];
sx q[0];
rz(-1.491703) q[0];
sx q[0];
rz(2.3344759) q[0];
x q[1];
rz(1.9867861) q[2];
sx q[2];
rz(-1.0495397) q[2];
sx q[2];
rz(2.398587) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5036936) q[1];
sx q[1];
rz(-2.854373) q[1];
sx q[1];
rz(2.2566608) q[1];
rz(-pi) q[2];
rz(1.0200649) q[3];
sx q[3];
rz(-1.2675084) q[3];
sx q[3];
rz(2.1316949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18315135) q[2];
sx q[2];
rz(-0.30854598) q[2];
sx q[2];
rz(1.4754254) q[2];
rz(-0.685855) q[3];
sx q[3];
rz(-2.1619022) q[3];
sx q[3];
rz(0.53387749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.1118065) q[0];
sx q[0];
rz(-0.15203467) q[0];
sx q[0];
rz(-1.6492122) q[0];
rz(-2.1624508) q[1];
sx q[1];
rz(-1.6056332) q[1];
sx q[1];
rz(-1.2243366) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6394854) q[0];
sx q[0];
rz(-1.7687135) q[0];
sx q[0];
rz(-1.5396126) q[0];
rz(-2.5923205) q[2];
sx q[2];
rz(-1.7584929) q[2];
sx q[2];
rz(2.930738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30019125) q[1];
sx q[1];
rz(-0.61392861) q[1];
sx q[1];
rz(1.5924686) q[1];
x q[2];
rz(-2.79412) q[3];
sx q[3];
rz(-2.1712448) q[3];
sx q[3];
rz(2.1640409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7225723) q[2];
sx q[2];
rz(-1.4973065) q[2];
sx q[2];
rz(-2.2255619) q[2];
rz(1.3845059) q[3];
sx q[3];
rz(-1.2576831) q[3];
sx q[3];
rz(2.4345583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0093805669) q[0];
sx q[0];
rz(-0.051932422) q[0];
sx q[0];
rz(-0.95426553) q[0];
rz(0.43680278) q[1];
sx q[1];
rz(-1.5780459) q[1];
sx q[1];
rz(0.11016914) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018579114) q[0];
sx q[0];
rz(-1.8948104) q[0];
sx q[0];
rz(1.3186245) q[0];
rz(2.4442441) q[2];
sx q[2];
rz(-1.9983872) q[2];
sx q[2];
rz(-0.22874895) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4272144) q[1];
sx q[1];
rz(-2.7958779) q[1];
sx q[1];
rz(-1.5580672) q[1];
rz(-pi) q[2];
rz(-2.8875868) q[3];
sx q[3];
rz(-2.0585103) q[3];
sx q[3];
rz(-0.28764492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7097912) q[2];
sx q[2];
rz(-0.0063889901) q[2];
sx q[2];
rz(0.94376454) q[2];
rz(-0.56378311) q[3];
sx q[3];
rz(-1.0958593) q[3];
sx q[3];
rz(-1.110466) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3540045) q[0];
sx q[0];
rz(-2.8546951) q[0];
sx q[0];
rz(-0.34550825) q[0];
rz(-1.0954789) q[1];
sx q[1];
rz(-1.4778719) q[1];
sx q[1];
rz(-1.5617721) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0559383) q[0];
sx q[0];
rz(-1.602218) q[0];
sx q[0];
rz(-2.3532193) q[0];
rz(-pi) q[1];
rz(1.6510294) q[2];
sx q[2];
rz(-0.65872619) q[2];
sx q[2];
rz(-1.6541437) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6553538) q[1];
sx q[1];
rz(-1.2062688) q[1];
sx q[1];
rz(-2.2097387) q[1];
rz(-pi) q[2];
rz(-0.49786477) q[3];
sx q[3];
rz(-1.5427054) q[3];
sx q[3];
rz(1.3974578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76274189) q[2];
sx q[2];
rz(-2.8801535) q[2];
sx q[2];
rz(1.4307107) q[2];
rz(2.6070969) q[3];
sx q[3];
rz(-1.675324) q[3];
sx q[3];
rz(-0.9526332) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5296103) q[0];
sx q[0];
rz(-3.072325) q[0];
sx q[0];
rz(-1.3105422) q[0];
rz(-0.88690859) q[1];
sx q[1];
rz(-0.84019089) q[1];
sx q[1];
rz(-1.8720253) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5909605) q[0];
sx q[0];
rz(-0.3754339) q[0];
sx q[0];
rz(1.9272789) q[0];
x q[1];
rz(-1.194095) q[2];
sx q[2];
rz(-1.4262876) q[2];
sx q[2];
rz(-1.2928177) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7901526) q[1];
sx q[1];
rz(-1.8428486) q[1];
sx q[1];
rz(2.9012009) q[1];
x q[2];
rz(-0.19143243) q[3];
sx q[3];
rz(-0.20729724) q[3];
sx q[3];
rz(-1.3216922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5294007) q[2];
sx q[2];
rz(-1.3357013) q[2];
sx q[2];
rz(2.9086435) q[2];
rz(-1.285078) q[3];
sx q[3];
rz(-0.67684567) q[3];
sx q[3];
rz(2.7555833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9832298) q[0];
sx q[0];
rz(-0.60929275) q[0];
sx q[0];
rz(3.0443211) q[0];
rz(-2.3376047) q[1];
sx q[1];
rz(-2.4318047) q[1];
sx q[1];
rz(2.9248617) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6410206) q[0];
sx q[0];
rz(-1.7413229) q[0];
sx q[0];
rz(-3.0855623) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9998776) q[2];
sx q[2];
rz(-1.6167621) q[2];
sx q[2];
rz(-1.3675929) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35516675) q[1];
sx q[1];
rz(-2.2192419) q[1];
sx q[1];
rz(2.4207958) q[1];
x q[2];
rz(1.9166462) q[3];
sx q[3];
rz(-2.2035363) q[3];
sx q[3];
rz(-0.94055292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29297605) q[2];
sx q[2];
rz(-2.8605707) q[2];
sx q[2];
rz(-0.49523735) q[2];
rz(1.5589335) q[3];
sx q[3];
rz(-2.0571183) q[3];
sx q[3];
rz(2.9787279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0611298) q[0];
sx q[0];
rz(-1.5123788) q[0];
sx q[0];
rz(2.6176591) q[0];
rz(-1.0864661) q[1];
sx q[1];
rz(-2.6230984) q[1];
sx q[1];
rz(-2.7883504) q[1];
rz(-3.1343717) q[2];
sx q[2];
rz(-1.6025958) q[2];
sx q[2];
rz(-0.64000426) q[2];
rz(-1.6752406) q[3];
sx q[3];
rz(-2.4267195) q[3];
sx q[3];
rz(0.60703312) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
