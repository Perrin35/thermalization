OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(-2.2280333) q[0];
sx q[0];
rz(1.7295184) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(5.6875416) q[1];
sx q[1];
rz(7.7653801) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1887814) q[0];
sx q[0];
rz(-0.45438284) q[0];
sx q[0];
rz(0.36034583) q[0];
x q[1];
rz(2.8843845) q[2];
sx q[2];
rz(-1.4423941) q[2];
sx q[2];
rz(0.53127015) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2130148) q[1];
sx q[1];
rz(-1.1482114) q[1];
sx q[1];
rz(1.0282474) q[1];
x q[2];
rz(-1.4115303) q[3];
sx q[3];
rz(-2.5730238) q[3];
sx q[3];
rz(1.1390613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1564864) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(2.1872897) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(-1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(-1.5401309) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(0.96347934) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.119495) q[0];
sx q[0];
rz(-2.1848328) q[0];
sx q[0];
rz(3.1387781) q[0];
x q[1];
rz(-1.0785525) q[2];
sx q[2];
rz(-0.62765593) q[2];
sx q[2];
rz(-0.17222675) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.023167921) q[1];
sx q[1];
rz(-1.2512565) q[1];
sx q[1];
rz(2.097514) q[1];
x q[2];
rz(2.1727932) q[3];
sx q[3];
rz(-0.41799823) q[3];
sx q[3];
rz(-2.9434162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(2.3965805) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(-2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(-2.0939317) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(0.55999666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3968351) q[0];
sx q[0];
rz(-2.7042537) q[0];
sx q[0];
rz(-2.0646981) q[0];
rz(-3.0736018) q[2];
sx q[2];
rz(-2.8376841) q[2];
sx q[2];
rz(1.4246724) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33369505) q[1];
sx q[1];
rz(-1.9296608) q[1];
sx q[1];
rz(1.7829423) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3476944) q[3];
sx q[3];
rz(-2.5069935) q[3];
sx q[3];
rz(-1.8397699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3893163) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(2.9690572) q[2];
rz(-2.1595188) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1383706) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(0.28451434) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(1.2987312) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5643698) q[0];
sx q[0];
rz(-1.3457314) q[0];
sx q[0];
rz(1.0610915) q[0];
rz(-pi) q[1];
rz(-1.5721333) q[2];
sx q[2];
rz(-2.3751811) q[2];
sx q[2];
rz(0.02174755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.270077) q[1];
sx q[1];
rz(-1.3844826) q[1];
sx q[1];
rz(2.1427076) q[1];
x q[2];
rz(1.2426162) q[3];
sx q[3];
rz(-1.6771852) q[3];
sx q[3];
rz(-2.1387517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(2.7569125) q[2];
rz(2.3875333) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(1.4543021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6032747) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(-2.8447661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9223698) q[0];
sx q[0];
rz(-0.57225675) q[0];
sx q[0];
rz(0.072364307) q[0];
rz(-pi) q[1];
rz(-2.9752521) q[2];
sx q[2];
rz(-2.3918249) q[2];
sx q[2];
rz(-1.2321842) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0592812) q[1];
sx q[1];
rz(-2.0356405) q[1];
sx q[1];
rz(2.6962198) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4279168) q[3];
sx q[3];
rz(-1.5034961) q[3];
sx q[3];
rz(2.0590559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(0.39247593) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(-2.8241482) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1257989) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(2.3902067) q[0];
rz(-1.8136576) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(-0.60633916) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6316815) q[0];
sx q[0];
rz(-1.5874377) q[0];
sx q[0];
rz(3.095754) q[0];
x q[1];
rz(-1.0688071) q[2];
sx q[2];
rz(-2.398218) q[2];
sx q[2];
rz(2.5943287) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6341083) q[1];
sx q[1];
rz(-1.023523) q[1];
sx q[1];
rz(2.9299111) q[1];
x q[2];
rz(1.243152) q[3];
sx q[3];
rz(-2.9122105) q[3];
sx q[3];
rz(2.4250507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63885826) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(-2.0022557) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6095603) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(-2.6334921) q[0];
rz(-1.5787026) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(2.3513444) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7471874) q[0];
sx q[0];
rz(-0.98441511) q[0];
sx q[0];
rz(1.2140843) q[0];
rz(-0.21364084) q[2];
sx q[2];
rz(-1.5822615) q[2];
sx q[2];
rz(-1.2917047) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0223169) q[1];
sx q[1];
rz(-1.9609309) q[1];
sx q[1];
rz(-1.8632061) q[1];
x q[2];
rz(-1.1931476) q[3];
sx q[3];
rz(-2.1834063) q[3];
sx q[3];
rz(2.2539504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.053085176) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(-1.6648071) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(-3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168032) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(-2.0943663) q[0];
rz(-2.5324902) q[1];
sx q[1];
rz(-1.4139688) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7097276) q[0];
sx q[0];
rz(-2.1301529) q[0];
sx q[0];
rz(3.1064242) q[0];
x q[1];
rz(0.65865626) q[2];
sx q[2];
rz(-0.63630644) q[2];
sx q[2];
rz(3.0898526) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3125004) q[1];
sx q[1];
rz(-1.6626076) q[1];
sx q[1];
rz(-1.5822259) q[1];
rz(-pi) q[2];
rz(-2.8453313) q[3];
sx q[3];
rz(-2.1564266) q[3];
sx q[3];
rz(2.0010009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(0.88225538) q[2];
rz(1.7404209) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27154487) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(-1.4260938) q[0];
rz(-0.081461279) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(-2.5833599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0840069) q[0];
sx q[0];
rz(-1.1450197) q[0];
sx q[0];
rz(-2.7140679) q[0];
rz(-1.4633281) q[2];
sx q[2];
rz(-1.5170013) q[2];
sx q[2];
rz(1.6966284) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12802943) q[1];
sx q[1];
rz(-2.3131144) q[1];
sx q[1];
rz(0.086704266) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34096034) q[3];
sx q[3];
rz(-0.51135671) q[3];
sx q[3];
rz(-2.2380059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50487173) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.5378392) q[0];
rz(2.3161855) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-2.6182981) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16738811) q[0];
sx q[0];
rz(-0.61199576) q[0];
sx q[0];
rz(-0.1083072) q[0];
rz(0.9988437) q[2];
sx q[2];
rz(-1.6076644) q[2];
sx q[2];
rz(-1.8429304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20269468) q[1];
sx q[1];
rz(-1.885186) q[1];
sx q[1];
rz(-0.37757341) q[1];
rz(2.8966122) q[3];
sx q[3];
rz(-1.3462726) q[3];
sx q[3];
rz(-2.6905439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.837073) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(-2.8722897) q[2];
rz(0.4942016) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(2.0555029) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3257278) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(-1.4670463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(0.011209839) q[2];
sx q[2];
rz(-1.8335473) q[2];
sx q[2];
rz(2.0469472) q[2];
rz(2.3446463) q[3];
sx q[3];
rz(-1.7399825) q[3];
sx q[3];
rz(0.83132838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
