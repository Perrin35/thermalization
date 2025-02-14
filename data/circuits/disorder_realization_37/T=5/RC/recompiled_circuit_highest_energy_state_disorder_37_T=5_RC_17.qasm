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
rz(-0.77684075) q[0];
sx q[0];
rz(2.2651894) q[0];
sx q[0];
rz(9.6776008) q[0];
rz(1.1480992) q[1];
sx q[1];
rz(-2.2263081) q[1];
sx q[1];
rz(-0.80438703) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1348477) q[0];
sx q[0];
rz(-0.50323707) q[0];
sx q[0];
rz(-2.2630967) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6635799) q[2];
sx q[2];
rz(-0.4768663) q[2];
sx q[2];
rz(-3.06682) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8411257) q[1];
sx q[1];
rz(-1.9411095) q[1];
sx q[1];
rz(1.9410067) q[1];
x q[2];
rz(-0.34396307) q[3];
sx q[3];
rz(-1.4620145) q[3];
sx q[3];
rz(0.0057084486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.09314166) q[2];
sx q[2];
rz(-2.1790049) q[2];
sx q[2];
rz(0.87240458) q[2];
rz(-0.33556542) q[3];
sx q[3];
rz(-1.7458956) q[3];
sx q[3];
rz(-0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63808477) q[0];
sx q[0];
rz(-0.26569772) q[0];
sx q[0];
rz(-2.9124394) q[0];
rz(0.19042641) q[1];
sx q[1];
rz(-1.3399905) q[1];
sx q[1];
rz(-2.4280039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3772954) q[0];
sx q[0];
rz(-0.72749253) q[0];
sx q[0];
rz(1.1614252) q[0];
rz(-pi) q[1];
rz(-0.83723703) q[2];
sx q[2];
rz(-1.1568489) q[2];
sx q[2];
rz(-0.6809823) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5428764) q[1];
sx q[1];
rz(-2.2440254) q[1];
sx q[1];
rz(1.190541) q[1];
rz(-pi) q[2];
rz(2.9715367) q[3];
sx q[3];
rz(-0.77203686) q[3];
sx q[3];
rz(2.7626729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68636346) q[2];
sx q[2];
rz(-2.8296622) q[2];
sx q[2];
rz(0.4757821) q[2];
rz(1.0698211) q[3];
sx q[3];
rz(-1.0082303) q[3];
sx q[3];
rz(-0.76655918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0692724) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(-1.9092165) q[0];
rz(2.6909289) q[1];
sx q[1];
rz(-1.9602937) q[1];
sx q[1];
rz(-2.252069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4159195) q[0];
sx q[0];
rz(-1.1664205) q[0];
sx q[0];
rz(-0.63979353) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1289775) q[2];
sx q[2];
rz(-1.106749) q[2];
sx q[2];
rz(-1.9249664) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0625713) q[1];
sx q[1];
rz(-1.7815184) q[1];
sx q[1];
rz(-0.91075588) q[1];
rz(-pi) q[2];
rz(1.9185478) q[3];
sx q[3];
rz(-1.9257716) q[3];
sx q[3];
rz(0.53670151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4543317) q[2];
sx q[2];
rz(-0.4250409) q[2];
sx q[2];
rz(-1.8827776) q[2];
rz(1.2339633) q[3];
sx q[3];
rz(-1.1326658) q[3];
sx q[3];
rz(3.1335355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71037978) q[0];
sx q[0];
rz(-2.2358535) q[0];
sx q[0];
rz(1.4696962) q[0];
rz(-1.9944893) q[1];
sx q[1];
rz(-1.0018145) q[1];
sx q[1];
rz(2.240644) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6916136) q[0];
sx q[0];
rz(-1.8834582) q[0];
sx q[0];
rz(1.6603966) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9587165) q[2];
sx q[2];
rz(-1.0462073) q[2];
sx q[2];
rz(-0.9768578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9228153) q[1];
sx q[1];
rz(-0.73000693) q[1];
sx q[1];
rz(-0.14132146) q[1];
x q[2];
rz(-2.4086558) q[3];
sx q[3];
rz(-1.6896491) q[3];
sx q[3];
rz(0.87866966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49742302) q[2];
sx q[2];
rz(-1.3966509) q[2];
sx q[2];
rz(-2.4141342) q[2];
rz(-0.71508956) q[3];
sx q[3];
rz(-0.46052027) q[3];
sx q[3];
rz(-2.6434744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6775976) q[0];
sx q[0];
rz(-1.7514739) q[0];
sx q[0];
rz(0.736262) q[0];
rz(2.5212506) q[1];
sx q[1];
rz(-0.54519975) q[1];
sx q[1];
rz(-1.5249407) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0686091) q[0];
sx q[0];
rz(-1.5155547) q[0];
sx q[0];
rz(1.5609571) q[0];
rz(-2.1651046) q[2];
sx q[2];
rz(-0.39820489) q[2];
sx q[2];
rz(-1.9081685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8680537) q[1];
sx q[1];
rz(-2.2333849) q[1];
sx q[1];
rz(2.6602547) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3297263) q[3];
sx q[3];
rz(-0.97211876) q[3];
sx q[3];
rz(-2.6177399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5357431) q[2];
sx q[2];
rz(-2.3652786) q[2];
sx q[2];
rz(2.1898451) q[2];
rz(2.7239299) q[3];
sx q[3];
rz(-0.2499191) q[3];
sx q[3];
rz(1.1550268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8010537) q[0];
sx q[0];
rz(-1.8908353) q[0];
sx q[0];
rz(-2.387555) q[0];
rz(-0.89318371) q[1];
sx q[1];
rz(-1.040753) q[1];
sx q[1];
rz(0.16597861) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18046019) q[0];
sx q[0];
rz(-1.6749973) q[0];
sx q[0];
rz(-0.68485884) q[0];
rz(-pi) q[1];
rz(-0.31678172) q[2];
sx q[2];
rz(-2.6133279) q[2];
sx q[2];
rz(-1.007391) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2365723) q[1];
sx q[1];
rz(-2.330594) q[1];
sx q[1];
rz(-0.6375957) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2955722) q[3];
sx q[3];
rz(-1.7050192) q[3];
sx q[3];
rz(1.272066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37504998) q[2];
sx q[2];
rz(-1.1376209) q[2];
sx q[2];
rz(-3.1362015) q[2];
rz(3.052875) q[3];
sx q[3];
rz(-1.5327449) q[3];
sx q[3];
rz(-0.79571342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8832815) q[0];
sx q[0];
rz(-1.2092051) q[0];
sx q[0];
rz(0.2051556) q[0];
rz(0.908665) q[1];
sx q[1];
rz(-0.87037194) q[1];
sx q[1];
rz(-0.044513449) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7617247) q[0];
sx q[0];
rz(-1.4141448) q[0];
sx q[0];
rz(3.0101315) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1037285) q[2];
sx q[2];
rz(-1.711039) q[2];
sx q[2];
rz(-0.87272206) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14340885) q[1];
sx q[1];
rz(-1.0034195) q[1];
sx q[1];
rz(-0.25556775) q[1];
x q[2];
rz(-1.6904852) q[3];
sx q[3];
rz(-1.5482248) q[3];
sx q[3];
rz(-2.165497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8160416) q[2];
sx q[2];
rz(-0.52246919) q[2];
sx q[2];
rz(2.6204056) q[2];
rz(0.19736396) q[3];
sx q[3];
rz(-1.3857931) q[3];
sx q[3];
rz(0.50281966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40081438) q[0];
sx q[0];
rz(-0.1939119) q[0];
sx q[0];
rz(0.25522301) q[0];
rz(-0.1420282) q[1];
sx q[1];
rz(-1.7250926) q[1];
sx q[1];
rz(-0.35370383) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58388761) q[0];
sx q[0];
rz(-1.2209799) q[0];
sx q[0];
rz(2.7932998) q[0];
rz(-pi) q[1];
rz(0.92112215) q[2];
sx q[2];
rz(-0.73474681) q[2];
sx q[2];
rz(0.91780797) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4434017) q[1];
sx q[1];
rz(-2.576722) q[1];
sx q[1];
rz(-2.2400212) q[1];
rz(-pi) q[2];
rz(-0.74905101) q[3];
sx q[3];
rz(-1.1370249) q[3];
sx q[3];
rz(-1.449079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3069309) q[2];
sx q[2];
rz(-2.8262409) q[2];
sx q[2];
rz(1.6174512) q[2];
rz(-2.2751685) q[3];
sx q[3];
rz(-1.6774991) q[3];
sx q[3];
rz(2.5518937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50188142) q[0];
sx q[0];
rz(-0.15637936) q[0];
sx q[0];
rz(-1.7891275) q[0];
rz(-1.2347219) q[1];
sx q[1];
rz(-0.23209485) q[1];
sx q[1];
rz(-0.51188767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3054661) q[0];
sx q[0];
rz(-0.16634596) q[0];
sx q[0];
rz(1.2438828) q[0];
rz(-pi) q[1];
rz(1.7281248) q[2];
sx q[2];
rz(-0.81153389) q[2];
sx q[2];
rz(3.1026156) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8066387) q[1];
sx q[1];
rz(-2.6158164) q[1];
sx q[1];
rz(-2.5277565) q[1];
rz(-2.1549822) q[3];
sx q[3];
rz(-1.2704777) q[3];
sx q[3];
rz(-2.9943313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1341683) q[2];
sx q[2];
rz(-1.3205426) q[2];
sx q[2];
rz(0.39719886) q[2];
rz(1.0154137) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(-0.5526244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39887244) q[0];
sx q[0];
rz(-1.9341368) q[0];
sx q[0];
rz(-0.43750986) q[0];
rz(-2.4028026) q[1];
sx q[1];
rz(-0.80016017) q[1];
sx q[1];
rz(-1.4791666) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6300059) q[0];
sx q[0];
rz(-1.5151205) q[0];
sx q[0];
rz(-2.6540709) q[0];
rz(-1.2703367) q[2];
sx q[2];
rz(-2.766937) q[2];
sx q[2];
rz(-0.38196358) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4829172) q[1];
sx q[1];
rz(-1.5086996) q[1];
sx q[1];
rz(0.10644043) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71120925) q[3];
sx q[3];
rz(-2.2235302) q[3];
sx q[3];
rz(-2.0026596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2207569) q[2];
sx q[2];
rz(-1.3584542) q[2];
sx q[2];
rz(-2.9662568) q[2];
rz(1.0026503) q[3];
sx q[3];
rz(-0.23313871) q[3];
sx q[3];
rz(-1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6709082) q[0];
sx q[0];
rz(-1.4574454) q[0];
sx q[0];
rz(-1.1048143) q[0];
rz(0.15057527) q[1];
sx q[1];
rz(-0.87599788) q[1];
sx q[1];
rz(-1.8465975) q[1];
rz(-1.5631093) q[2];
sx q[2];
rz(-1.769071) q[2];
sx q[2];
rz(-1.789112) q[2];
rz(-0.2507052) q[3];
sx q[3];
rz(-1.7964994) q[3];
sx q[3];
rz(-2.9356706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
