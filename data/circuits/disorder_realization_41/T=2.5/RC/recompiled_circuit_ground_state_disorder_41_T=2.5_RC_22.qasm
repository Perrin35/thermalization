OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95076686) q[0];
sx q[0];
rz(-1.8523676) q[0];
sx q[0];
rz(-3.0486795) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(-0.73268259) q[1];
sx q[1];
rz(-1.2394152) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3071741) q[0];
sx q[0];
rz(-2.7980013) q[0];
sx q[0];
rz(0.8309721) q[0];
rz(-pi) q[1];
rz(-1.7669414) q[2];
sx q[2];
rz(-1.3117547) q[2];
sx q[2];
rz(2.1262622) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.69529205) q[1];
sx q[1];
rz(-1.2406893) q[1];
sx q[1];
rz(1.6020726) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0650571) q[3];
sx q[3];
rz(-2.0856033) q[3];
sx q[3];
rz(2.6635287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6713082) q[2];
sx q[2];
rz(-1.4802063) q[2];
sx q[2];
rz(1.7926463) q[2];
rz(0.74364439) q[3];
sx q[3];
rz(-2.8642004) q[3];
sx q[3];
rz(-2.9644137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0002366) q[0];
sx q[0];
rz(-1.5445222) q[0];
sx q[0];
rz(2.8479688) q[0];
rz(2.1108421) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(0.34293276) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021504952) q[0];
sx q[0];
rz(-1.3653127) q[0];
sx q[0];
rz(1.0387102) q[0];
rz(-pi) q[1];
rz(1.8937102) q[2];
sx q[2];
rz(-0.39034778) q[2];
sx q[2];
rz(2.5935612) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1520815) q[1];
sx q[1];
rz(-2.40675) q[1];
sx q[1];
rz(-2.8192725) q[1];
x q[2];
rz(0.25094752) q[3];
sx q[3];
rz(-0.95635575) q[3];
sx q[3];
rz(0.80048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0127516) q[2];
sx q[2];
rz(-1.7495456) q[2];
sx q[2];
rz(-0.7592321) q[2];
rz(3.1208842) q[3];
sx q[3];
rz(-1.8755251) q[3];
sx q[3];
rz(-0.43829632) q[3];
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
rz(-0.52221209) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(0.0044599175) q[0];
rz(-2.4941173) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(-0.78537816) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42703907) q[0];
sx q[0];
rz(-1.6511962) q[0];
sx q[0];
rz(-0.15958448) q[0];
rz(-pi) q[1];
rz(-1.8846604) q[2];
sx q[2];
rz(-1.7376657) q[2];
sx q[2];
rz(0.97970574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0177893) q[1];
sx q[1];
rz(-0.54517704) q[1];
sx q[1];
rz(-1.4383609) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.249751) q[3];
sx q[3];
rz(-1.0247318) q[3];
sx q[3];
rz(0.31072703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.37041) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(1.7837589) q[2];
rz(2.8010098) q[3];
sx q[3];
rz(-0.95652306) q[3];
sx q[3];
rz(-0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99712813) q[0];
sx q[0];
rz(-1.0839394) q[0];
sx q[0];
rz(-0.88515627) q[0];
rz(-2.3947233) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(-1.0964099) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089450739) q[0];
sx q[0];
rz(-0.47596395) q[0];
sx q[0];
rz(2.3030998) q[0];
rz(-1.095216) q[2];
sx q[2];
rz(-3.0747483) q[2];
sx q[2];
rz(0.16193709) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4064286) q[1];
sx q[1];
rz(-2.6044835) q[1];
sx q[1];
rz(-2.6178611) q[1];
x q[2];
rz(1.0616333) q[3];
sx q[3];
rz(-2.235755) q[3];
sx q[3];
rz(0.10275118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5235644) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(0.024624126) q[2];
rz(0.63557449) q[3];
sx q[3];
rz(-0.98819757) q[3];
sx q[3];
rz(-2.2583101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6898952) q[0];
sx q[0];
rz(-0.30245936) q[0];
sx q[0];
rz(0.46689335) q[0];
rz(-3.0310071) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(-2.0707524) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8798053) q[0];
sx q[0];
rz(-2.6222485) q[0];
sx q[0];
rz(2.065361) q[0];
rz(1.7532639) q[2];
sx q[2];
rz(-1.8844205) q[2];
sx q[2];
rz(-0.30280606) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.80761284) q[1];
sx q[1];
rz(-0.65704221) q[1];
sx q[1];
rz(0.64589898) q[1];
rz(-pi) q[2];
rz(-2.4399906) q[3];
sx q[3];
rz(-1.0402586) q[3];
sx q[3];
rz(-0.70487937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0232627) q[2];
sx q[2];
rz(-1.7795965) q[2];
sx q[2];
rz(-2.2892717) q[2];
rz(-0.037633745) q[3];
sx q[3];
rz(-1.9084385) q[3];
sx q[3];
rz(2.5467303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.0090050176) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(-1.2888541) q[0];
rz(1.6291078) q[1];
sx q[1];
rz(-1.1237203) q[1];
sx q[1];
rz(1.1423133) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70002551) q[0];
sx q[0];
rz(-1.6684932) q[0];
sx q[0];
rz(-2.7872681) q[0];
rz(-pi) q[1];
rz(-1.8577447) q[2];
sx q[2];
rz(-2.9252238) q[2];
sx q[2];
rz(-1.8277825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2494275) q[1];
sx q[1];
rz(-1.9122423) q[1];
sx q[1];
rz(-1.3072827) q[1];
rz(-pi) q[2];
rz(0.50636815) q[3];
sx q[3];
rz(-0.35249235) q[3];
sx q[3];
rz(-1.6329721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(-3.0211871) q[2];
rz(-1.985792) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(-0.22404484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8026546) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(1.6876203) q[0];
rz(1.5441719) q[1];
sx q[1];
rz(-1.495196) q[1];
sx q[1];
rz(2.6905751) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16674834) q[0];
sx q[0];
rz(-1.5112425) q[0];
sx q[0];
rz(-0.3015428) q[0];
rz(0.23239179) q[2];
sx q[2];
rz(-1.2133779) q[2];
sx q[2];
rz(-3.0136303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6587131) q[1];
sx q[1];
rz(-0.45156839) q[1];
sx q[1];
rz(2.8612479) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5805832) q[3];
sx q[3];
rz(-0.76044816) q[3];
sx q[3];
rz(-1.994054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.94577998) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(-1.3607402) q[2];
rz(-2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(1.7520693) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1239531) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(-0.70607591) q[0];
rz(1.5977244) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(-2.2854038) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34205758) q[0];
sx q[0];
rz(-1.5634131) q[0];
sx q[0];
rz(1.5513585) q[0];
rz(-pi) q[1];
rz(-2.0308308) q[2];
sx q[2];
rz(-0.80321124) q[2];
sx q[2];
rz(0.75424657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1284898) q[1];
sx q[1];
rz(-1.0016514) q[1];
sx q[1];
rz(1.8123958) q[1];
x q[2];
rz(-2.059465) q[3];
sx q[3];
rz(-1.6889945) q[3];
sx q[3];
rz(-1.4335872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2989444) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(0.16743463) q[2];
rz(1.5542479) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(-2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7635968) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(-0.3717306) q[0];
rz(1.4338214) q[1];
sx q[1];
rz(-1.6038409) q[1];
sx q[1];
rz(-2.8415714) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2903295) q[0];
sx q[0];
rz(-1.6870878) q[0];
sx q[0];
rz(2.854748) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57711592) q[2];
sx q[2];
rz(-1.6047928) q[2];
sx q[2];
rz(0.94552065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3206818) q[1];
sx q[1];
rz(-1.6947985) q[1];
sx q[1];
rz(-2.4895666) q[1];
x q[2];
rz(0.12436538) q[3];
sx q[3];
rz(-0.87166407) q[3];
sx q[3];
rz(-1.82064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6649449) q[2];
sx q[2];
rz(-2.6740394) q[2];
sx q[2];
rz(-2.7139968) q[2];
rz(-1.556373) q[3];
sx q[3];
rz(-1.1170324) q[3];
sx q[3];
rz(0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72117358) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(-0.46947259) q[0];
rz(0.92533127) q[1];
sx q[1];
rz(-1.0877437) q[1];
sx q[1];
rz(0.023177711) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0231004) q[0];
sx q[0];
rz(-1.4064624) q[0];
sx q[0];
rz(-2.611155) q[0];
x q[1];
rz(-0.9018578) q[2];
sx q[2];
rz(-2.0624071) q[2];
sx q[2];
rz(0.41504809) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7851023) q[1];
sx q[1];
rz(-1.7447326) q[1];
sx q[1];
rz(2.9746303) q[1];
x q[2];
rz(-2.9523207) q[3];
sx q[3];
rz(-2.2918275) q[3];
sx q[3];
rz(1.7060351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2716486) q[2];
sx q[2];
rz(-1.2056489) q[2];
sx q[2];
rz(-0.49087697) q[2];
rz(1.0902181) q[3];
sx q[3];
rz(-0.96929437) q[3];
sx q[3];
rz(-2.8276665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8563817) q[0];
sx q[0];
rz(-2.2638392) q[0];
sx q[0];
rz(1.0765156) q[0];
rz(-0.27676997) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(-1.8770915) q[2];
sx q[2];
rz(-2.2363792) q[2];
sx q[2];
rz(-2.1909942) q[2];
rz(1.5908949) q[3];
sx q[3];
rz(-1.0363967) q[3];
sx q[3];
rz(-3.1283621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
