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
rz(0.68248263) q[0];
sx q[0];
rz(-2.493296) q[0];
sx q[0];
rz(-0.034870235) q[0];
rz(1.3999445) q[1];
sx q[1];
rz(6.016461) q[1];
sx q[1];
rz(7.4363554) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54915743) q[0];
sx q[0];
rz(-1.2197291) q[0];
sx q[0];
rz(0.017819509) q[0];
rz(-1.7628402) q[2];
sx q[2];
rz(-1.7720601) q[2];
sx q[2];
rz(1.3167238) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64434796) q[1];
sx q[1];
rz(-1.3887857) q[1];
sx q[1];
rz(-2.2031511) q[1];
rz(-pi) q[2];
rz(0.014433666) q[3];
sx q[3];
rz(-1.1199411) q[3];
sx q[3];
rz(-2.9026287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1109041) q[2];
sx q[2];
rz(-1.8762174) q[2];
sx q[2];
rz(0.14070025) q[2];
rz(-1.5594907) q[3];
sx q[3];
rz(-2.9086106) q[3];
sx q[3];
rz(-1.8039907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2577308) q[0];
sx q[0];
rz(-0.95153874) q[0];
sx q[0];
rz(-2.0146433) q[0];
rz(1.1832712) q[1];
sx q[1];
rz(-1.1401581) q[1];
sx q[1];
rz(2.8347051) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94102717) q[0];
sx q[0];
rz(-0.38492003) q[0];
sx q[0];
rz(-0.12080111) q[0];
rz(-pi) q[1];
rz(0.65019239) q[2];
sx q[2];
rz(-0.99417328) q[2];
sx q[2];
rz(1.2675831) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3813016) q[1];
sx q[1];
rz(-0.25549251) q[1];
sx q[1];
rz(-0.5762655) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3653698) q[3];
sx q[3];
rz(-2.6449124) q[3];
sx q[3];
rz(-1.5154509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41944501) q[2];
sx q[2];
rz(-2.751613) q[2];
sx q[2];
rz(-1.5004213) q[2];
rz(1.8309343) q[3];
sx q[3];
rz(-2.1597517) q[3];
sx q[3];
rz(2.380044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0618884) q[0];
sx q[0];
rz(-0.56489983) q[0];
sx q[0];
rz(-1.0796219) q[0];
rz(-2.9234746) q[1];
sx q[1];
rz(-1.7354542) q[1];
sx q[1];
rz(1.1880147) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2750888) q[0];
sx q[0];
rz(-1.554477) q[0];
sx q[0];
rz(2.5191576) q[0];
rz(-0.57680371) q[2];
sx q[2];
rz(-1.7110023) q[2];
sx q[2];
rz(1.9847676) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0674141) q[1];
sx q[1];
rz(-1.6700831) q[1];
sx q[1];
rz(1.5746142) q[1];
x q[2];
rz(-0.22150877) q[3];
sx q[3];
rz(-1.5762957) q[3];
sx q[3];
rz(1.4842141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5334566) q[2];
sx q[2];
rz(-0.44508219) q[2];
sx q[2];
rz(-0.76370561) q[2];
rz(-2.8869827) q[3];
sx q[3];
rz(-0.92146858) q[3];
sx q[3];
rz(2.7931131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3449645) q[0];
sx q[0];
rz(-0.94709713) q[0];
sx q[0];
rz(3.1297041) q[0];
rz(-2.1628926) q[1];
sx q[1];
rz(-1.782676) q[1];
sx q[1];
rz(-1.1241283) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55292144) q[0];
sx q[0];
rz(-1.8055989) q[0];
sx q[0];
rz(0.42490439) q[0];
rz(-pi) q[1];
rz(-2.0072859) q[2];
sx q[2];
rz(-1.0352897) q[2];
sx q[2];
rz(2.8447163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.939399) q[1];
sx q[1];
rz(-1.2244818) q[1];
sx q[1];
rz(-1.1119197) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5683764) q[3];
sx q[3];
rz(-2.3251136) q[3];
sx q[3];
rz(-0.029820853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0594788) q[2];
sx q[2];
rz(-0.79136005) q[2];
sx q[2];
rz(0.7938844) q[2];
rz(-1.7700899) q[3];
sx q[3];
rz(-2.3243148) q[3];
sx q[3];
rz(-0.98384682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7612979) q[0];
sx q[0];
rz(-1.5341606) q[0];
sx q[0];
rz(2.5979331) q[0];
rz(1.1477973) q[1];
sx q[1];
rz(-2.5416083) q[1];
sx q[1];
rz(-2.0701087) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5543361) q[0];
sx q[0];
rz(-1.8556459) q[0];
sx q[0];
rz(1.9946804) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29419326) q[2];
sx q[2];
rz(-1.718866) q[2];
sx q[2];
rz(1.8235056) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0746911) q[1];
sx q[1];
rz(-1.2558535) q[1];
sx q[1];
rz(-0.47547646) q[1];
rz(-pi) q[2];
rz(2.6390063) q[3];
sx q[3];
rz(-0.87257507) q[3];
sx q[3];
rz(-1.0909441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.063283) q[2];
sx q[2];
rz(-2.6368243) q[2];
sx q[2];
rz(-2.2341991) q[2];
rz(-0.0026957938) q[3];
sx q[3];
rz(-1.4009652) q[3];
sx q[3];
rz(1.2292954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9232848) q[0];
sx q[0];
rz(-1.2907499) q[0];
sx q[0];
rz(1.9521335) q[0];
rz(0.16667357) q[1];
sx q[1];
rz(-2.1107626) q[1];
sx q[1];
rz(1.6698242) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6643065) q[0];
sx q[0];
rz(-2.1843231) q[0];
sx q[0];
rz(-1.0856347) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.010021539) q[2];
sx q[2];
rz(-1.3424266) q[2];
sx q[2];
rz(-0.10579338) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7547524) q[1];
sx q[1];
rz(-2.0943145) q[1];
sx q[1];
rz(2.0493815) q[1];
x q[2];
rz(1.8182033) q[3];
sx q[3];
rz(-0.76843211) q[3];
sx q[3];
rz(-0.79977126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2344096) q[2];
sx q[2];
rz(-3.0295591) q[2];
sx q[2];
rz(-0.33667931) q[2];
rz(1.178721) q[3];
sx q[3];
rz(-1.7095292) q[3];
sx q[3];
rz(2.7746576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1692899) q[0];
sx q[0];
rz(-2.4827835) q[0];
sx q[0];
rz(-1.4038053) q[0];
rz(-0.16109578) q[1];
sx q[1];
rz(-2.1889071) q[1];
sx q[1];
rz(-0.62532997) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0172796) q[0];
sx q[0];
rz(-1.0758022) q[0];
sx q[0];
rz(-0.87484579) q[0];
rz(0.9383049) q[2];
sx q[2];
rz(-1.9646346) q[2];
sx q[2];
rz(2.1932067) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.13994234) q[1];
sx q[1];
rz(-1.6403908) q[1];
sx q[1];
rz(-1.0535295) q[1];
rz(-pi) q[2];
rz(1.9670385) q[3];
sx q[3];
rz(-1.696664) q[3];
sx q[3];
rz(-2.270442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.27955851) q[2];
sx q[2];
rz(-1.6509849) q[2];
sx q[2];
rz(0.030869182) q[2];
rz(-0.75012642) q[3];
sx q[3];
rz(-0.98837367) q[3];
sx q[3];
rz(-1.3039024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6136762) q[0];
sx q[0];
rz(-2.2247217) q[0];
sx q[0];
rz(3.0314639) q[0];
rz(-0.88422173) q[1];
sx q[1];
rz(-1.607211) q[1];
sx q[1];
rz(2.0933188) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78335242) q[0];
sx q[0];
rz(-1.8723346) q[0];
sx q[0];
rz(-1.3025955) q[0];
rz(-pi) q[1];
rz(-3.0763392) q[2];
sx q[2];
rz(-1.3056697) q[2];
sx q[2];
rz(-1.2733151) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.058052651) q[1];
sx q[1];
rz(-1.3370635) q[1];
sx q[1];
rz(-1.5104945) q[1];
x q[2];
rz(-0.061125867) q[3];
sx q[3];
rz(-1.8009652) q[3];
sx q[3];
rz(1.5445316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5732164) q[2];
sx q[2];
rz(-1.6171075) q[2];
sx q[2];
rz(-3.1061843) q[2];
rz(-3.0951485) q[3];
sx q[3];
rz(-1.2848264) q[3];
sx q[3];
rz(-3.0845643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34920084) q[0];
sx q[0];
rz(-2.2430895) q[0];
sx q[0];
rz(3.0080556) q[0];
rz(-1.4117433) q[1];
sx q[1];
rz(-1.1213419) q[1];
sx q[1];
rz(-2.0917361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8658757) q[0];
sx q[0];
rz(-1.1932853) q[0];
sx q[0];
rz(-1.4488245) q[0];
rz(-pi) q[1];
rz(2.1112415) q[2];
sx q[2];
rz(-1.5771416) q[2];
sx q[2];
rz(-0.66101569) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9304815) q[1];
sx q[1];
rz(-0.74277987) q[1];
sx q[1];
rz(-1.5106661) q[1];
x q[2];
rz(2.7383835) q[3];
sx q[3];
rz(-1.8397188) q[3];
sx q[3];
rz(-1.9115684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.92020804) q[2];
sx q[2];
rz(-0.68156663) q[2];
sx q[2];
rz(1.3725012) q[2];
rz(1.3051322) q[3];
sx q[3];
rz(-1.2499481) q[3];
sx q[3];
rz(2.5398152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16804279) q[0];
sx q[0];
rz(-1.9815227) q[0];
sx q[0];
rz(0.89657441) q[0];
rz(-0.98094455) q[1];
sx q[1];
rz(-0.70471057) q[1];
sx q[1];
rz(-3.0006345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6870689) q[0];
sx q[0];
rz(-2.087283) q[0];
sx q[0];
rz(0.27134925) q[0];
rz(-pi) q[1];
rz(-3.1257764) q[2];
sx q[2];
rz(-0.86456233) q[2];
sx q[2];
rz(0.13705239) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.856784) q[1];
sx q[1];
rz(-1.5601381) q[1];
sx q[1];
rz(2.6066236) q[1];
x q[2];
rz(-2.6550547) q[3];
sx q[3];
rz(-1.7953904) q[3];
sx q[3];
rz(2.3096245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7732546) q[2];
sx q[2];
rz(-2.1835623) q[2];
sx q[2];
rz(-0.21993302) q[2];
rz(0.95647612) q[3];
sx q[3];
rz(-2.317231) q[3];
sx q[3];
rz(0.94455702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7664976) q[0];
sx q[0];
rz(-1.5909593) q[0];
sx q[0];
rz(2.5084556) q[0];
rz(-1.3491058) q[1];
sx q[1];
rz(-0.92020412) q[1];
sx q[1];
rz(-0.79457582) q[1];
rz(2.6146087) q[2];
sx q[2];
rz(-2.5415642) q[2];
sx q[2];
rz(2.1733282) q[2];
rz(2.293855) q[3];
sx q[3];
rz(-1.6796057) q[3];
sx q[3];
rz(-2.3484505) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
