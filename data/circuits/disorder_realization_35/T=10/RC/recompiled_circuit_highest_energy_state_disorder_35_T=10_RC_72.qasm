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
rz(-2.45911) q[0];
sx q[0];
rz(-0.64829666) q[0];
sx q[0];
rz(0.034870235) q[0];
rz(1.3999445) q[1];
sx q[1];
rz(-0.26672426) q[1];
sx q[1];
rz(-1.9884225) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0277676) q[0];
sx q[0];
rz(-1.5875289) q[0];
sx q[0];
rz(-1.9219148) q[0];
x q[1];
rz(-1.3787525) q[2];
sx q[2];
rz(-1.3695326) q[2];
sx q[2];
rz(-1.8248688) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.68426484) q[1];
sx q[1];
rz(-2.4870291) q[1];
sx q[1];
rz(1.2689204) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0216924) q[3];
sx q[3];
rz(-1.557805) q[3];
sx q[3];
rz(-1.338122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.030688588) q[2];
sx q[2];
rz(-1.2653753) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.2577308) q[0];
sx q[0];
rz(-0.95153874) q[0];
sx q[0];
rz(1.1269493) q[0];
rz(-1.9583215) q[1];
sx q[1];
rz(-1.1401581) q[1];
sx q[1];
rz(-0.3068876) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81079262) q[0];
sx q[0];
rz(-1.188827) q[0];
sx q[0];
rz(-1.522014) q[0];
rz(-pi) q[1];
rz(-0.82123001) q[2];
sx q[2];
rz(-2.3013287) q[2];
sx q[2];
rz(0.3187547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3908071) q[1];
sx q[1];
rz(-1.7089426) q[1];
sx q[1];
rz(-2.9259794) q[1];
rz(-1.7762228) q[3];
sx q[3];
rz(-2.6449124) q[3];
sx q[3];
rz(1.6261418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7221476) q[2];
sx q[2];
rz(-0.38997969) q[2];
sx q[2];
rz(1.6411714) q[2];
rz(1.8309343) q[3];
sx q[3];
rz(-2.1597517) q[3];
sx q[3];
rz(-0.76154861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(0.2181181) q[1];
sx q[1];
rz(-1.7354542) q[1];
sx q[1];
rz(1.1880147) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145633) q[0];
sx q[0];
rz(-0.62262056) q[0];
sx q[0];
rz(-0.027986469) q[0];
rz(-pi) q[1];
rz(1.4039878) q[2];
sx q[2];
rz(-2.1412289) q[2];
sx q[2];
rz(-0.32329971) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0289171) q[1];
sx q[1];
rz(-0.099359902) q[1];
sx q[1];
rz(3.1032853) q[1];
rz(-3.1165666) q[3];
sx q[3];
rz(-0.2215759) q[3];
sx q[3];
rz(3.0794249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5334566) q[2];
sx q[2];
rz(-0.44508219) q[2];
sx q[2];
rz(-0.76370561) q[2];
rz(-2.8869827) q[3];
sx q[3];
rz(-2.2201241) q[3];
sx q[3];
rz(0.34847954) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3449645) q[0];
sx q[0];
rz(-0.94709713) q[0];
sx q[0];
rz(-3.1297041) q[0];
rz(-2.1628926) q[1];
sx q[1];
rz(-1.3589166) q[1];
sx q[1];
rz(1.1241283) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91299401) q[0];
sx q[0];
rz(-1.9833195) q[0];
sx q[0];
rz(-1.8275632) q[0];
x q[1];
rz(-2.5619036) q[2];
sx q[2];
rz(-1.9429029) q[2];
sx q[2];
rz(1.6339782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1710137) q[1];
sx q[1];
rz(-2.5742324) q[1];
sx q[1];
rz(2.2544276) q[1];
rz(-pi) q[2];
rz(-0.57321623) q[3];
sx q[3];
rz(-2.3251136) q[3];
sx q[3];
rz(3.1117718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0594788) q[2];
sx q[2];
rz(-0.79136005) q[2];
sx q[2];
rz(0.7938844) q[2];
rz(-1.7700899) q[3];
sx q[3];
rz(-0.81727782) q[3];
sx q[3];
rz(-2.1577458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3802948) q[0];
sx q[0];
rz(-1.5341606) q[0];
sx q[0];
rz(0.5436596) q[0];
rz(-1.9937953) q[1];
sx q[1];
rz(-0.59998435) q[1];
sx q[1];
rz(2.0701087) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5872566) q[0];
sx q[0];
rz(-1.2859467) q[0];
sx q[0];
rz(1.1469123) q[0];
x q[1];
rz(-1.7254099) q[2];
sx q[2];
rz(-1.8616759) q[2];
sx q[2];
rz(2.9335528) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0669015) q[1];
sx q[1];
rz(-1.2558535) q[1];
sx q[1];
rz(2.6661162) q[1];
rz(-pi) q[2];
rz(-0.50258639) q[3];
sx q[3];
rz(-2.2690176) q[3];
sx q[3];
rz(1.0909441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.063283) q[2];
sx q[2];
rz(-0.50476837) q[2];
sx q[2];
rz(2.2341991) q[2];
rz(0.0026957938) q[3];
sx q[3];
rz(-1.4009652) q[3];
sx q[3];
rz(1.9122972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9232848) q[0];
sx q[0];
rz(-1.2907499) q[0];
sx q[0];
rz(1.9521335) q[0];
rz(-0.16667357) q[1];
sx q[1];
rz(-1.0308301) q[1];
sx q[1];
rz(1.6698242) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4772861) q[0];
sx q[0];
rz(-0.95726958) q[0];
sx q[0];
rz(1.0856347) q[0];
rz(-1.6138863) q[2];
sx q[2];
rz(-2.913007) q[2];
sx q[2];
rz(3.0800386) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1907726) q[1];
sx q[1];
rz(-0.69389486) q[1];
sx q[1];
rz(2.4681952) q[1];
x q[2];
rz(2.9091524) q[3];
sx q[3];
rz(-0.83139426) q[3];
sx q[3];
rz(-1.1375994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2344096) q[2];
sx q[2];
rz(-0.11203354) q[2];
sx q[2];
rz(0.33667931) q[2];
rz(1.178721) q[3];
sx q[3];
rz(-1.7095292) q[3];
sx q[3];
rz(2.7746576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1692899) q[0];
sx q[0];
rz(-0.65880913) q[0];
sx q[0];
rz(1.4038053) q[0];
rz(-2.9804969) q[1];
sx q[1];
rz(-0.95268551) q[1];
sx q[1];
rz(-0.62532997) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1243131) q[0];
sx q[0];
rz(-1.0758022) q[0];
sx q[0];
rz(0.87484579) q[0];
rz(-2.2032878) q[2];
sx q[2];
rz(-1.9646346) q[2];
sx q[2];
rz(2.1932067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8323601) q[1];
sx q[1];
rz(-0.52150583) q[1];
sx q[1];
rz(1.4307561) q[1];
rz(0.13631374) q[3];
sx q[3];
rz(-1.1778616) q[3];
sx q[3];
rz(-0.64717347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8620341) q[2];
sx q[2];
rz(-1.4906078) q[2];
sx q[2];
rz(3.1107235) q[2];
rz(-0.75012642) q[3];
sx q[3];
rz(-2.153219) q[3];
sx q[3];
rz(-1.8376902) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6136762) q[0];
sx q[0];
rz(-0.91687098) q[0];
sx q[0];
rz(3.0314639) q[0];
rz(2.2573709) q[1];
sx q[1];
rz(-1.607211) q[1];
sx q[1];
rz(-1.0482739) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78335242) q[0];
sx q[0];
rz(-1.269258) q[0];
sx q[0];
rz(-1.8389971) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8364622) q[2];
sx q[2];
rz(-1.507826) q[2];
sx q[2];
rz(2.8612325) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.058052651) q[1];
sx q[1];
rz(-1.3370635) q[1];
sx q[1];
rz(-1.5104945) q[1];
x q[2];
rz(-3.0804668) q[3];
sx q[3];
rz(-1.3406274) q[3];
sx q[3];
rz(-1.597061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5732164) q[2];
sx q[2];
rz(-1.6171075) q[2];
sx q[2];
rz(-0.035408346) q[2];
rz(0.046444166) q[3];
sx q[3];
rz(-1.2848264) q[3];
sx q[3];
rz(-3.0845643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7923918) q[0];
sx q[0];
rz(-0.89850315) q[0];
sx q[0];
rz(-3.0080556) q[0];
rz(1.4117433) q[1];
sx q[1];
rz(-1.1213419) q[1];
sx q[1];
rz(-1.0498566) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8916666) q[0];
sx q[0];
rz(-1.6841411) q[0];
sx q[0];
rz(0.38007315) q[0];
rz(-pi) q[1];
x q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7375917) q[1];
sx q[1];
rz(-1.5301413) q[1];
sx q[1];
rz(2.3126751) q[1];
rz(-pi) q[2];
rz(2.7383835) q[3];
sx q[3];
rz(-1.8397188) q[3];
sx q[3];
rz(1.2300242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.92020804) q[2];
sx q[2];
rz(-0.68156663) q[2];
sx q[2];
rz(-1.3725012) q[2];
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
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16804279) q[0];
sx q[0];
rz(-1.9815227) q[0];
sx q[0];
rz(-2.2450182) q[0];
rz(0.98094455) q[1];
sx q[1];
rz(-2.4368821) q[1];
sx q[1];
rz(-3.0006345) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94146848) q[0];
sx q[0];
rz(-0.57766047) q[0];
sx q[0];
rz(2.0117652) q[0];
rz(-pi) q[1];
rz(0.015816255) q[2];
sx q[2];
rz(-2.2770303) q[2];
sx q[2];
rz(-0.13705239) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.28480865) q[1];
sx q[1];
rz(-1.5601381) q[1];
sx q[1];
rz(-2.6066236) q[1];
rz(2.6871164) q[3];
sx q[3];
rz(-2.6095172) q[3];
sx q[3];
rz(-0.34031351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.36833802) q[2];
sx q[2];
rz(-2.1835623) q[2];
sx q[2];
rz(-2.9216596) q[2];
rz(2.1851165) q[3];
sx q[3];
rz(-0.82436162) q[3];
sx q[3];
rz(-2.1970356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3750951) q[0];
sx q[0];
rz(-1.5506333) q[0];
sx q[0];
rz(-0.63313708) q[0];
rz(1.7924869) q[1];
sx q[1];
rz(-0.92020412) q[1];
sx q[1];
rz(-0.79457582) q[1];
rz(2.6146087) q[2];
sx q[2];
rz(-2.5415642) q[2];
sx q[2];
rz(2.1733282) q[2];
rz(1.7344162) q[3];
sx q[3];
rz(-2.4118578) q[3];
sx q[3];
rz(2.2414997) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
