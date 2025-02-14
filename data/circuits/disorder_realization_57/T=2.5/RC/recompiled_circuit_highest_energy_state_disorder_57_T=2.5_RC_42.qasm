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
rz(1.3996742) q[0];
sx q[0];
rz(-1.6784486) q[0];
sx q[0];
rz(-1.4611257) q[0];
rz(2.6598016) q[1];
sx q[1];
rz(-2.3269589) q[1];
sx q[1];
rz(2.2021267) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6326673) q[0];
sx q[0];
rz(-1.5611751) q[0];
sx q[0];
rz(1.5271355) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5599776) q[2];
sx q[2];
rz(-2.5293859) q[2];
sx q[2];
rz(-0.78919995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2880683) q[1];
sx q[1];
rz(-1.6109835) q[1];
sx q[1];
rz(0.010577135) q[1];
x q[2];
rz(1.3529196) q[3];
sx q[3];
rz(-1.1718201) q[3];
sx q[3];
rz(-0.50673317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7751094) q[2];
sx q[2];
rz(-1.6747687) q[2];
sx q[2];
rz(-0.62466204) q[2];
rz(-2.1810253) q[3];
sx q[3];
rz(-1.4550236) q[3];
sx q[3];
rz(-0.87489405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4082773) q[0];
sx q[0];
rz(-0.69077078) q[0];
sx q[0];
rz(-2.9498003) q[0];
rz(-2.5674112) q[1];
sx q[1];
rz(-1.6185113) q[1];
sx q[1];
rz(-0.40723732) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2977266) q[0];
sx q[0];
rz(-0.17619421) q[0];
sx q[0];
rz(-1.5017444) q[0];
x q[1];
rz(0.064795061) q[2];
sx q[2];
rz(-1.7470834) q[2];
sx q[2];
rz(0.48688146) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6169238) q[1];
sx q[1];
rz(-1.7546931) q[1];
sx q[1];
rz(-2.6016629) q[1];
rz(3.0519227) q[3];
sx q[3];
rz(-2.0715791) q[3];
sx q[3];
rz(1.9865685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9521956) q[2];
sx q[2];
rz(-1.1789097) q[2];
sx q[2];
rz(-0.65323812) q[2];
rz(-2.2630528) q[3];
sx q[3];
rz(-1.0983175) q[3];
sx q[3];
rz(-3.0285335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1427796) q[0];
sx q[0];
rz(-1.7072562) q[0];
sx q[0];
rz(1.0502195) q[0];
rz(-0.6160008) q[1];
sx q[1];
rz(-1.7237677) q[1];
sx q[1];
rz(0.91316694) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8693059) q[0];
sx q[0];
rz(-1.7401631) q[0];
sx q[0];
rz(-3.0336551) q[0];
x q[1];
rz(0.11868015) q[2];
sx q[2];
rz(-1.2859697) q[2];
sx q[2];
rz(2.8800346) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2985067) q[1];
sx q[1];
rz(-1.8131327) q[1];
sx q[1];
rz(1.0964811) q[1];
rz(-pi) q[2];
rz(-0.3766587) q[3];
sx q[3];
rz(-0.91729255) q[3];
sx q[3];
rz(2.3953826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.22127613) q[2];
sx q[2];
rz(-2.2344799) q[2];
sx q[2];
rz(3.082412) q[2];
rz(-0.83909285) q[3];
sx q[3];
rz(-1.9234761) q[3];
sx q[3];
rz(0.86735094) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7729823) q[0];
sx q[0];
rz(-1.0091857) q[0];
sx q[0];
rz(0.34977812) q[0];
rz(1.7971669) q[1];
sx q[1];
rz(-2.5549922) q[1];
sx q[1];
rz(-0.1127359) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37575133) q[0];
sx q[0];
rz(-2.021702) q[0];
sx q[0];
rz(-1.1896053) q[0];
rz(0.91103986) q[2];
sx q[2];
rz(-1.3415114) q[2];
sx q[2];
rz(-0.89067543) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60242059) q[1];
sx q[1];
rz(-0.93257553) q[1];
sx q[1];
rz(-2.763487) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2318683) q[3];
sx q[3];
rz(-2.1436757) q[3];
sx q[3];
rz(-2.6343144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.747644) q[2];
sx q[2];
rz(-1.5564206) q[2];
sx q[2];
rz(-0.057223884) q[2];
rz(-1.8852437) q[3];
sx q[3];
rz(-2.5689954) q[3];
sx q[3];
rz(-2.4109139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8771862) q[0];
sx q[0];
rz(-1.1671966) q[0];
sx q[0];
rz(2.4315244) q[0];
rz(0.26142985) q[1];
sx q[1];
rz(-1.5855887) q[1];
sx q[1];
rz(0.95103055) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3426845) q[0];
sx q[0];
rz(-2.5747188) q[0];
sx q[0];
rz(-0.64115144) q[0];
x q[1];
rz(-2.0567953) q[2];
sx q[2];
rz(-2.1446781) q[2];
sx q[2];
rz(1.6465555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3465865) q[1];
sx q[1];
rz(-1.7730496) q[1];
sx q[1];
rz(0.70960416) q[1];
x q[2];
rz(-1.4360436) q[3];
sx q[3];
rz(-0.6373848) q[3];
sx q[3];
rz(-2.4911043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6422673) q[2];
sx q[2];
rz(-0.32411164) q[2];
sx q[2];
rz(2.8322241) q[2];
rz(2.046509) q[3];
sx q[3];
rz(-2.0131922) q[3];
sx q[3];
rz(-0.40422082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22769044) q[0];
sx q[0];
rz(-2.9588283) q[0];
sx q[0];
rz(-0.90232724) q[0];
rz(0.062571049) q[1];
sx q[1];
rz(-1.2689271) q[1];
sx q[1];
rz(1.2455469) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7880873) q[0];
sx q[0];
rz(-2.2416032) q[0];
sx q[0];
rz(1.3735347) q[0];
x q[1];
rz(-1.0311216) q[2];
sx q[2];
rz(-1.5226242) q[2];
sx q[2];
rz(-1.1205412) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77038232) q[1];
sx q[1];
rz(-2.4744445) q[1];
sx q[1];
rz(-2.4333111) q[1];
x q[2];
rz(0.67288237) q[3];
sx q[3];
rz(-1.1816634) q[3];
sx q[3];
rz(-1.0997888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7237599) q[2];
sx q[2];
rz(-1.2699026) q[2];
sx q[2];
rz(1.7977686) q[2];
rz(1.2835245) q[3];
sx q[3];
rz(-0.44461861) q[3];
sx q[3];
rz(-2.4166687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.035324) q[0];
sx q[0];
rz(-3.1043053) q[0];
sx q[0];
rz(2.7653747) q[0];
rz(0.47209921) q[1];
sx q[1];
rz(-0.68196982) q[1];
sx q[1];
rz(0.45904407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7721586) q[0];
sx q[0];
rz(-0.87352814) q[0];
sx q[0];
rz(-0.56241547) q[0];
rz(-pi) q[1];
rz(2.6920175) q[2];
sx q[2];
rz(-1.4130985) q[2];
sx q[2];
rz(-1.1166935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9320786) q[1];
sx q[1];
rz(-2.9649295) q[1];
sx q[1];
rz(-2.7448369) q[1];
rz(-pi) q[2];
rz(-0.59765069) q[3];
sx q[3];
rz(-1.581218) q[3];
sx q[3];
rz(2.8437993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0692733) q[2];
sx q[2];
rz(-0.71162144) q[2];
sx q[2];
rz(2.7834328) q[2];
rz(0.88851309) q[3];
sx q[3];
rz(-2.1202959) q[3];
sx q[3];
rz(-2.2108938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096980378) q[0];
sx q[0];
rz(-2.8428069) q[0];
sx q[0];
rz(0.97691798) q[0];
rz(-0.75374976) q[1];
sx q[1];
rz(-1.6981373) q[1];
sx q[1];
rz(-1.5026198) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.174265) q[0];
sx q[0];
rz(-0.87192029) q[0];
sx q[0];
rz(-1.4682659) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5112652) q[2];
sx q[2];
rz(-2.3663372) q[2];
sx q[2];
rz(1.1182736) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21587196) q[1];
sx q[1];
rz(-1.873107) q[1];
sx q[1];
rz(0.69958256) q[1];
rz(-pi) q[2];
rz(0.23076337) q[3];
sx q[3];
rz(-2.181859) q[3];
sx q[3];
rz(-2.3045674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94540191) q[2];
sx q[2];
rz(-1.4620898) q[2];
sx q[2];
rz(0.84575829) q[2];
rz(2.3614007) q[3];
sx q[3];
rz(-1.3581685) q[3];
sx q[3];
rz(-2.5030524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82275259) q[0];
sx q[0];
rz(-1.0718811) q[0];
sx q[0];
rz(2.7986797) q[0];
rz(-2.137939) q[1];
sx q[1];
rz(-1.2316848) q[1];
sx q[1];
rz(2.4249605) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4591689) q[0];
sx q[0];
rz(-0.83324106) q[0];
sx q[0];
rz(0.47969819) q[0];
x q[1];
rz(-0.54729531) q[2];
sx q[2];
rz(-2.5792053) q[2];
sx q[2];
rz(-2.9862816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3874599) q[1];
sx q[1];
rz(-2.4841822) q[1];
sx q[1];
rz(-0.44771835) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4695806) q[3];
sx q[3];
rz(-1.8895961) q[3];
sx q[3];
rz(2.4254951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0922682) q[2];
sx q[2];
rz(-1.5169787) q[2];
sx q[2];
rz(-2.1410904) q[2];
rz(2.3246824) q[3];
sx q[3];
rz(-1.8701845) q[3];
sx q[3];
rz(-0.36417714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3889121) q[0];
sx q[0];
rz(-2.2081544) q[0];
sx q[0];
rz(2.3626784) q[0];
rz(2.3498416) q[1];
sx q[1];
rz(-1.2261032) q[1];
sx q[1];
rz(-0.34599272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7900171) q[0];
sx q[0];
rz(-1.6540961) q[0];
sx q[0];
rz(2.6935355) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6093639) q[2];
sx q[2];
rz(-2.8380605) q[2];
sx q[2];
rz(0.26798778) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4311667) q[1];
sx q[1];
rz(-1.0252608) q[1];
sx q[1];
rz(-3.1005076) q[1];
rz(-pi) q[2];
rz(2.2947512) q[3];
sx q[3];
rz(-1.8454285) q[3];
sx q[3];
rz(1.3534897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9927696) q[2];
sx q[2];
rz(-1.6042446) q[2];
sx q[2];
rz(2.0796622) q[2];
rz(-2.3289833) q[3];
sx q[3];
rz(-1.9989697) q[3];
sx q[3];
rz(0.44212166) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84520311) q[0];
sx q[0];
rz(-0.82798216) q[0];
sx q[0];
rz(1.6430586) q[0];
rz(-2.8614112) q[1];
sx q[1];
rz(-1.9277086) q[1];
sx q[1];
rz(2.3452506) q[1];
rz(1.610757) q[2];
sx q[2];
rz(-1.9794977) q[2];
sx q[2];
rz(-0.13499311) q[2];
rz(-0.82874684) q[3];
sx q[3];
rz(-1.256905) q[3];
sx q[3];
rz(-1.4245413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
