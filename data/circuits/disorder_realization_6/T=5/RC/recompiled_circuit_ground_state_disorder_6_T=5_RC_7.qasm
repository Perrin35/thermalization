OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6239983) q[0];
sx q[0];
rz(-0.52596337) q[0];
sx q[0];
rz(-2.9232803) q[0];
rz(1.4766308) q[1];
sx q[1];
rz(-2.6638439) q[1];
sx q[1];
rz(-0.55396095) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71611878) q[0];
sx q[0];
rz(-2.0332575) q[0];
sx q[0];
rz(1.0008971) q[0];
rz(1.5114097) q[2];
sx q[2];
rz(-2.0989387) q[2];
sx q[2];
rz(2.8092761) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3240149) q[1];
sx q[1];
rz(-2.2715873) q[1];
sx q[1];
rz(0.27591095) q[1];
x q[2];
rz(1.9945108) q[3];
sx q[3];
rz(-2.0689575) q[3];
sx q[3];
rz(2.8776282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7654984) q[2];
sx q[2];
rz(-0.29355294) q[2];
sx q[2];
rz(1.2676839) q[2];
rz(0.82967657) q[3];
sx q[3];
rz(-1.6206348) q[3];
sx q[3];
rz(2.7177496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053452881) q[0];
sx q[0];
rz(-1.3351048) q[0];
sx q[0];
rz(-1.7470737) q[0];
rz(2.246619) q[1];
sx q[1];
rz(-2.112969) q[1];
sx q[1];
rz(-1.5825533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6532063) q[0];
sx q[0];
rz(-1.0307285) q[0];
sx q[0];
rz(-2.3652943) q[0];
rz(0.3560654) q[2];
sx q[2];
rz(-1.967397) q[2];
sx q[2];
rz(0.44677904) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6239418) q[1];
sx q[1];
rz(-1.2576767) q[1];
sx q[1];
rz(-2.4409358) q[1];
rz(-pi) q[2];
rz(0.96554324) q[3];
sx q[3];
rz(-0.50600921) q[3];
sx q[3];
rz(0.11270302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6643657) q[2];
sx q[2];
rz(-1.619728) q[2];
sx q[2];
rz(-0.030755432) q[2];
rz(-2.6886046) q[3];
sx q[3];
rz(-2.9121297) q[3];
sx q[3];
rz(-2.9636813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40193108) q[0];
sx q[0];
rz(-0.10098305) q[0];
sx q[0];
rz(0.82823753) q[0];
rz(0.047686934) q[1];
sx q[1];
rz(-0.86367718) q[1];
sx q[1];
rz(-1.2275068) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9102893) q[0];
sx q[0];
rz(-2.4132015) q[0];
sx q[0];
rz(-2.4475218) q[0];
rz(-pi) q[1];
rz(2.9256608) q[2];
sx q[2];
rz(-1.1732444) q[2];
sx q[2];
rz(0.88653195) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.28042284) q[1];
sx q[1];
rz(-0.52045989) q[1];
sx q[1];
rz(-2.1153482) q[1];
x q[2];
rz(-1.5121721) q[3];
sx q[3];
rz(-2.4120962) q[3];
sx q[3];
rz(-0.54704715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.092827169) q[2];
sx q[2];
rz(-2.2266677) q[2];
sx q[2];
rz(-1.1009334) q[2];
rz(-0.01384211) q[3];
sx q[3];
rz(-1.3661386) q[3];
sx q[3];
rz(-0.83465105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58532995) q[0];
sx q[0];
rz(-1.7381373) q[0];
sx q[0];
rz(0.334326) q[0];
rz(-2.4090134) q[1];
sx q[1];
rz(-2.1926447) q[1];
sx q[1];
rz(0.11071959) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3667612) q[0];
sx q[0];
rz(-0.80935055) q[0];
sx q[0];
rz(1.5917718) q[0];
rz(-pi) q[1];
rz(-0.16558318) q[2];
sx q[2];
rz(-1.6628169) q[2];
sx q[2];
rz(1.657287) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9423805) q[1];
sx q[1];
rz(-0.95006493) q[1];
sx q[1];
rz(-2.278028) q[1];
rz(0.48583123) q[3];
sx q[3];
rz(-2.3049092) q[3];
sx q[3];
rz(-2.776718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.683814) q[2];
sx q[2];
rz(-2.8595698) q[2];
sx q[2];
rz(-1.4078183) q[2];
rz(1.1553361) q[3];
sx q[3];
rz(-1.9852394) q[3];
sx q[3];
rz(-2.9467764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4487576) q[0];
sx q[0];
rz(-2.223707) q[0];
sx q[0];
rz(-0.99739972) q[0];
rz(-1.6150486) q[1];
sx q[1];
rz(-2.5020182) q[1];
sx q[1];
rz(-1.7220727) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2888694) q[0];
sx q[0];
rz(-2.8670681) q[0];
sx q[0];
rz(2.2122266) q[0];
x q[1];
rz(-2.9827955) q[2];
sx q[2];
rz(-2.3613075) q[2];
sx q[2];
rz(-2.4312388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43806009) q[1];
sx q[1];
rz(-2.4171481) q[1];
sx q[1];
rz(-2.0634335) q[1];
rz(3.0180198) q[3];
sx q[3];
rz(-1.6032748) q[3];
sx q[3];
rz(2.6952254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2898966) q[2];
sx q[2];
rz(-0.09446129) q[2];
sx q[2];
rz(-0.32290253) q[2];
rz(-1.1139392) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(-0.41306257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36393976) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(-0.80663484) q[0];
rz(-0.58397645) q[1];
sx q[1];
rz(-1.1176502) q[1];
sx q[1];
rz(-1.1899828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89878824) q[0];
sx q[0];
rz(-1.6938279) q[0];
sx q[0];
rz(-1.2205628) q[0];
rz(-pi) q[1];
rz(0.34752589) q[2];
sx q[2];
rz(-0.56171562) q[2];
sx q[2];
rz(0.39081854) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2417702) q[1];
sx q[1];
rz(-0.17647753) q[1];
sx q[1];
rz(-1.8679669) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7403931) q[3];
sx q[3];
rz(-0.84753321) q[3];
sx q[3];
rz(0.71454988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40818647) q[2];
sx q[2];
rz(-1.5032282) q[2];
sx q[2];
rz(2.9564986) q[2];
rz(1.5271651) q[3];
sx q[3];
rz(-1.7388758) q[3];
sx q[3];
rz(2.4534498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9516893) q[0];
sx q[0];
rz(-2.1898495) q[0];
sx q[0];
rz(0.21251799) q[0];
rz(-2.3566133) q[1];
sx q[1];
rz(-1.4947944) q[1];
sx q[1];
rz(2.3627538) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76701421) q[0];
sx q[0];
rz(-2.2901669) q[0];
sx q[0];
rz(2.3423751) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5539407) q[2];
sx q[2];
rz(-1.8853123) q[2];
sx q[2];
rz(-1.8314198) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96152119) q[1];
sx q[1];
rz(-0.33726443) q[1];
sx q[1];
rz(2.2757761) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8605729) q[3];
sx q[3];
rz(-1.1906151) q[3];
sx q[3];
rz(0.036503867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51398858) q[2];
sx q[2];
rz(-2.6476761) q[2];
sx q[2];
rz(-2.7331875) q[2];
rz(-2.4397395) q[3];
sx q[3];
rz(-0.93188325) q[3];
sx q[3];
rz(-1.2592038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34847611) q[0];
sx q[0];
rz(-1.2154673) q[0];
sx q[0];
rz(-0.066019639) q[0];
rz(-1.6325715) q[1];
sx q[1];
rz(-1.0450109) q[1];
sx q[1];
rz(0.95796934) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45487472) q[0];
sx q[0];
rz(-1.2871847) q[0];
sx q[0];
rz(-0.4768178) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1195378) q[2];
sx q[2];
rz(-2.3081452) q[2];
sx q[2];
rz(1.1481783) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.66855318) q[1];
sx q[1];
rz(-0.91382342) q[1];
sx q[1];
rz(1.1725575) q[1];
x q[2];
rz(-2.1470966) q[3];
sx q[3];
rz(-2.9589783) q[3];
sx q[3];
rz(1.2507358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8336739) q[2];
sx q[2];
rz(-1.5235528) q[2];
sx q[2];
rz(1.4109122) q[2];
rz(0.69502568) q[3];
sx q[3];
rz(-1.4796939) q[3];
sx q[3];
rz(-0.01785774) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0787635) q[0];
sx q[0];
rz(-1.48209) q[0];
sx q[0];
rz(0.98989809) q[0];
rz(-0.46317378) q[1];
sx q[1];
rz(-1.4521234) q[1];
sx q[1];
rz(1.90082) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1027066) q[0];
sx q[0];
rz(-0.30065824) q[0];
sx q[0];
rz(-1.7666398) q[0];
rz(-pi) q[1];
rz(-1.1616917) q[2];
sx q[2];
rz(-0.4767524) q[2];
sx q[2];
rz(-2.7331405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.040673469) q[1];
sx q[1];
rz(-1.8534436) q[1];
sx q[1];
rz(1.2987178) q[1];
x q[2];
rz(-2.1911591) q[3];
sx q[3];
rz(-1.9470805) q[3];
sx q[3];
rz(2.2971414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3339633) q[2];
sx q[2];
rz(-1.3796076) q[2];
sx q[2];
rz(-2.4228952) q[2];
rz(-1.1527609) q[3];
sx q[3];
rz(-1.4216239) q[3];
sx q[3];
rz(1.0144455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.5929247) q[0];
sx q[0];
rz(-2.7935226) q[0];
sx q[0];
rz(0.80192178) q[0];
rz(-2.0536664) q[1];
sx q[1];
rz(-1.8049003) q[1];
sx q[1];
rz(0.71802872) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73488) q[0];
sx q[0];
rz(-2.3588484) q[0];
sx q[0];
rz(-0.46371622) q[0];
rz(-pi) q[1];
rz(2.6177216) q[2];
sx q[2];
rz(-1.1351368) q[2];
sx q[2];
rz(-2.4031529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7321328) q[1];
sx q[1];
rz(-1.902305) q[1];
sx q[1];
rz(0.055268754) q[1];
rz(-0.78022782) q[3];
sx q[3];
rz(-2.2426668) q[3];
sx q[3];
rz(-1.5742009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.7903018) q[2];
sx q[2];
rz(-0.21778926) q[2];
sx q[2];
rz(-0.51266074) q[2];
rz(2.3971108) q[3];
sx q[3];
rz(-2.1148041) q[3];
sx q[3];
rz(2.3775533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9164593) q[0];
sx q[0];
rz(-1.8712578) q[0];
sx q[0];
rz(1.901392) q[0];
rz(2.4254639) q[1];
sx q[1];
rz(-0.54812535) q[1];
sx q[1];
rz(-2.5352238) q[1];
rz(-3.0307583) q[2];
sx q[2];
rz(-1.7946984) q[2];
sx q[2];
rz(2.6738965) q[2];
rz(1.8567139) q[3];
sx q[3];
rz(-2.7919522) q[3];
sx q[3];
rz(-1.5408564) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
