OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.53192294) q[0];
sx q[0];
rz(1.6165531) q[0];
sx q[0];
rz(11.308148) q[0];
rz(7.8426709) q[1];
sx q[1];
rz(5.7962228) q[1];
sx q[1];
rz(16.143057) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4420051) q[0];
sx q[0];
rz(-2.3908092) q[0];
sx q[0];
rz(2.4155099) q[0];
rz(-2.6003709) q[2];
sx q[2];
rz(-2.5337484) q[2];
sx q[2];
rz(1.7153049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5158029) q[1];
sx q[1];
rz(-2.8608659) q[1];
sx q[1];
rz(1.2899542) q[1];
rz(-pi) q[2];
rz(-0.55032309) q[3];
sx q[3];
rz(-1.9903515) q[3];
sx q[3];
rz(1.5997831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9293999) q[2];
sx q[2];
rz(-1.8636899) q[2];
sx q[2];
rz(-2.0983992) q[2];
rz(-2.7405401) q[3];
sx q[3];
rz(-1.4570823) q[3];
sx q[3];
rz(0.57798398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9143518) q[0];
sx q[0];
rz(-0.14509097) q[0];
sx q[0];
rz(-2.5703854) q[0];
rz(-2.1696443) q[1];
sx q[1];
rz(-0.74058878) q[1];
sx q[1];
rz(-2.0170225) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1069431) q[0];
sx q[0];
rz(-0.9616344) q[0];
sx q[0];
rz(-2.8657495) q[0];
x q[1];
rz(0.45887453) q[2];
sx q[2];
rz(-0.70415184) q[2];
sx q[2];
rz(1.5062694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3414866) q[1];
sx q[1];
rz(-2.3012495) q[1];
sx q[1];
rz(-3.0154702) q[1];
x q[2];
rz(1.589682) q[3];
sx q[3];
rz(-1.7715997) q[3];
sx q[3];
rz(2.7310179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9198415) q[2];
sx q[2];
rz(-1.5566748) q[2];
sx q[2];
rz(2.705503) q[2];
rz(-0.74887216) q[3];
sx q[3];
rz(-2.336453) q[3];
sx q[3];
rz(0.82912904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1083168) q[0];
sx q[0];
rz(-1.3949787) q[0];
sx q[0];
rz(-0.088031553) q[0];
rz(-2.2875359) q[1];
sx q[1];
rz(-0.66103926) q[1];
sx q[1];
rz(2.0565654) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2449834) q[0];
sx q[0];
rz(-1.5036146) q[0];
sx q[0];
rz(-1.2981197) q[0];
rz(2.1846071) q[2];
sx q[2];
rz(-0.26273604) q[2];
sx q[2];
rz(-1.764591) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2287187) q[1];
sx q[1];
rz(-1.7342028) q[1];
sx q[1];
rz(-2.0151369) q[1];
rz(2.964193) q[3];
sx q[3];
rz(-2.7466603) q[3];
sx q[3];
rz(1.4876786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4714841) q[2];
sx q[2];
rz(-1.4205616) q[2];
sx q[2];
rz(-2.3954083) q[2];
rz(0.21607312) q[3];
sx q[3];
rz(-1.8862628) q[3];
sx q[3];
rz(-2.9358673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3716607) q[0];
sx q[0];
rz(-3.0682204) q[0];
sx q[0];
rz(1.368847) q[0];
rz(-0.61028496) q[1];
sx q[1];
rz(-1.7758324) q[1];
sx q[1];
rz(2.761421) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1606522) q[0];
sx q[0];
rz(-1.9987172) q[0];
sx q[0];
rz(1.9141657) q[0];
rz(-pi) q[1];
rz(-1.3743782) q[2];
sx q[2];
rz(-2.0140208) q[2];
sx q[2];
rz(0.049510844) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7138765) q[1];
sx q[1];
rz(-0.52449924) q[1];
sx q[1];
rz(-1.8668411) q[1];
rz(-pi) q[2];
rz(0.065297619) q[3];
sx q[3];
rz(-1.9717798) q[3];
sx q[3];
rz(-2.9915265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88468203) q[2];
sx q[2];
rz(-2.1941049) q[2];
sx q[2];
rz(1.0379637) q[2];
rz(2.7773618) q[3];
sx q[3];
rz(-2.3528152) q[3];
sx q[3];
rz(0.68575931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6895741) q[0];
sx q[0];
rz(-1.4301393) q[0];
sx q[0];
rz(0.83056393) q[0];
rz(-3.0914302) q[1];
sx q[1];
rz(-3.0144033) q[1];
sx q[1];
rz(-0.62279472) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1019331) q[0];
sx q[0];
rz(-0.23124157) q[0];
sx q[0];
rz(2.7949265) q[0];
x q[1];
rz(2.7217614) q[2];
sx q[2];
rz(-1.3908885) q[2];
sx q[2];
rz(0.77984389) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.894815) q[1];
sx q[1];
rz(-2.3843002) q[1];
sx q[1];
rz(-1.0411109) q[1];
x q[2];
rz(1.5708617) q[3];
sx q[3];
rz(-1.5483466) q[3];
sx q[3];
rz(-0.27014905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.043776) q[2];
sx q[2];
rz(-2.706683) q[2];
sx q[2];
rz(-2.3338976) q[2];
rz(1.5378599) q[3];
sx q[3];
rz(-1.5065498) q[3];
sx q[3];
rz(0.22411552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1767126) q[0];
sx q[0];
rz(-1.8105312) q[0];
sx q[0];
rz(1.4759195) q[0];
rz(-1.2337947) q[1];
sx q[1];
rz(-1.7314311) q[1];
sx q[1];
rz(0.15353157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4013256) q[0];
sx q[0];
rz(-3.0318878) q[0];
sx q[0];
rz(-0.033880635) q[0];
rz(-pi) q[1];
rz(2.0664882) q[2];
sx q[2];
rz(-1.9736991) q[2];
sx q[2];
rz(-0.35225454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.081806) q[1];
sx q[1];
rz(-2.012552) q[1];
sx q[1];
rz(-0.55952832) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5028822) q[3];
sx q[3];
rz(-1.494246) q[3];
sx q[3];
rz(-2.9552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4393356) q[2];
sx q[2];
rz(-0.80436891) q[2];
sx q[2];
rz(0.55629936) q[2];
rz(3.124681) q[3];
sx q[3];
rz(-2.1664679) q[3];
sx q[3];
rz(2.4771966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.21338129) q[0];
sx q[0];
rz(-3.0688372) q[0];
sx q[0];
rz(-0.67759204) q[0];
rz(-1.821359) q[1];
sx q[1];
rz(-1.0878539) q[1];
sx q[1];
rz(-2.5755612) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2293733) q[0];
sx q[0];
rz(-1.1347596) q[0];
sx q[0];
rz(-1.9422117) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2957057) q[2];
sx q[2];
rz(-3.0873211) q[2];
sx q[2];
rz(1.4473297) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2905419) q[1];
sx q[1];
rz(-2.4220963) q[1];
sx q[1];
rz(-1.8012723) q[1];
rz(1.612408) q[3];
sx q[3];
rz(-2.0609612) q[3];
sx q[3];
rz(-0.13652882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90938202) q[2];
sx q[2];
rz(-2.7558694) q[2];
sx q[2];
rz(-2.5717226) q[2];
rz(1.0424403) q[3];
sx q[3];
rz(-1.180155) q[3];
sx q[3];
rz(-1.0038143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7637699) q[0];
sx q[0];
rz(-0.42709392) q[0];
sx q[0];
rz(1.0231934) q[0];
rz(-0.91776735) q[1];
sx q[1];
rz(-2.387391) q[1];
sx q[1];
rz(2.2775211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3765434) q[0];
sx q[0];
rz(-1.6934388) q[0];
sx q[0];
rz(-1.3255694) q[0];
rz(-pi) q[1];
rz(2.239924) q[2];
sx q[2];
rz(-2.5732917) q[2];
sx q[2];
rz(-1.2327884) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0297599) q[1];
sx q[1];
rz(-2.2297292) q[1];
sx q[1];
rz(-0.62739851) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.731063) q[3];
sx q[3];
rz(-0.77068146) q[3];
sx q[3];
rz(0.16134027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0920022) q[2];
sx q[2];
rz(-1.2668173) q[2];
sx q[2];
rz(2.3769296) q[2];
rz(2.2406254) q[3];
sx q[3];
rz(-2.4668906) q[3];
sx q[3];
rz(-1.0990134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6242999) q[0];
sx q[0];
rz(-2.7462672) q[0];
sx q[0];
rz(-0.98315352) q[0];
rz(0.82955018) q[1];
sx q[1];
rz(-1.7770276) q[1];
sx q[1];
rz(0.57873908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.549141) q[0];
sx q[0];
rz(-1.5835174) q[0];
sx q[0];
rz(-1.1770269) q[0];
x q[1];
rz(-0.77909536) q[2];
sx q[2];
rz(-2.4367174) q[2];
sx q[2];
rz(1.5388956) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14565258) q[1];
sx q[1];
rz(-2.4697587) q[1];
sx q[1];
rz(2.2493659) q[1];
rz(-pi) q[2];
rz(1.6074568) q[3];
sx q[3];
rz(-1.9158773) q[3];
sx q[3];
rz(-0.35547977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65043727) q[2];
sx q[2];
rz(-1.2885685) q[2];
sx q[2];
rz(-3.0299661) q[2];
rz(-0.95587436) q[3];
sx q[3];
rz(-0.99146944) q[3];
sx q[3];
rz(1.8808232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2535506) q[0];
sx q[0];
rz(-2.084806) q[0];
sx q[0];
rz(3.1274617) q[0];
rz(1.1125394) q[1];
sx q[1];
rz(-0.91347778) q[1];
sx q[1];
rz(-1.6003476) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3936716) q[0];
sx q[0];
rz(-1.4427516) q[0];
sx q[0];
rz(-2.9927954) q[0];
rz(-0.92629536) q[2];
sx q[2];
rz(-1.4699114) q[2];
sx q[2];
rz(2.6580194) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83859119) q[1];
sx q[1];
rz(-1.934774) q[1];
sx q[1];
rz(1.6345657) q[1];
x q[2];
rz(-2.3722911) q[3];
sx q[3];
rz(-0.23450867) q[3];
sx q[3];
rz(-2.086139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9666226) q[2];
sx q[2];
rz(-1.9113767) q[2];
sx q[2];
rz(-0.54565412) q[2];
rz(-0.086006554) q[3];
sx q[3];
rz(-1.6274446) q[3];
sx q[3];
rz(2.8470794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94588146) q[0];
sx q[0];
rz(-2.0989037) q[0];
sx q[0];
rz(1.209191) q[0];
rz(-0.60278268) q[1];
sx q[1];
rz(-0.61059112) q[1];
sx q[1];
rz(0.75377725) q[1];
rz(1.3134342) q[2];
sx q[2];
rz(-0.69862294) q[2];
sx q[2];
rz(-0.84702484) q[2];
rz(-2.9908441) q[3];
sx q[3];
rz(-0.6282395) q[3];
sx q[3];
rz(0.30368526) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
