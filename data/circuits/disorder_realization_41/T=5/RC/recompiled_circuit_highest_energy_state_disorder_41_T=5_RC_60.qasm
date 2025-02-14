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
rz(-0.77223414) q[0];
sx q[0];
rz(5.0679591) q[0];
sx q[0];
rz(8.2740347) q[0];
rz(2.8590705) q[1];
sx q[1];
rz(-1.1257659) q[1];
sx q[1];
rz(-1.1123302) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5939358) q[0];
sx q[0];
rz(-2.8540683) q[0];
sx q[0];
rz(2.7178746) q[0];
x q[1];
rz(3.0289828) q[2];
sx q[2];
rz(-2.1499584) q[2];
sx q[2];
rz(0.65161588) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9291275) q[1];
sx q[1];
rz(-2.0501137) q[1];
sx q[1];
rz(2.1838837) q[1];
x q[2];
rz(-0.78948094) q[3];
sx q[3];
rz(-1.0320569) q[3];
sx q[3];
rz(2.1110502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34899601) q[2];
sx q[2];
rz(-2.0670321) q[2];
sx q[2];
rz(1.3953155) q[2];
rz(-3.0830749) q[3];
sx q[3];
rz(-1.4165712) q[3];
sx q[3];
rz(-1.831656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7175452) q[0];
sx q[0];
rz(-2.738364) q[0];
sx q[0];
rz(-0.44977093) q[0];
rz(-0.46478081) q[1];
sx q[1];
rz(-1.8315146) q[1];
sx q[1];
rz(-2.8890077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9802482) q[0];
sx q[0];
rz(-0.35156116) q[0];
sx q[0];
rz(-1.5241966) q[0];
x q[1];
rz(-0.54754642) q[2];
sx q[2];
rz(-0.73807588) q[2];
sx q[2];
rz(2.944327) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8014268) q[1];
sx q[1];
rz(-2.8599842) q[1];
sx q[1];
rz(1.5605218) q[1];
rz(-pi) q[2];
rz(-2.3303495) q[3];
sx q[3];
rz(-1.1170804) q[3];
sx q[3];
rz(0.32920255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8978591) q[2];
sx q[2];
rz(-2.2458138) q[2];
sx q[2];
rz(-0.016156999) q[2];
rz(-2.2195418) q[3];
sx q[3];
rz(-2.1478896) q[3];
sx q[3];
rz(-0.13185681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.048412662) q[0];
sx q[0];
rz(-1.484363) q[0];
sx q[0];
rz(-0.68159252) q[0];
rz(-2.5438578) q[1];
sx q[1];
rz(-1.455541) q[1];
sx q[1];
rz(0.32424232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3675243) q[0];
sx q[0];
rz(-1.2307967) q[0];
sx q[0];
rz(2.6794479) q[0];
rz(-pi) q[1];
rz(0.50679548) q[2];
sx q[2];
rz(-1.5118216) q[2];
sx q[2];
rz(0.61117327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5983008) q[1];
sx q[1];
rz(-2.9678015) q[1];
sx q[1];
rz(1.325118) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3310166) q[3];
sx q[3];
rz(-0.62285715) q[3];
sx q[3];
rz(-1.9641124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9165667) q[2];
sx q[2];
rz(-2.7498507) q[2];
sx q[2];
rz(2.0908835) q[2];
rz(-0.74058908) q[3];
sx q[3];
rz(-1.2935484) q[3];
sx q[3];
rz(0.061323969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4968313) q[0];
sx q[0];
rz(-0.83503857) q[0];
sx q[0];
rz(0.9541676) q[0];
rz(2.0046115) q[1];
sx q[1];
rz(-1.515772) q[1];
sx q[1];
rz(2.383393) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0643142) q[0];
sx q[0];
rz(-1.2998733) q[0];
sx q[0];
rz(-3.0999557) q[0];
rz(-pi) q[1];
rz(1.3761282) q[2];
sx q[2];
rz(-2.3050781) q[2];
sx q[2];
rz(2.0912974) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41685795) q[1];
sx q[1];
rz(-1.6685889) q[1];
sx q[1];
rz(-0.87585249) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5810892) q[3];
sx q[3];
rz(-2.309628) q[3];
sx q[3];
rz(0.21496274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.053146426) q[2];
sx q[2];
rz(-0.55067486) q[2];
sx q[2];
rz(-1.3147563) q[2];
rz(2.8407319) q[3];
sx q[3];
rz(-1.8010537) q[3];
sx q[3];
rz(-0.68527591) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9984197) q[0];
sx q[0];
rz(-1.5851861) q[0];
sx q[0];
rz(1.5013303) q[0];
rz(0.38652626) q[1];
sx q[1];
rz(-1.0685579) q[1];
sx q[1];
rz(1.5949257) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37171897) q[0];
sx q[0];
rz(-1.3981016) q[0];
sx q[0];
rz(0.70575373) q[0];
rz(-pi) q[1];
rz(1.1207259) q[2];
sx q[2];
rz(-2.483027) q[2];
sx q[2];
rz(-0.77641247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1569231) q[1];
sx q[1];
rz(-0.44141967) q[1];
sx q[1];
rz(2.943931) q[1];
rz(-pi) q[2];
rz(0.076885496) q[3];
sx q[3];
rz(-2.5664483) q[3];
sx q[3];
rz(0.98371668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3950562) q[2];
sx q[2];
rz(-0.48607963) q[2];
sx q[2];
rz(-1.0699832) q[2];
rz(0.88548958) q[3];
sx q[3];
rz(-1.0642137) q[3];
sx q[3];
rz(1.1522256) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2162061) q[0];
sx q[0];
rz(-2.8031271) q[0];
sx q[0];
rz(2.0881407) q[0];
rz(2.9523051) q[1];
sx q[1];
rz(-2.3129251) q[1];
sx q[1];
rz(-1.4922356) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5976335) q[0];
sx q[0];
rz(-1.2043556) q[0];
sx q[0];
rz(-2.9935915) q[0];
x q[1];
rz(-1.8771421) q[2];
sx q[2];
rz(-1.7736846) q[2];
sx q[2];
rz(0.63837793) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.63000667) q[1];
sx q[1];
rz(-1.6701856) q[1];
sx q[1];
rz(-2.6376827) q[1];
rz(-1.6118649) q[3];
sx q[3];
rz(-2.3784935) q[3];
sx q[3];
rz(1.3535318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3377127) q[2];
sx q[2];
rz(-2.1264075) q[2];
sx q[2];
rz(-1.052915) q[2];
rz(3.0247011) q[3];
sx q[3];
rz(-2.0785073) q[3];
sx q[3];
rz(-1.6803928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18148024) q[0];
sx q[0];
rz(-0.57992613) q[0];
sx q[0];
rz(2.5970698) q[0];
rz(-2.8879884) q[1];
sx q[1];
rz(-2.8896152) q[1];
sx q[1];
rz(0.25467083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38886759) q[0];
sx q[0];
rz(-1.6427186) q[0];
sx q[0];
rz(-3.1353463) q[0];
rz(-pi) q[1];
rz(-2.7484863) q[2];
sx q[2];
rz(-1.4908893) q[2];
sx q[2];
rz(-0.42482947) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54714845) q[1];
sx q[1];
rz(-2.3871506) q[1];
sx q[1];
rz(-1.6625089) q[1];
x q[2];
rz(-2.5768877) q[3];
sx q[3];
rz(-0.89313358) q[3];
sx q[3];
rz(1.1260179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9165245) q[2];
sx q[2];
rz(-0.93776408) q[2];
sx q[2];
rz(0.41898215) q[2];
rz(-0.032912832) q[3];
sx q[3];
rz(-1.2337647) q[3];
sx q[3];
rz(2.029443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58881775) q[0];
sx q[0];
rz(-0.7426312) q[0];
sx q[0];
rz(1.8137929) q[0];
rz(1.0325507) q[1];
sx q[1];
rz(-1.3464709) q[1];
sx q[1];
rz(0.58951497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96846554) q[0];
sx q[0];
rz(-1.0778946) q[0];
sx q[0];
rz(2.1987134) q[0];
rz(-1.3016745) q[2];
sx q[2];
rz(-1.8523714) q[2];
sx q[2];
rz(-2.9219404) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8679598) q[1];
sx q[1];
rz(-2.3936749) q[1];
sx q[1];
rz(0.91435097) q[1];
rz(-pi) q[2];
rz(0.39535268) q[3];
sx q[3];
rz(-1.1046003) q[3];
sx q[3];
rz(2.6564997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0716268) q[2];
sx q[2];
rz(-0.67487851) q[2];
sx q[2];
rz(0.91342941) q[2];
rz(2.6895798) q[3];
sx q[3];
rz(-1.2830696) q[3];
sx q[3];
rz(1.3445331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2695059) q[0];
sx q[0];
rz(-1.0259314) q[0];
sx q[0];
rz(-1.5103229) q[0];
rz(1.3144846) q[1];
sx q[1];
rz(-1.038082) q[1];
sx q[1];
rz(0.41183919) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.340419) q[0];
sx q[0];
rz(-0.95511758) q[0];
sx q[0];
rz(2.8007617) q[0];
rz(-pi) q[1];
rz(-2.9963295) q[2];
sx q[2];
rz(-0.32215873) q[2];
sx q[2];
rz(1.8069428) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.50166124) q[1];
sx q[1];
rz(-1.2218214) q[1];
sx q[1];
rz(0.56569143) q[1];
x q[2];
rz(0.69782902) q[3];
sx q[3];
rz(-0.48664722) q[3];
sx q[3];
rz(-0.4781639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1371548) q[2];
sx q[2];
rz(-1.6496481) q[2];
sx q[2];
rz(2.906666) q[2];
rz(2.2233985) q[3];
sx q[3];
rz(-0.72489679) q[3];
sx q[3];
rz(1.0880067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10201564) q[0];
sx q[0];
rz(-0.35440847) q[0];
sx q[0];
rz(-1.8812195) q[0];
rz(1.4541939) q[1];
sx q[1];
rz(-1.6834384) q[1];
sx q[1];
rz(1.6967324) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4871976) q[0];
sx q[0];
rz(-0.62614589) q[0];
sx q[0];
rz(0.74855133) q[0];
x q[1];
rz(0.52839519) q[2];
sx q[2];
rz(-0.95167347) q[2];
sx q[2];
rz(0.88232679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9493172) q[1];
sx q[1];
rz(-2.2819464) q[1];
sx q[1];
rz(-1.7356731) q[1];
rz(-pi) q[2];
rz(1.904018) q[3];
sx q[3];
rz(-1.929627) q[3];
sx q[3];
rz(2.8976208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5626278) q[2];
sx q[2];
rz(-2.8426888) q[2];
sx q[2];
rz(-0.12358269) q[2];
rz(0.28892162) q[3];
sx q[3];
rz(-1.2762504) q[3];
sx q[3];
rz(2.6821274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4041483) q[0];
sx q[0];
rz(-1.7973719) q[0];
sx q[0];
rz(1.6112882) q[0];
rz(-0.96724802) q[1];
sx q[1];
rz(-0.9728685) q[1];
sx q[1];
rz(-2.2546993) q[1];
rz(2.5423302) q[2];
sx q[2];
rz(-1.5256554) q[2];
sx q[2];
rz(-1.2427606) q[2];
rz(-2.5206) q[3];
sx q[3];
rz(-0.98747333) q[3];
sx q[3];
rz(2.3458521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
