OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(-2.6053083) q[0];
sx q[0];
rz(0.94776881) q[0];
rz(-1.3287969) q[1];
sx q[1];
rz(4.4089945) q[1];
sx q[1];
rz(10.452527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9337024) q[0];
sx q[0];
rz(-0.37651248) q[0];
sx q[0];
rz(-0.062113751) q[0];
x q[1];
rz(1.1711636) q[2];
sx q[2];
rz(-2.6174449) q[2];
sx q[2];
rz(-1.5607967) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9103968) q[1];
sx q[1];
rz(-0.80363552) q[1];
sx q[1];
rz(2.5956144) q[1];
rz(2.5146033) q[3];
sx q[3];
rz(-1.1999745) q[3];
sx q[3];
rz(0.8644608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1296922) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(-1.0502846) q[2];
rz(-2.0283279) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691836) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(-0.29775277) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(2.0334977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1381123) q[0];
sx q[0];
rz(-1.7797884) q[0];
sx q[0];
rz(-0.94603993) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96252243) q[2];
sx q[2];
rz(-1.9307185) q[2];
sx q[2];
rz(-0.72812176) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88156466) q[1];
sx q[1];
rz(-1.2177055) q[1];
sx q[1];
rz(-0.41207037) q[1];
rz(-pi) q[2];
rz(-1.3012582) q[3];
sx q[3];
rz(-2.2552935) q[3];
sx q[3];
rz(1.0052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46253282) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(-0.97529808) q[2];
rz(-0.9179999) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179203) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(-0.54291022) q[0];
rz(-2.2593598) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-0.96484819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1436413) q[0];
sx q[0];
rz(-2.7086341) q[0];
sx q[0];
rz(2.4233682) q[0];
rz(-0.13871128) q[2];
sx q[2];
rz(-2.841946) q[2];
sx q[2];
rz(-1.2219929) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45005349) q[1];
sx q[1];
rz(-1.4496526) q[1];
sx q[1];
rz(-0.064784592) q[1];
x q[2];
rz(1.996702) q[3];
sx q[3];
rz(-0.45255462) q[3];
sx q[3];
rz(-1.9354613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0470011) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(-1.1331406) q[2];
rz(-2.9099693) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8621181) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(-1.7650771) q[0];
rz(2.6230985) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(-0.24212295) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69232363) q[0];
sx q[0];
rz(-1.4081435) q[0];
sx q[0];
rz(2.4641894) q[0];
rz(2.7119615) q[2];
sx q[2];
rz(-1.5589082) q[2];
sx q[2];
rz(2.5846543) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8072847) q[1];
sx q[1];
rz(-2.0172999) q[1];
sx q[1];
rz(1.3825033) q[1];
rz(1.4297156) q[3];
sx q[3];
rz(-1.9851079) q[3];
sx q[3];
rz(2.8431818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1233998) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(0.094853178) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78385329) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(-2.7807996) q[0];
rz(-1.7533253) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(2.0070019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.40588) q[0];
sx q[0];
rz(-2.1942733) q[0];
sx q[0];
rz(2.1924125) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1410494) q[2];
sx q[2];
rz(-1.2417214) q[2];
sx q[2];
rz(0.82440257) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6280917) q[1];
sx q[1];
rz(-0.47422945) q[1];
sx q[1];
rz(0.34631108) q[1];
rz(-pi) q[2];
rz(1.6547104) q[3];
sx q[3];
rz(-2.0062371) q[3];
sx q[3];
rz(-2.8403789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1363042) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(-0.57265442) q[2];
rz(0.92875656) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(-1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(2.8144433) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(-1.2231187) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(1.9893601) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1359065) q[0];
sx q[0];
rz(-1.5581157) q[0];
sx q[0];
rz(1.6760875) q[0];
x q[1];
rz(1.1587028) q[2];
sx q[2];
rz(-2.0372143) q[2];
sx q[2];
rz(-0.42883021) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8730948) q[1];
sx q[1];
rz(-2.9151306) q[1];
sx q[1];
rz(-0.48392673) q[1];
rz(-0.96529393) q[3];
sx q[3];
rz(-0.98635841) q[3];
sx q[3];
rz(-0.30129978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46889177) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836442) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(0.062285475) q[0];
rz(-2.9557872) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(2.7468162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5236854) q[0];
sx q[0];
rz(-0.46645188) q[0];
sx q[0];
rz(-2.594069) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6997509) q[2];
sx q[2];
rz(-2.664898) q[2];
sx q[2];
rz(1.9112019) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3697606) q[1];
sx q[1];
rz(-2.6409915) q[1];
sx q[1];
rz(-0.53497772) q[1];
x q[2];
rz(1.6077605) q[3];
sx q[3];
rz(-1.4174403) q[3];
sx q[3];
rz(1.6085898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-0.33005565) q[2];
sx q[2];
rz(0.27302343) q[2];
rz(-1.8388883) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(2.8295512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39772314) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(1.2851108) q[0];
rz(1.5015191) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1947945) q[0];
sx q[0];
rz(-1.1800982) q[0];
sx q[0];
rz(0.35238738) q[0];
rz(-pi) q[1];
rz(0.060872002) q[2];
sx q[2];
rz(-0.69673046) q[2];
sx q[2];
rz(-2.3827162) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8074639) q[1];
sx q[1];
rz(-0.51528105) q[1];
sx q[1];
rz(2.3433102) q[1];
rz(2.1741381) q[3];
sx q[3];
rz(-2.6199322) q[3];
sx q[3];
rz(1.1286917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0354707) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(-1.0236053) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(1.5203083) q[0];
rz(-0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(0.83713371) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9595813) q[0];
sx q[0];
rz(-0.89996979) q[0];
sx q[0];
rz(-1.3120033) q[0];
x q[1];
rz(1.2741954) q[2];
sx q[2];
rz(-1.0146078) q[2];
sx q[2];
rz(-1.1013168) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3763807) q[1];
sx q[1];
rz(-1.7227255) q[1];
sx q[1];
rz(0.10581776) q[1];
rz(-pi) q[2];
rz(2.8505441) q[3];
sx q[3];
rz(-1.0927199) q[3];
sx q[3];
rz(-2.2136798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5403486) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(-2.1458697) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(1.3109591) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(1.0572222) q[0];
rz(-3.0341042) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(-2.1616139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9419132) q[0];
sx q[0];
rz(-1.4451888) q[0];
sx q[0];
rz(0.051254) q[0];
x q[1];
rz(0.51771848) q[2];
sx q[2];
rz(-1.1283518) q[2];
sx q[2];
rz(-0.25416086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6085538) q[1];
sx q[1];
rz(-1.8095008) q[1];
sx q[1];
rz(-2.0134316) q[1];
x q[2];
rz(1.8944593) q[3];
sx q[3];
rz(-1.4753046) q[3];
sx q[3];
rz(3.132706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3048627) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(2.1137962) q[2];
rz(1.7547539) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(-2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733611) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(0.83203075) q[1];
sx q[1];
rz(-0.46453005) q[1];
sx q[1];
rz(0.66418905) q[1];
rz(-2.1297395) q[2];
sx q[2];
rz(-1.7218628) q[2];
sx q[2];
rz(-2.0414258) q[2];
rz(0.50073033) q[3];
sx q[3];
rz(-1.2320319) q[3];
sx q[3];
rz(-2.4945955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
