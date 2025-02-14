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
rz(1.2430159) q[0];
sx q[0];
rz(-1.1216811) q[0];
sx q[0];
rz(-2.6774874) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(2.2321489) q[1];
sx q[1];
rz(7.5997054) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1027604) q[0];
sx q[0];
rz(-1.6359207) q[0];
sx q[0];
rz(-1.536157) q[0];
rz(-3.0010536) q[2];
sx q[2];
rz(-0.34374434) q[2];
sx q[2];
rz(3.112971) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39151788) q[1];
sx q[1];
rz(-0.98528359) q[1];
sx q[1];
rz(-1.2445009) q[1];
rz(0.85342225) q[3];
sx q[3];
rz(-2.7629921) q[3];
sx q[3];
rz(-0.84747696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2509649) q[2];
sx q[2];
rz(-0.35718063) q[2];
sx q[2];
rz(-2.8446021) q[2];
rz(0.65303981) q[3];
sx q[3];
rz(-1.5584471) q[3];
sx q[3];
rz(-2.0853341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62946573) q[0];
sx q[0];
rz(-2.1194206) q[0];
sx q[0];
rz(-0.28489354) q[0];
rz(0.051305436) q[1];
sx q[1];
rz(-2.3785794) q[1];
sx q[1];
rz(0.6368534) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1276217) q[0];
sx q[0];
rz(-0.68249615) q[0];
sx q[0];
rz(-2.3861814) q[0];
x q[1];
rz(-0.35164386) q[2];
sx q[2];
rz(-1.571901) q[2];
sx q[2];
rz(2.0948834) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11733774) q[1];
sx q[1];
rz(-1.7452612) q[1];
sx q[1];
rz(1.4937449) q[1];
rz(-0.70592441) q[3];
sx q[3];
rz(-1.4691938) q[3];
sx q[3];
rz(1.3639579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9597943) q[2];
sx q[2];
rz(-2.2836298) q[2];
sx q[2];
rz(2.174343) q[2];
rz(2.8773384) q[3];
sx q[3];
rz(-1.5033009) q[3];
sx q[3];
rz(1.9056457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018175) q[0];
sx q[0];
rz(-1.5622666) q[0];
sx q[0];
rz(-1.1199957) q[0];
rz(2.6657875) q[1];
sx q[1];
rz(-2.0640524) q[1];
sx q[1];
rz(-2.8376104) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8037) q[0];
sx q[0];
rz(-2.5753382) q[0];
sx q[0];
rz(1.5192274) q[0];
x q[1];
rz(1.6725704) q[2];
sx q[2];
rz(-1.8074805) q[2];
sx q[2];
rz(-0.19381154) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5604104) q[1];
sx q[1];
rz(-0.85542233) q[1];
sx q[1];
rz(-0.13611273) q[1];
rz(1.7695565) q[3];
sx q[3];
rz(-2.244368) q[3];
sx q[3];
rz(-2.4186132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1659282) q[2];
sx q[2];
rz(-2.1037481) q[2];
sx q[2];
rz(-0.69592875) q[2];
rz(-1.1337229) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(-0.55293647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7021779) q[0];
sx q[0];
rz(-1.004383) q[0];
sx q[0];
rz(2.0942005) q[0];
rz(-1.9678496) q[1];
sx q[1];
rz(-2.4672697) q[1];
sx q[1];
rz(1.5210927) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58529638) q[0];
sx q[0];
rz(-1.4504787) q[0];
sx q[0];
rz(-0.5043982) q[0];
x q[1];
rz(1.1999454) q[2];
sx q[2];
rz(-2.4003138) q[2];
sx q[2];
rz(1.0560869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2167705) q[1];
sx q[1];
rz(-2.3433279) q[1];
sx q[1];
rz(-2.5096748) q[1];
rz(-pi) q[2];
rz(-1.2866045) q[3];
sx q[3];
rz(-1.2582558) q[3];
sx q[3];
rz(-0.72060637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.095470458) q[2];
sx q[2];
rz(-2.0746524) q[2];
sx q[2];
rz(-2.6037237) q[2];
rz(-1.4872023) q[3];
sx q[3];
rz(-2.7345149) q[3];
sx q[3];
rz(-2.4552104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9293514) q[0];
sx q[0];
rz(-2.90726) q[0];
sx q[0];
rz(0.60337639) q[0];
rz(0.76250184) q[1];
sx q[1];
rz(-1.1926032) q[1];
sx q[1];
rz(0.71294436) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0423316) q[0];
sx q[0];
rz(-1.5047899) q[0];
sx q[0];
rz(1.1724654) q[0];
rz(-pi) q[1];
rz(2.0335401) q[2];
sx q[2];
rz(-1.6959091) q[2];
sx q[2];
rz(-0.81415983) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.566943) q[1];
sx q[1];
rz(-1.7716367) q[1];
sx q[1];
rz(1.0667136) q[1];
rz(-pi) q[2];
rz(-0.6176881) q[3];
sx q[3];
rz(-1.0782982) q[3];
sx q[3];
rz(-2.4503051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1684299) q[2];
sx q[2];
rz(-0.91732401) q[2];
sx q[2];
rz(1.6707576) q[2];
rz(-0.51703185) q[3];
sx q[3];
rz(-1.0443338) q[3];
sx q[3];
rz(1.8487336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1626728) q[0];
sx q[0];
rz(-0.51386583) q[0];
sx q[0];
rz(0.55147076) q[0];
rz(0.035471352) q[1];
sx q[1];
rz(-2.008581) q[1];
sx q[1];
rz(1.5303401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0198675) q[0];
sx q[0];
rz(-0.95917738) q[0];
sx q[0];
rz(-0.83329552) q[0];
rz(-pi) q[1];
rz(-0.82357652) q[2];
sx q[2];
rz(-1.608143) q[2];
sx q[2];
rz(-2.3970791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6833466) q[1];
sx q[1];
rz(-0.85696942) q[1];
sx q[1];
rz(-0.48269646) q[1];
rz(-1.6812165) q[3];
sx q[3];
rz(-2.4856119) q[3];
sx q[3];
rz(-0.15925281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0772721) q[2];
sx q[2];
rz(-0.52017009) q[2];
sx q[2];
rz(1.8989829) q[2];
rz(0.51314917) q[3];
sx q[3];
rz(-0.41392252) q[3];
sx q[3];
rz(-1.3750403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4921017) q[0];
sx q[0];
rz(-2.212337) q[0];
sx q[0];
rz(3.0249366) q[0];
rz(-2.3684582) q[1];
sx q[1];
rz(-1.1583068) q[1];
sx q[1];
rz(-1.6366417) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3706544) q[0];
sx q[0];
rz(-0.8621093) q[0];
sx q[0];
rz(-2.7518726) q[0];
x q[1];
rz(-3.0819986) q[2];
sx q[2];
rz(-1.182654) q[2];
sx q[2];
rz(-1.3898894) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22681397) q[1];
sx q[1];
rz(-2.249301) q[1];
sx q[1];
rz(2.6752276) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7724593) q[3];
sx q[3];
rz(-2.4551486) q[3];
sx q[3];
rz(2.4001207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2152805) q[2];
sx q[2];
rz(-2.8947688) q[2];
sx q[2];
rz(-2.2543294) q[2];
rz(0.67982802) q[3];
sx q[3];
rz(-0.66893783) q[3];
sx q[3];
rz(0.63290709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.55650869) q[0];
sx q[0];
rz(-1.8483138) q[0];
sx q[0];
rz(0.65628091) q[0];
rz(0.20690021) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(-1.9237178) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6164056) q[0];
sx q[0];
rz(-1.5357247) q[0];
sx q[0];
rz(-3.0588018) q[0];
rz(-pi) q[1];
rz(1.0704667) q[2];
sx q[2];
rz(-0.64979759) q[2];
sx q[2];
rz(-1.4184679) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43133538) q[1];
sx q[1];
rz(-1.4568304) q[1];
sx q[1];
rz(2.9821787) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37533203) q[3];
sx q[3];
rz(-0.64556827) q[3];
sx q[3];
rz(-0.69824346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8021585) q[2];
sx q[2];
rz(-1.755244) q[2];
sx q[2];
rz(2.459724) q[2];
rz(-1.1784461) q[3];
sx q[3];
rz(-2.384187) q[3];
sx q[3];
rz(1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7449529) q[0];
sx q[0];
rz(-1.6093901) q[0];
sx q[0];
rz(-0.21016453) q[0];
rz(-0.22178966) q[1];
sx q[1];
rz(-0.91786018) q[1];
sx q[1];
rz(3.0992357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6307459) q[0];
sx q[0];
rz(-1.5486985) q[0];
sx q[0];
rz(0.18760292) q[0];
rz(-0.75889672) q[2];
sx q[2];
rz(-0.79115552) q[2];
sx q[2];
rz(0.10221419) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2376304) q[1];
sx q[1];
rz(-1.018486) q[1];
sx q[1];
rz(-1.2321074) q[1];
rz(2.6082392) q[3];
sx q[3];
rz(-1.1660398) q[3];
sx q[3];
rz(0.08431708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28045851) q[2];
sx q[2];
rz(-1.0886085) q[2];
sx q[2];
rz(-0.66656485) q[2];
rz(-0.76513964) q[3];
sx q[3];
rz(-2.9355526) q[3];
sx q[3];
rz(-1.7407181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6078981) q[0];
sx q[0];
rz(-0.38132897) q[0];
sx q[0];
rz(-2.4545942) q[0];
rz(-1.2835361) q[1];
sx q[1];
rz(-1.1524408) q[1];
sx q[1];
rz(-0.57354617) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4592996) q[0];
sx q[0];
rz(-0.98381616) q[0];
sx q[0];
rz(-0.94265509) q[0];
rz(-1.1185455) q[2];
sx q[2];
rz(-0.49511038) q[2];
sx q[2];
rz(-0.64436382) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0189218) q[1];
sx q[1];
rz(-1.9098567) q[1];
sx q[1];
rz(-1.0215205) q[1];
rz(0.76400842) q[3];
sx q[3];
rz(-2.2413018) q[3];
sx q[3];
rz(1.7830199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0428697) q[2];
sx q[2];
rz(-1.1770153) q[2];
sx q[2];
rz(-0.086611835) q[2];
rz(-3.0232271) q[3];
sx q[3];
rz(-1.5425073) q[3];
sx q[3];
rz(-0.043244403) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31556986) q[0];
sx q[0];
rz(-1.0181027) q[0];
sx q[0];
rz(-2.3636567) q[0];
rz(0.28620537) q[1];
sx q[1];
rz(-0.83120167) q[1];
sx q[1];
rz(1.3298159) q[1];
rz(-0.61997531) q[2];
sx q[2];
rz(-2.740553) q[2];
sx q[2];
rz(0.16410826) q[2];
rz(0.56747464) q[3];
sx q[3];
rz(-2.6226061) q[3];
sx q[3];
rz(-2.0593986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
