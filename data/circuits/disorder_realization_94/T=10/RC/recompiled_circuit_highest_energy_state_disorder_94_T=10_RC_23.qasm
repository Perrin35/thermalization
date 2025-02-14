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
rz(0.46410528) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(2.2321489) q[1];
sx q[1];
rz(7.5997054) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6911294) q[0];
sx q[0];
rz(-3.0678406) q[0];
sx q[0];
rz(2.6534257) q[0];
rz(-pi) q[1];
rz(-0.34061247) q[2];
sx q[2];
rz(-1.6180218) q[2];
sx q[2];
rz(-1.7318341) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2006372) q[1];
sx q[1];
rz(-0.66087729) q[1];
sx q[1];
rz(-0.4502859) q[1];
x q[2];
rz(1.2795696) q[3];
sx q[3];
rz(-1.816245) q[3];
sx q[3];
rz(-3.0994741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2509649) q[2];
sx q[2];
rz(-2.784412) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62946573) q[0];
sx q[0];
rz(-2.1194206) q[0];
sx q[0];
rz(-2.8566991) q[0];
rz(-0.051305436) q[1];
sx q[1];
rz(-0.76301328) q[1];
sx q[1];
rz(0.6368534) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1276217) q[0];
sx q[0];
rz(-2.4590965) q[0];
sx q[0];
rz(-2.3861814) q[0];
x q[1];
rz(0.35164386) q[2];
sx q[2];
rz(-1.571901) q[2];
sx q[2];
rz(-2.0948834) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6057558) q[1];
sx q[1];
rz(-2.9510289) q[1];
sx q[1];
rz(2.7298353) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15588197) q[3];
sx q[3];
rz(-0.71195275) q[3];
sx q[3];
rz(-0.32526325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9597943) q[2];
sx q[2];
rz(-0.85796285) q[2];
sx q[2];
rz(-2.174343) q[2];
rz(-0.26425427) q[3];
sx q[3];
rz(-1.5033009) q[3];
sx q[3];
rz(-1.2359469) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018175) q[0];
sx q[0];
rz(-1.5793261) q[0];
sx q[0];
rz(2.0215969) q[0];
rz(-2.6657875) q[1];
sx q[1];
rz(-2.0640524) q[1];
sx q[1];
rz(2.8376104) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7426152) q[0];
sx q[0];
rz(-2.1362059) q[0];
sx q[0];
rz(3.1088367) q[0];
x q[1];
rz(-0.39865785) q[2];
sx q[2];
rz(-0.2572607) q[2];
sx q[2];
rz(-0.21695732) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3545565) q[1];
sx q[1];
rz(-0.72595412) q[1];
sx q[1];
rz(1.7257084) q[1];
rz(-pi) q[2];
rz(1.7695565) q[3];
sx q[3];
rz(-0.89722465) q[3];
sx q[3];
rz(2.4186132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9756644) q[2];
sx q[2];
rz(-2.1037481) q[2];
sx q[2];
rz(2.4456639) q[2];
rz(1.1337229) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(0.55293647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021779) q[0];
sx q[0];
rz(-2.1372097) q[0];
sx q[0];
rz(-1.0473921) q[0];
rz(1.9678496) q[1];
sx q[1];
rz(-0.67432299) q[1];
sx q[1];
rz(1.5210927) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0516617) q[0];
sx q[0];
rz(-1.0703846) q[0];
sx q[0];
rz(1.7080281) q[0];
rz(-pi) q[1];
rz(1.1999454) q[2];
sx q[2];
rz(-2.4003138) q[2];
sx q[2];
rz(-2.0855057) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4075267) q[1];
sx q[1];
rz(-2.1868949) q[1];
sx q[1];
rz(-1.0259088) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4269104) q[3];
sx q[3];
rz(-0.41926786) q[3];
sx q[3];
rz(-1.4803606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.095470458) q[2];
sx q[2];
rz(-2.0746524) q[2];
sx q[2];
rz(-0.53786892) q[2];
rz(1.6543903) q[3];
sx q[3];
rz(-0.40707773) q[3];
sx q[3];
rz(-0.68638221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122413) q[0];
sx q[0];
rz(-0.23433267) q[0];
sx q[0];
rz(-2.5382163) q[0];
rz(0.76250184) q[1];
sx q[1];
rz(-1.9489894) q[1];
sx q[1];
rz(2.4286483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.514587) q[0];
sx q[0];
rz(-2.7381185) q[0];
sx q[0];
rz(1.4019985) q[0];
x q[1];
rz(1.1080526) q[2];
sx q[2];
rz(-1.6959091) q[2];
sx q[2];
rz(0.81415983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34316396) q[1];
sx q[1];
rz(-2.6021932) q[1];
sx q[1];
rz(-1.1718962) q[1];
rz(-0.98865786) q[3];
sx q[3];
rz(-2.1063559) q[3];
sx q[3];
rz(-2.5861507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1684299) q[2];
sx q[2];
rz(-0.91732401) q[2];
sx q[2];
rz(1.6707576) q[2];
rz(0.51703185) q[3];
sx q[3];
rz(-1.0443338) q[3];
sx q[3];
rz(-1.8487336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1626728) q[0];
sx q[0];
rz(-0.51386583) q[0];
sx q[0];
rz(-2.5901219) q[0];
rz(0.035471352) q[1];
sx q[1];
rz(-1.1330117) q[1];
sx q[1];
rz(1.6112526) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96824232) q[0];
sx q[0];
rz(-0.98778546) q[0];
sx q[0];
rz(0.75847404) q[0];
rz(0.82357652) q[2];
sx q[2];
rz(-1.5334497) q[2];
sx q[2];
rz(0.74451358) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0083654005) q[1];
sx q[1];
rz(-2.3043345) q[1];
sx q[1];
rz(-1.0788697) q[1];
x q[2];
rz(2.2238268) q[3];
sx q[3];
rz(-1.6380596) q[3];
sx q[3];
rz(-1.3239087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0643206) q[2];
sx q[2];
rz(-0.52017009) q[2];
sx q[2];
rz(1.8989829) q[2];
rz(-0.51314917) q[3];
sx q[3];
rz(-0.41392252) q[3];
sx q[3];
rz(-1.7665524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4921017) q[0];
sx q[0];
rz(-0.92925564) q[0];
sx q[0];
rz(0.11665601) q[0];
rz(2.3684582) q[1];
sx q[1];
rz(-1.9832858) q[1];
sx q[1];
rz(-1.6366417) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3706544) q[0];
sx q[0];
rz(-0.8621093) q[0];
sx q[0];
rz(-0.38972008) q[0];
rz(-pi) q[1];
rz(0.059594056) q[2];
sx q[2];
rz(-1.182654) q[2];
sx q[2];
rz(1.7517032) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9147787) q[1];
sx q[1];
rz(-0.89229167) q[1];
sx q[1];
rz(-0.46636502) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7724593) q[3];
sx q[3];
rz(-0.68644409) q[3];
sx q[3];
rz(0.74147195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2152805) q[2];
sx q[2];
rz(-0.24682385) q[2];
sx q[2];
rz(0.8872633) q[2];
rz(2.4617646) q[3];
sx q[3];
rz(-2.4726548) q[3];
sx q[3];
rz(0.63290709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585084) q[0];
sx q[0];
rz(-1.8483138) q[0];
sx q[0];
rz(0.65628091) q[0];
rz(0.20690021) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(1.2178749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.525187) q[0];
sx q[0];
rz(-1.6058679) q[0];
sx q[0];
rz(-0.082790815) q[0];
rz(2.0711259) q[2];
sx q[2];
rz(-0.64979759) q[2];
sx q[2];
rz(1.4184679) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6177442) q[1];
sx q[1];
rz(-2.9459125) q[1];
sx q[1];
rz(2.5168672) q[1];
rz(-2.5303305) q[3];
sx q[3];
rz(-1.3484133) q[3];
sx q[3];
rz(-2.5739447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3394341) q[2];
sx q[2];
rz(-1.3863486) q[2];
sx q[2];
rz(0.68186861) q[2];
rz(-1.1784461) q[3];
sx q[3];
rz(-2.384187) q[3];
sx q[3];
rz(1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39663974) q[0];
sx q[0];
rz(-1.5322026) q[0];
sx q[0];
rz(-0.21016453) q[0];
rz(0.22178966) q[1];
sx q[1];
rz(-0.91786018) q[1];
sx q[1];
rz(0.04235696) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51084679) q[0];
sx q[0];
rz(-1.5928942) q[0];
sx q[0];
rz(0.18760292) q[0];
rz(2.1788939) q[2];
sx q[2];
rz(-2.1129932) q[2];
sx q[2];
rz(-1.0350943) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48407979) q[1];
sx q[1];
rz(-1.8575605) q[1];
sx q[1];
rz(-2.5628255) q[1];
x q[2];
rz(-0.5333535) q[3];
sx q[3];
rz(-1.9755529) q[3];
sx q[3];
rz(-0.08431708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.28045851) q[2];
sx q[2];
rz(-2.0529842) q[2];
sx q[2];
rz(2.4750278) q[2];
rz(-2.376453) q[3];
sx q[3];
rz(-2.9355526) q[3];
sx q[3];
rz(-1.4008745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53369451) q[0];
sx q[0];
rz(-0.38132897) q[0];
sx q[0];
rz(0.68699849) q[0];
rz(-1.8580565) q[1];
sx q[1];
rz(-1.1524408) q[1];
sx q[1];
rz(-2.5680465) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6476558) q[0];
sx q[0];
rz(-2.082061) q[0];
sx q[0];
rz(2.4535116) q[0];
x q[1];
rz(1.1186662) q[2];
sx q[2];
rz(-1.3616478) q[2];
sx q[2];
rz(2.6191408) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0189218) q[1];
sx q[1];
rz(-1.9098567) q[1];
sx q[1];
rz(1.0215205) q[1];
rz(-pi) q[2];
rz(-0.76400842) q[3];
sx q[3];
rz(-2.2413018) q[3];
sx q[3];
rz(1.3585728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0428697) q[2];
sx q[2];
rz(-1.9645773) q[2];
sx q[2];
rz(-3.0549808) q[2];
rz(-0.11836554) q[3];
sx q[3];
rz(-1.5990853) q[3];
sx q[3];
rz(-0.043244403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260228) q[0];
sx q[0];
rz(-2.12349) q[0];
sx q[0];
rz(0.77793599) q[0];
rz(-2.8553873) q[1];
sx q[1];
rz(-0.83120167) q[1];
sx q[1];
rz(1.3298159) q[1];
rz(0.61997531) q[2];
sx q[2];
rz(-0.4010396) q[2];
sx q[2];
rz(-2.9774844) q[2];
rz(1.2729011) q[3];
sx q[3];
rz(-2.0023228) q[3];
sx q[3];
rz(1.7154233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
