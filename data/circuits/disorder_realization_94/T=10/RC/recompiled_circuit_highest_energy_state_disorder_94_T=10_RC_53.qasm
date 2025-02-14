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
rz(2.0199116) q[0];
sx q[0];
rz(12.102265) q[0];
rz(-0.75587455) q[1];
sx q[1];
rz(-2.2321489) q[1];
sx q[1];
rz(1.3165201) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038832207) q[0];
sx q[0];
rz(-1.6359207) q[0];
sx q[0];
rz(-1.6054356) q[0];
rz(-pi) q[1];
rz(3.0010536) q[2];
sx q[2];
rz(-2.7978483) q[2];
sx q[2];
rz(3.112971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7500748) q[1];
sx q[1];
rz(-2.1563091) q[1];
sx q[1];
rz(1.2445009) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2795696) q[3];
sx q[3];
rz(-1.816245) q[3];
sx q[3];
rz(0.042118532) q[3];
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
rz(0.29699057) q[2];
rz(0.65303981) q[3];
sx q[3];
rz(-1.5831455) q[3];
sx q[3];
rz(-1.0562586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5121269) q[0];
sx q[0];
rz(-1.0221721) q[0];
sx q[0];
rz(2.8566991) q[0];
rz(3.0902872) q[1];
sx q[1];
rz(-2.3785794) q[1];
sx q[1];
rz(-0.6368534) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013971) q[0];
sx q[0];
rz(-0.68249615) q[0];
sx q[0];
rz(-0.75541124) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0032071438) q[2];
sx q[2];
rz(-0.35164552) q[2];
sx q[2];
rz(2.6144947) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4400582) q[1];
sx q[1];
rz(-1.4949168) q[1];
sx q[1];
rz(-0.17497356) q[1];
rz(-pi) q[2];
rz(-2.4356682) q[3];
sx q[3];
rz(-1.6723989) q[3];
sx q[3];
rz(1.3639579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.18179831) q[2];
sx q[2];
rz(-0.85796285) q[2];
sx q[2];
rz(-0.96724969) q[2];
rz(2.8773384) q[3];
sx q[3];
rz(-1.5033009) q[3];
sx q[3];
rz(1.9056457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13977519) q[0];
sx q[0];
rz(-1.5622666) q[0];
sx q[0];
rz(-1.1199957) q[0];
rz(-0.4758052) q[1];
sx q[1];
rz(-2.0640524) q[1];
sx q[1];
rz(-2.8376104) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39897746) q[0];
sx q[0];
rz(-2.1362059) q[0];
sx q[0];
rz(0.032755927) q[0];
rz(-1.6725704) q[2];
sx q[2];
rz(-1.3341122) q[2];
sx q[2];
rz(2.9477811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5811823) q[1];
sx q[1];
rz(-0.85542233) q[1];
sx q[1];
rz(0.13611273) q[1];
x q[2];
rz(-1.3720361) q[3];
sx q[3];
rz(-2.244368) q[3];
sx q[3];
rz(-2.4186132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1659282) q[2];
sx q[2];
rz(-1.0378446) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43941471) q[0];
sx q[0];
rz(-1.004383) q[0];
sx q[0];
rz(-1.0473921) q[0];
rz(-1.1737431) q[1];
sx q[1];
rz(-0.67432299) q[1];
sx q[1];
rz(1.5210927) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0516617) q[0];
sx q[0];
rz(-2.0712081) q[0];
sx q[0];
rz(1.4335645) q[0];
x q[1];
rz(0.86444433) q[2];
sx q[2];
rz(-1.3235759) q[2];
sx q[2];
rz(0.79402393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4075267) q[1];
sx q[1];
rz(-2.1868949) q[1];
sx q[1];
rz(2.1156838) q[1];
rz(-pi) q[2];
rz(-1.8549881) q[3];
sx q[3];
rz(-1.8833369) q[3];
sx q[3];
rz(-0.72060637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0461222) q[2];
sx q[2];
rz(-2.0746524) q[2];
sx q[2];
rz(2.6037237) q[2];
rz(1.6543903) q[3];
sx q[3];
rz(-2.7345149) q[3];
sx q[3];
rz(-2.4552104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(1.2122413) q[0];
sx q[0];
rz(-0.23433267) q[0];
sx q[0];
rz(0.60337639) q[0];
rz(0.76250184) q[1];
sx q[1];
rz(-1.9489894) q[1];
sx q[1];
rz(-0.71294436) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099261053) q[0];
sx q[0];
rz(-1.6368027) q[0];
sx q[0];
rz(1.1724654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.845417) q[2];
sx q[2];
rz(-2.6634187) q[2];
sx q[2];
rz(-0.51152767) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0281369) q[1];
sx q[1];
rz(-1.077768) q[1];
sx q[1];
rz(-2.9131469) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3943123) q[3];
sx q[3];
rz(-2.3722014) q[3];
sx q[3];
rz(1.6748493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1684299) q[2];
sx q[2];
rz(-2.2242686) q[2];
sx q[2];
rz(-1.6707576) q[2];
rz(2.6245608) q[3];
sx q[3];
rz(-2.0972589) q[3];
sx q[3];
rz(-1.8487336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9789199) q[0];
sx q[0];
rz(-2.6277268) q[0];
sx q[0];
rz(-0.55147076) q[0];
rz(-3.1061213) q[1];
sx q[1];
rz(-2.008581) q[1];
sx q[1];
rz(-1.6112526) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96824232) q[0];
sx q[0];
rz(-2.1538072) q[0];
sx q[0];
rz(-0.75847404) q[0];
rz(-1.5158723) q[2];
sx q[2];
rz(-0.74797219) q[2];
sx q[2];
rz(0.8665646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.0083654005) q[1];
sx q[1];
rz(-0.83725819) q[1];
sx q[1];
rz(-2.0627229) q[1];
rz(-pi) q[2];
rz(-0.084613581) q[3];
sx q[3];
rz(-0.91949465) q[3];
sx q[3];
rz(2.8433329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0772721) q[2];
sx q[2];
rz(-2.6214226) q[2];
sx q[2];
rz(-1.2426097) q[2];
rz(-2.6284435) q[3];
sx q[3];
rz(-0.41392252) q[3];
sx q[3];
rz(1.7665524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4921017) q[0];
sx q[0];
rz(-0.92925564) q[0];
sx q[0];
rz(-0.11665601) q[0];
rz(-0.77313441) q[1];
sx q[1];
rz(-1.9832858) q[1];
sx q[1];
rz(-1.6366417) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7709383) q[0];
sx q[0];
rz(-2.2794834) q[0];
sx q[0];
rz(-0.38972008) q[0];
rz(-3.0819986) q[2];
sx q[2];
rz(-1.9589387) q[2];
sx q[2];
rz(1.3898894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6924062) q[1];
sx q[1];
rz(-0.80185651) q[1];
sx q[1];
rz(1.0620326) q[1];
rz(-2.4890763) q[3];
sx q[3];
rz(-1.8015141) q[3];
sx q[3];
rz(-2.0215061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2152805) q[2];
sx q[2];
rz(-2.8947688) q[2];
sx q[2];
rz(0.8872633) q[2];
rz(0.67982802) q[3];
sx q[3];
rz(-0.66893783) q[3];
sx q[3];
rz(-2.5086856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585084) q[0];
sx q[0];
rz(-1.2932788) q[0];
sx q[0];
rz(2.4853117) q[0];
rz(2.9346924) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(-1.2178749) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.098893) q[0];
sx q[0];
rz(-1.4880565) q[0];
sx q[0];
rz(1.5356043) q[0];
rz(-pi) q[1];
rz(2.0711259) q[2];
sx q[2];
rz(-0.64979759) q[2];
sx q[2];
rz(1.4184679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7102573) q[1];
sx q[1];
rz(-1.4568304) q[1];
sx q[1];
rz(2.9821787) q[1];
rz(1.3013873) q[3];
sx q[3];
rz(-0.97668934) q[3];
sx q[3];
rz(1.9850933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8021585) q[2];
sx q[2];
rz(-1.3863486) q[2];
sx q[2];
rz(0.68186861) q[2];
rz(1.9631466) q[3];
sx q[3];
rz(-2.384187) q[3];
sx q[3];
rz(-1.2371548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39663974) q[0];
sx q[0];
rz(-1.6093901) q[0];
sx q[0];
rz(0.21016453) q[0];
rz(-2.919803) q[1];
sx q[1];
rz(-2.2237325) q[1];
sx q[1];
rz(-0.04235696) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1975178) q[0];
sx q[0];
rz(-2.952708) q[0];
sx q[0];
rz(0.11795363) q[0];
rz(-0.63318166) q[2];
sx q[2];
rz(-1.0594308) q[2];
sx q[2];
rz(-0.88054576) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90396229) q[1];
sx q[1];
rz(-1.018486) q[1];
sx q[1];
rz(1.9094853) q[1];
rz(-pi) q[2];
rz(-0.70019763) q[3];
sx q[3];
rz(-0.65749107) q[3];
sx q[3];
rz(-2.0746865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28045851) q[2];
sx q[2];
rz(-2.0529842) q[2];
sx q[2];
rz(0.66656485) q[2];
rz(-0.76513964) q[3];
sx q[3];
rz(-0.20604006) q[3];
sx q[3];
rz(1.7407181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53369451) q[0];
sx q[0];
rz(-2.7602637) q[0];
sx q[0];
rz(-2.4545942) q[0];
rz(-1.8580565) q[1];
sx q[1];
rz(-1.1524408) q[1];
sx q[1];
rz(-2.5680465) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.682293) q[0];
sx q[0];
rz(-2.1577765) q[0];
sx q[0];
rz(-0.94265509) q[0];
rz(-1.1186662) q[2];
sx q[2];
rz(-1.3616478) q[2];
sx q[2];
rz(-2.6191408) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4926238) q[1];
sx q[1];
rz(-2.0855806) q[1];
sx q[1];
rz(0.3920946) q[1];
rz(2.2881094) q[3];
sx q[3];
rz(-0.9694582) q[3];
sx q[3];
rz(2.7784612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0428697) q[2];
sx q[2];
rz(-1.9645773) q[2];
sx q[2];
rz(0.086611835) q[2];
rz(0.11836554) q[3];
sx q[3];
rz(-1.5425073) q[3];
sx q[3];
rz(-0.043244403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8260228) q[0];
sx q[0];
rz(-1.0181027) q[0];
sx q[0];
rz(-2.3636567) q[0];
rz(-2.8553873) q[1];
sx q[1];
rz(-0.83120167) q[1];
sx q[1];
rz(1.3298159) q[1];
rz(-1.8123476) q[2];
sx q[2];
rz(-1.2474682) q[2];
sx q[2];
rz(-0.4954485) q[2];
rz(-1.8686915) q[3];
sx q[3];
rz(-2.0023228) q[3];
sx q[3];
rz(1.7154233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
