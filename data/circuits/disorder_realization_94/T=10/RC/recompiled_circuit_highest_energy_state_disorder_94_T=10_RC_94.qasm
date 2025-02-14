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
rz(-1.8985768) q[0];
sx q[0];
rz(-2.0199116) q[0];
sx q[0];
rz(-0.46410528) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(-0.90944374) q[1];
sx q[1];
rz(-1.3165201) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45046321) q[0];
sx q[0];
rz(-0.073752068) q[0];
sx q[0];
rz(-2.6534257) q[0];
rz(-pi) q[1];
rz(-2.8009802) q[2];
sx q[2];
rz(-1.6180218) q[2];
sx q[2];
rz(1.7318341) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7500748) q[1];
sx q[1];
rz(-0.98528359) q[1];
sx q[1];
rz(1.8970918) q[1];
rz(-pi) q[2];
rz(-1.2795696) q[3];
sx q[3];
rz(-1.3253477) q[3];
sx q[3];
rz(-3.0994741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89062771) q[2];
sx q[2];
rz(-2.784412) q[2];
sx q[2];
rz(0.29699057) q[2];
rz(-2.4885528) q[3];
sx q[3];
rz(-1.5831455) q[3];
sx q[3];
rz(-1.0562586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62946573) q[0];
sx q[0];
rz(-2.1194206) q[0];
sx q[0];
rz(-2.8566991) q[0];
rz(-3.0902872) q[1];
sx q[1];
rz(-0.76301328) q[1];
sx q[1];
rz(2.5047393) q[1];
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
rz(1.5696196) q[2];
sx q[2];
rz(-1.92244) q[2];
sx q[2];
rz(-0.52368173) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11733774) q[1];
sx q[1];
rz(-1.7452612) q[1];
sx q[1];
rz(-1.4937449) q[1];
rz(-2.4356682) q[3];
sx q[3];
rz(-1.6723989) q[3];
sx q[3];
rz(1.3639579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9597943) q[2];
sx q[2];
rz(-0.85796285) q[2];
sx q[2];
rz(-0.96724969) q[2];
rz(2.8773384) q[3];
sx q[3];
rz(-1.6382917) q[3];
sx q[3];
rz(-1.9056457) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13977519) q[0];
sx q[0];
rz(-1.5793261) q[0];
sx q[0];
rz(1.1199957) q[0];
rz(2.6657875) q[1];
sx q[1];
rz(-1.0775403) q[1];
sx q[1];
rz(2.8376104) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7426152) q[0];
sx q[0];
rz(-1.0053867) q[0];
sx q[0];
rz(-0.032755927) q[0];
x q[1];
rz(0.39865785) q[2];
sx q[2];
rz(-2.884332) q[2];
sx q[2];
rz(-0.21695732) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3545565) q[1];
sx q[1];
rz(-2.4156385) q[1];
sx q[1];
rz(-1.7257084) q[1];
rz(2.899053) q[3];
sx q[3];
rz(-0.69787301) q[3];
sx q[3];
rz(-1.0353116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9756644) q[2];
sx q[2];
rz(-1.0378446) q[2];
sx q[2];
rz(0.69592875) q[2];
rz(-2.0078697) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(-2.5886562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021779) q[0];
sx q[0];
rz(-1.004383) q[0];
sx q[0];
rz(-2.0942005) q[0];
rz(-1.9678496) q[1];
sx q[1];
rz(-0.67432299) q[1];
sx q[1];
rz(1.6204999) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3701909) q[0];
sx q[0];
rz(-2.6242497) q[0];
sx q[0];
rz(0.24513729) q[0];
rz(-pi) q[1];
rz(-2.8212566) q[2];
sx q[2];
rz(-2.2514859) q[2];
sx q[2];
rz(-2.5706511) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17361372) q[1];
sx q[1];
rz(-1.1340177) q[1];
sx q[1];
rz(-0.69154214) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8168746) q[3];
sx q[3];
rz(-1.8408662) q[3];
sx q[3];
rz(-2.3809759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0461222) q[2];
sx q[2];
rz(-1.0669402) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2122413) q[0];
sx q[0];
rz(-2.90726) q[0];
sx q[0];
rz(-2.5382163) q[0];
rz(0.76250184) q[1];
sx q[1];
rz(-1.9489894) q[1];
sx q[1];
rz(2.4286483) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6978074) q[0];
sx q[0];
rz(-1.968211) q[0];
sx q[0];
rz(0.071594588) q[0];
rz(-pi) q[1];
x q[1];
rz(3.001956) q[2];
sx q[2];
rz(-1.1119482) q[2];
sx q[2];
rz(2.322784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11345574) q[1];
sx q[1];
rz(-1.077768) q[1];
sx q[1];
rz(-2.9131469) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6176881) q[3];
sx q[3];
rz(-2.0632944) q[3];
sx q[3];
rz(0.69128752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1684299) q[2];
sx q[2];
rz(-2.2242686) q[2];
sx q[2];
rz(1.470835) q[2];
rz(-0.51703185) q[3];
sx q[3];
rz(-2.0972589) q[3];
sx q[3];
rz(1.2928591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1626728) q[0];
sx q[0];
rz(-0.51386583) q[0];
sx q[0];
rz(-2.5901219) q[0];
rz(3.1061213) q[1];
sx q[1];
rz(-1.1330117) q[1];
sx q[1];
rz(-1.6112526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0198675) q[0];
sx q[0];
rz(-2.1824153) q[0];
sx q[0];
rz(-2.3082971) q[0];
x q[1];
rz(0.82357652) q[2];
sx q[2];
rz(-1.5334497) q[2];
sx q[2];
rz(0.74451358) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0083654005) q[1];
sx q[1];
rz(-0.83725819) q[1];
sx q[1];
rz(1.0788697) q[1];
x q[2];
rz(3.0569791) q[3];
sx q[3];
rz(-0.91949465) q[3];
sx q[3];
rz(-0.29825975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0643206) q[2];
sx q[2];
rz(-0.52017009) q[2];
sx q[2];
rz(-1.8989829) q[2];
rz(-2.6284435) q[3];
sx q[3];
rz(-2.7276701) q[3];
sx q[3];
rz(-1.7665524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-3.0249366) q[0];
rz(-0.77313441) q[1];
sx q[1];
rz(-1.1583068) q[1];
sx q[1];
rz(1.6366417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0610675) q[0];
sx q[0];
rz(-1.8634029) q[0];
sx q[0];
rz(2.3181897) q[0];
x q[1];
rz(0.059594056) q[2];
sx q[2];
rz(-1.9589387) q[2];
sx q[2];
rz(-1.7517032) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6924062) q[1];
sx q[1];
rz(-0.80185651) q[1];
sx q[1];
rz(1.0620326) q[1];
rz(2.4890763) q[3];
sx q[3];
rz(-1.8015141) q[3];
sx q[3];
rz(-1.1200865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2152805) q[2];
sx q[2];
rz(-2.8947688) q[2];
sx q[2];
rz(-0.8872633) q[2];
rz(-0.67982802) q[3];
sx q[3];
rz(-0.66893783) q[3];
sx q[3];
rz(-0.63290709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585084) q[0];
sx q[0];
rz(-1.2932788) q[0];
sx q[0];
rz(2.4853117) q[0];
rz(0.20690021) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(-1.9237178) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6961795) q[0];
sx q[0];
rz(-0.089897297) q[0];
sx q[0];
rz(-0.40125664) q[0];
rz(0.98274173) q[2];
sx q[2];
rz(-1.2763192) q[2];
sx q[2];
rz(0.5628995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52384842) q[1];
sx q[1];
rz(-2.9459125) q[1];
sx q[1];
rz(0.62472549) q[1];
rz(0.61126216) q[3];
sx q[3];
rz(-1.3484133) q[3];
sx q[3];
rz(-2.5739447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8021585) q[2];
sx q[2];
rz(-1.755244) q[2];
sx q[2];
rz(0.68186861) q[2];
rz(1.9631466) q[3];
sx q[3];
rz(-2.384187) q[3];
sx q[3];
rz(1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39663974) q[0];
sx q[0];
rz(-1.5322026) q[0];
sx q[0];
rz(-2.9314281) q[0];
rz(-0.22178966) q[1];
sx q[1];
rz(-2.2237325) q[1];
sx q[1];
rz(-3.0992357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9440749) q[0];
sx q[0];
rz(-2.952708) q[0];
sx q[0];
rz(3.023639) q[0];
x q[1];
rz(2.1788939) q[2];
sx q[2];
rz(-2.1129932) q[2];
sx q[2];
rz(-1.0350943) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6463464) q[1];
sx q[1];
rz(-2.5030285) q[1];
sx q[1];
rz(2.6471443) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6082392) q[3];
sx q[3];
rz(-1.9755529) q[3];
sx q[3];
rz(-0.08431708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28045851) q[2];
sx q[2];
rz(-1.0886085) q[2];
sx q[2];
rz(0.66656485) q[2];
rz(-2.376453) q[3];
sx q[3];
rz(-0.20604006) q[3];
sx q[3];
rz(1.4008745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6078981) q[0];
sx q[0];
rz(-2.7602637) q[0];
sx q[0];
rz(0.68699849) q[0];
rz(-1.8580565) q[1];
sx q[1];
rz(-1.9891519) q[1];
sx q[1];
rz(2.5680465) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6015778) q[0];
sx q[0];
rz(-0.83160831) q[0];
sx q[0];
rz(2.4180146) q[0];
x q[1];
rz(-2.0229265) q[2];
sx q[2];
rz(-1.3616478) q[2];
sx q[2];
rz(2.6191408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0189218) q[1];
sx q[1];
rz(-1.9098567) q[1];
sx q[1];
rz(2.1200722) q[1];
x q[2];
rz(-2.403026) q[3];
sx q[3];
rz(-2.1436678) q[3];
sx q[3];
rz(0.74921879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0428697) q[2];
sx q[2];
rz(-1.9645773) q[2];
sx q[2];
rz(0.086611835) q[2];
rz(0.11836554) q[3];
sx q[3];
rz(-1.5990853) q[3];
sx q[3];
rz(0.043244403) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260228) q[0];
sx q[0];
rz(-2.12349) q[0];
sx q[0];
rz(0.77793599) q[0];
rz(-0.28620537) q[1];
sx q[1];
rz(-2.310391) q[1];
sx q[1];
rz(-1.8117767) q[1];
rz(-2.5216173) q[2];
sx q[2];
rz(-0.4010396) q[2];
sx q[2];
rz(-2.9774844) q[2];
rz(2.6927041) q[3];
sx q[3];
rz(-1.3009303) q[3];
sx q[3];
rz(-3.124685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
