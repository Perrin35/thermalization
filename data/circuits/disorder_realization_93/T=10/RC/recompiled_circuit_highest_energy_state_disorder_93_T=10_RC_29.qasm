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
rz(0.93208575) q[0];
sx q[0];
rz(4.1756364) q[0];
sx q[0];
rz(11.812165) q[0];
rz(1.7475313) q[1];
sx q[1];
rz(-0.59352195) q[1];
sx q[1];
rz(1.7041748) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93255723) q[0];
sx q[0];
rz(-0.13929312) q[0];
sx q[0];
rz(-0.4362696) q[0];
rz(-3.1033663) q[2];
sx q[2];
rz(-1.1817721) q[2];
sx q[2];
rz(-1.807275) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.61753786) q[1];
sx q[1];
rz(-1.0144985) q[1];
sx q[1];
rz(1.9281045) q[1];
rz(-1.1877443) q[3];
sx q[3];
rz(-1.446734) q[3];
sx q[3];
rz(-1.2318486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5022016) q[2];
sx q[2];
rz(-2.907967) q[2];
sx q[2];
rz(0.26546738) q[2];
rz(-1.4207077) q[3];
sx q[3];
rz(-1.9749494) q[3];
sx q[3];
rz(1.7336806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.079064) q[0];
sx q[0];
rz(-1.205235) q[0];
sx q[0];
rz(1.0519387) q[0];
rz(2.1417292) q[1];
sx q[1];
rz(-1.544416) q[1];
sx q[1];
rz(0.76905191) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6175943) q[0];
sx q[0];
rz(-0.72223653) q[0];
sx q[0];
rz(0.87498949) q[0];
x q[1];
rz(1.6300419) q[2];
sx q[2];
rz(-2.1476204) q[2];
sx q[2];
rz(1.249493) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6144784) q[1];
sx q[1];
rz(-1.1243781) q[1];
sx q[1];
rz(2.9698257) q[1];
rz(1.814428) q[3];
sx q[3];
rz(-1.9852178) q[3];
sx q[3];
rz(-0.72539402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.96180463) q[2];
sx q[2];
rz(-0.804681) q[2];
sx q[2];
rz(-1.556832) q[2];
rz(2.4368317) q[3];
sx q[3];
rz(-2.9289991) q[3];
sx q[3];
rz(-1.0889277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8162808) q[0];
sx q[0];
rz(-2.3820057) q[0];
sx q[0];
rz(-0.76667205) q[0];
rz(-2.7491772) q[1];
sx q[1];
rz(-2.2800443) q[1];
sx q[1];
rz(2.3426985) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63215374) q[0];
sx q[0];
rz(-1.5611794) q[0];
sx q[0];
rz(-2.0548178) q[0];
x q[1];
rz(2.8908752) q[2];
sx q[2];
rz(-1.483416) q[2];
sx q[2];
rz(-0.45044294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0392363) q[1];
sx q[1];
rz(-1.2627541) q[1];
sx q[1];
rz(-2.0935489) q[1];
rz(-0.17744448) q[3];
sx q[3];
rz(-1.6758306) q[3];
sx q[3];
rz(2.3326186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7400292) q[2];
sx q[2];
rz(-0.26569685) q[2];
sx q[2];
rz(-3.1166039) q[2];
rz(-1.3966903) q[3];
sx q[3];
rz(-1.9091505) q[3];
sx q[3];
rz(-0.11740824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(1.5806737) q[0];
sx q[0];
rz(-0.54045254) q[0];
sx q[0];
rz(2.8448291) q[0];
rz(-0.98776039) q[1];
sx q[1];
rz(-1.8901653) q[1];
sx q[1];
rz(-0.74772778) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4241705) q[0];
sx q[0];
rz(-1.352373) q[0];
sx q[0];
rz(-2.3481247) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.426366) q[2];
sx q[2];
rz(-1.919121) q[2];
sx q[2];
rz(2.0722318) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8750413) q[1];
sx q[1];
rz(-1.5822486) q[1];
sx q[1];
rz(2.7391866) q[1];
x q[2];
rz(-2.4676855) q[3];
sx q[3];
rz(-1.0847989) q[3];
sx q[3];
rz(2.7844519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1838386) q[2];
sx q[2];
rz(-0.87205333) q[2];
sx q[2];
rz(-0.5198861) q[2];
rz(2.1947491) q[3];
sx q[3];
rz(-1.1917453) q[3];
sx q[3];
rz(-0.051518353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2648322) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(1.9448036) q[0];
rz(-1.6668677) q[1];
sx q[1];
rz(-2.3366172) q[1];
sx q[1];
rz(-2.8428452) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0419436) q[0];
sx q[0];
rz(-1.22166) q[0];
sx q[0];
rz(0.70167244) q[0];
x q[1];
rz(-2.6514154) q[2];
sx q[2];
rz(-2.2818447) q[2];
sx q[2];
rz(-0.74029628) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1890002) q[1];
sx q[1];
rz(-0.77356427) q[1];
sx q[1];
rz(0.1632593) q[1];
x q[2];
rz(0.87033724) q[3];
sx q[3];
rz(-1.2559818) q[3];
sx q[3];
rz(0.39488068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8220736) q[2];
sx q[2];
rz(-2.3910429) q[2];
sx q[2];
rz(3.1263515) q[2];
rz(-3.0139253) q[3];
sx q[3];
rz(-0.36104194) q[3];
sx q[3];
rz(0.84967363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78758883) q[0];
sx q[0];
rz(-2.6241527) q[0];
sx q[0];
rz(-2.2678243) q[0];
rz(-1.9173701) q[1];
sx q[1];
rz(-2.1189549) q[1];
sx q[1];
rz(1.2194941) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6343459) q[0];
sx q[0];
rz(-1.1218881) q[0];
sx q[0];
rz(0.44621356) q[0];
rz(2.5415433) q[2];
sx q[2];
rz(-0.82227899) q[2];
sx q[2];
rz(-2.6278969) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.732033) q[1];
sx q[1];
rz(-2.8581706) q[1];
sx q[1];
rz(0.83614142) q[1];
rz(-0.55762724) q[3];
sx q[3];
rz(-0.67071299) q[3];
sx q[3];
rz(1.1190377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.445861) q[2];
sx q[2];
rz(-1.6709238) q[2];
sx q[2];
rz(-2.640558) q[2];
rz(-1.8028353) q[3];
sx q[3];
rz(-2.7300291) q[3];
sx q[3];
rz(-1.1494273) q[3];
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
rz(-2.5434791) q[0];
sx q[0];
rz(-1.9444332) q[0];
sx q[0];
rz(0.57394779) q[0];
rz(0.57169882) q[1];
sx q[1];
rz(-2.1015002) q[1];
sx q[1];
rz(-1.5572825) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2435574) q[0];
sx q[0];
rz(-1.9429617) q[0];
sx q[0];
rz(2.283147) q[0];
rz(-pi) q[1];
rz(1.1102067) q[2];
sx q[2];
rz(-0.71515036) q[2];
sx q[2];
rz(-2.2543224) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7148278) q[1];
sx q[1];
rz(-2.0457771) q[1];
sx q[1];
rz(0.91678859) q[1];
rz(-1.1787703) q[3];
sx q[3];
rz(-1.2268664) q[3];
sx q[3];
rz(2.1363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.456936) q[2];
sx q[2];
rz(-1.9921649) q[2];
sx q[2];
rz(-1.1535025) q[2];
rz(1.5270816) q[3];
sx q[3];
rz(-2.1187512) q[3];
sx q[3];
rz(-1.8089186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1077147) q[0];
sx q[0];
rz(-2.3623473) q[0];
sx q[0];
rz(-2.9039134) q[0];
rz(0.73860812) q[1];
sx q[1];
rz(-1.2146726) q[1];
sx q[1];
rz(3.1225263) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3349318) q[0];
sx q[0];
rz(-1.4981759) q[0];
sx q[0];
rz(3.1333092) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3735061) q[2];
sx q[2];
rz(-1.0407018) q[2];
sx q[2];
rz(-1.4370611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70266378) q[1];
sx q[1];
rz(-0.090320822) q[1];
sx q[1];
rz(-2.0144736) q[1];
rz(-pi) q[2];
rz(0.64225853) q[3];
sx q[3];
rz(-1.0805784) q[3];
sx q[3];
rz(2.1497906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1261403) q[2];
sx q[2];
rz(-1.3498053) q[2];
sx q[2];
rz(-2.9665185) q[2];
rz(-1.763688) q[3];
sx q[3];
rz(-0.62052369) q[3];
sx q[3];
rz(3.0231045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0835251) q[0];
sx q[0];
rz(-1.3578992) q[0];
sx q[0];
rz(-0.51496664) q[0];
rz(0.99634755) q[1];
sx q[1];
rz(-1.0641655) q[1];
sx q[1];
rz(1.6641585) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9290976) q[0];
sx q[0];
rz(-0.97508865) q[0];
sx q[0];
rz(-2.2736249) q[0];
x q[1];
rz(-1.5375877) q[2];
sx q[2];
rz(-2.1589212) q[2];
sx q[2];
rz(2.3606481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0465999) q[1];
sx q[1];
rz(-1.0674202) q[1];
sx q[1];
rz(0.51790389) q[1];
rz(-2.1636073) q[3];
sx q[3];
rz(-1.3372652) q[3];
sx q[3];
rz(1.9625036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0763863) q[2];
sx q[2];
rz(-1.4896769) q[2];
sx q[2];
rz(2.9822986) q[2];
rz(-1.4984891) q[3];
sx q[3];
rz(-0.67125932) q[3];
sx q[3];
rz(1.8026277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047073929) q[0];
sx q[0];
rz(-3.0975332) q[0];
sx q[0];
rz(2.2799168) q[0];
rz(-1.8623976) q[1];
sx q[1];
rz(-0.28293124) q[1];
sx q[1];
rz(0.18712015) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8409654) q[0];
sx q[0];
rz(-1.8442796) q[0];
sx q[0];
rz(-1.8086368) q[0];
rz(0.4942221) q[2];
sx q[2];
rz(-1.3307327) q[2];
sx q[2];
rz(-2.6721045) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71003733) q[1];
sx q[1];
rz(-1.3042684) q[1];
sx q[1];
rz(0.73898594) q[1];
x q[2];
rz(-2.7469603) q[3];
sx q[3];
rz(-1.2773716) q[3];
sx q[3];
rz(-2.54486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.44237915) q[2];
sx q[2];
rz(-2.1670161) q[2];
sx q[2];
rz(1.7787735) q[2];
rz(1.5022701) q[3];
sx q[3];
rz(-0.5539186) q[3];
sx q[3];
rz(-1.7207918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61158553) q[0];
sx q[0];
rz(-1.8066318) q[0];
sx q[0];
rz(-3.1351177) q[0];
rz(0.50637983) q[1];
sx q[1];
rz(-0.19229278) q[1];
sx q[1];
rz(-2.2813588) q[1];
rz(-1.2564332) q[2];
sx q[2];
rz(-2.6994575) q[2];
sx q[2];
rz(-0.14002249) q[2];
rz(-0.68001775) q[3];
sx q[3];
rz(-2.802805) q[3];
sx q[3];
rz(-1.203385) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
