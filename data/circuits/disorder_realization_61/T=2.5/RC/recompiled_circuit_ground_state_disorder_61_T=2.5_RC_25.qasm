OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1516079) q[0];
sx q[0];
rz(-0.94010544) q[0];
sx q[0];
rz(0.54036933) q[0];
rz(-2.6481533) q[1];
sx q[1];
rz(2.4210338) q[1];
sx q[1];
rz(8.8109206) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85640192) q[0];
sx q[0];
rz(-2.4840691) q[0];
sx q[0];
rz(2.8372) q[0];
rz(-pi) q[1];
rz(-3.1356642) q[2];
sx q[2];
rz(-1.8625755) q[2];
sx q[2];
rz(-0.31729441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39345523) q[1];
sx q[1];
rz(-0.84987133) q[1];
sx q[1];
rz(1.9490521) q[1];
rz(-pi) q[2];
rz(-0.085569445) q[3];
sx q[3];
rz(-0.24805476) q[3];
sx q[3];
rz(2.5906467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8453688) q[2];
sx q[2];
rz(-1.8989398) q[2];
sx q[2];
rz(-2.5073012) q[2];
rz(0.93572179) q[3];
sx q[3];
rz(-0.24565419) q[3];
sx q[3];
rz(-0.74265695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3212386) q[0];
sx q[0];
rz(-1.6004434) q[0];
sx q[0];
rz(0.4441922) q[0];
rz(-0.76002899) q[1];
sx q[1];
rz(-2.0030231) q[1];
sx q[1];
rz(-0.98145032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4430005) q[0];
sx q[0];
rz(-2.5410497) q[0];
sx q[0];
rz(-0.53437676) q[0];
rz(1.74238) q[2];
sx q[2];
rz(-1.7497471) q[2];
sx q[2];
rz(1.6323665) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3722575) q[1];
sx q[1];
rz(-0.44722873) q[1];
sx q[1];
rz(2.9370537) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6098374) q[3];
sx q[3];
rz(-0.49177836) q[3];
sx q[3];
rz(-0.39612202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0280219) q[2];
sx q[2];
rz(-2.5271723) q[2];
sx q[2];
rz(-1.9363972) q[2];
rz(-0.53660721) q[3];
sx q[3];
rz(-1.9034932) q[3];
sx q[3];
rz(-1.4046148) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2236915) q[0];
sx q[0];
rz(-1.8603928) q[0];
sx q[0];
rz(2.3028288) q[0];
rz(-2.5054848) q[1];
sx q[1];
rz(-1.4898224) q[1];
sx q[1];
rz(0.020523358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3394525) q[0];
sx q[0];
rz(-2.1198556) q[0];
sx q[0];
rz(-2.2744176) q[0];
rz(-pi) q[1];
rz(-2.5319935) q[2];
sx q[2];
rz(-0.72270279) q[2];
sx q[2];
rz(2.056207) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5387303) q[1];
sx q[1];
rz(-2.7285721) q[1];
sx q[1];
rz(-0.018750359) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2853363) q[3];
sx q[3];
rz(-1.4196287) q[3];
sx q[3];
rz(2.5976439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3402349) q[2];
sx q[2];
rz(-1.5284208) q[2];
sx q[2];
rz(2.197263) q[2];
rz(-1.0229735) q[3];
sx q[3];
rz(-1.9187656) q[3];
sx q[3];
rz(-2.003722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706547) q[0];
sx q[0];
rz(-1.174467) q[0];
sx q[0];
rz(-1.1573855) q[0];
rz(-1.4986787) q[1];
sx q[1];
rz(-1.4721556) q[1];
sx q[1];
rz(-1.9532983) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76283264) q[0];
sx q[0];
rz(-2.4597557) q[0];
sx q[0];
rz(2.3781611) q[0];
x q[1];
rz(-0.073578667) q[2];
sx q[2];
rz(-2.5917604) q[2];
sx q[2];
rz(0.64752196) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77515652) q[1];
sx q[1];
rz(-0.85298733) q[1];
sx q[1];
rz(0.78182533) q[1];
rz(-2.0075624) q[3];
sx q[3];
rz(-2.7705857) q[3];
sx q[3];
rz(-1.3001668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1184065) q[2];
sx q[2];
rz(-1.964317) q[2];
sx q[2];
rz(-2.4510621) q[2];
rz(-2.0265419) q[3];
sx q[3];
rz(-2.3670022) q[3];
sx q[3];
rz(1.5007277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0951776) q[0];
sx q[0];
rz(-1.7061808) q[0];
sx q[0];
rz(2.114356) q[0];
rz(-1.8401624) q[1];
sx q[1];
rz(-2.4899028) q[1];
sx q[1];
rz(2.8177736) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3394748) q[0];
sx q[0];
rz(-0.81011745) q[0];
sx q[0];
rz(1.0056061) q[0];
x q[1];
rz(0.01417966) q[2];
sx q[2];
rz(-2.1235596) q[2];
sx q[2];
rz(-1.9270037) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4869838) q[1];
sx q[1];
rz(-1.8909847) q[1];
sx q[1];
rz(0.57211188) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15669723) q[3];
sx q[3];
rz(-1.0644056) q[3];
sx q[3];
rz(1.2676257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7258437) q[2];
sx q[2];
rz(-2.9314633) q[2];
sx q[2];
rz(2.6591163) q[2];
rz(-1.1680565) q[3];
sx q[3];
rz(-1.9246293) q[3];
sx q[3];
rz(-2.4647958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2060858) q[0];
sx q[0];
rz(-0.85010234) q[0];
sx q[0];
rz(-2.2555943) q[0];
rz(0.63367263) q[1];
sx q[1];
rz(-0.88992563) q[1];
sx q[1];
rz(0.69127965) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38884866) q[0];
sx q[0];
rz(-0.92557478) q[0];
sx q[0];
rz(-1.6165125) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0646661) q[2];
sx q[2];
rz(-2.4483213) q[2];
sx q[2];
rz(-1.537854) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.012612494) q[1];
sx q[1];
rz(-0.75769934) q[1];
sx q[1];
rz(-1.6974259) q[1];
rz(-2.8915358) q[3];
sx q[3];
rz(-2.3365006) q[3];
sx q[3];
rz(-0.4522194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1401691) q[2];
sx q[2];
rz(-0.836335) q[2];
sx q[2];
rz(0.26958618) q[2];
rz(-0.94830281) q[3];
sx q[3];
rz(-1.6092665) q[3];
sx q[3];
rz(1.74291) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7798994) q[0];
sx q[0];
rz(-0.43710709) q[0];
sx q[0];
rz(1.0070739) q[0];
rz(-2.7359447) q[1];
sx q[1];
rz(-2.5463153) q[1];
sx q[1];
rz(-0.55353037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2562228) q[0];
sx q[0];
rz(-1.6281343) q[0];
sx q[0];
rz(-1.8256515) q[0];
rz(-pi) q[1];
rz(0.78041665) q[2];
sx q[2];
rz(-2.4274106) q[2];
sx q[2];
rz(1.2250021) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1055401) q[1];
sx q[1];
rz(-0.80784384) q[1];
sx q[1];
rz(2.4921662) q[1];
rz(-pi) q[2];
rz(1.5958435) q[3];
sx q[3];
rz(-1.1908997) q[3];
sx q[3];
rz(0.27835007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3125399) q[2];
sx q[2];
rz(-1.092814) q[2];
sx q[2];
rz(-1.2139758) q[2];
rz(-2.35516) q[3];
sx q[3];
rz(-1.395547) q[3];
sx q[3];
rz(3.1184375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.9455652) q[0];
sx q[0];
rz(-0.29226154) q[0];
sx q[0];
rz(-1.8096402) q[0];
rz(0.61344433) q[1];
sx q[1];
rz(-1.0204126) q[1];
sx q[1];
rz(2.6920998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17380781) q[0];
sx q[0];
rz(-1.0840084) q[0];
sx q[0];
rz(-3.0474902) q[0];
x q[1];
rz(-1.9067326) q[2];
sx q[2];
rz(-1.6523696) q[2];
sx q[2];
rz(2.5334266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66287884) q[1];
sx q[1];
rz(-0.7906853) q[1];
sx q[1];
rz(1.3776758) q[1];
rz(-pi) q[2];
rz(-0.2576377) q[3];
sx q[3];
rz(-1.6272568) q[3];
sx q[3];
rz(-0.27654058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8248262) q[2];
sx q[2];
rz(-0.68244857) q[2];
sx q[2];
rz(0.60834926) q[2];
rz(1.1634722) q[3];
sx q[3];
rz(-1.7721662) q[3];
sx q[3];
rz(-0.56330645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0332396) q[0];
sx q[0];
rz(-2.2305363) q[0];
sx q[0];
rz(-0.60229993) q[0];
rz(-2.1741518) q[1];
sx q[1];
rz(-0.93092218) q[1];
sx q[1];
rz(2.0379351) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7491584) q[0];
sx q[0];
rz(-1.1811678) q[0];
sx q[0];
rz(0.97831877) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5854563) q[2];
sx q[2];
rz(-0.93428627) q[2];
sx q[2];
rz(0.16600641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1608354) q[1];
sx q[1];
rz(-0.45113647) q[1];
sx q[1];
rz(-1.5680997) q[1];
rz(-pi) q[2];
rz(2.1488701) q[3];
sx q[3];
rz(-1.5141271) q[3];
sx q[3];
rz(-2.4538159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2063107) q[2];
sx q[2];
rz(-1.3056359) q[2];
sx q[2];
rz(-0.3375816) q[2];
rz(2.857699) q[3];
sx q[3];
rz(-2.4604535) q[3];
sx q[3];
rz(0.18686992) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5934481) q[0];
sx q[0];
rz(-1.633506) q[0];
sx q[0];
rz(-3.0737851) q[0];
rz(-1.1495122) q[1];
sx q[1];
rz(-1.5905453) q[1];
sx q[1];
rz(-2.6670719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.041458) q[0];
sx q[0];
rz(-2.9024419) q[0];
sx q[0];
rz(1.6275703) q[0];
rz(-pi) q[1];
rz(0.50751792) q[2];
sx q[2];
rz(-2.8675277) q[2];
sx q[2];
rz(3.0176891) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.671512) q[1];
sx q[1];
rz(-1.5776411) q[1];
sx q[1];
rz(1.4840625) q[1];
rz(-pi) q[2];
rz(-0.69863221) q[3];
sx q[3];
rz(-0.94480522) q[3];
sx q[3];
rz(2.0076942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6600251) q[2];
sx q[2];
rz(-2.2109172) q[2];
sx q[2];
rz(2.068326) q[2];
rz(-2.977071) q[3];
sx q[3];
rz(-1.3273032) q[3];
sx q[3];
rz(2.8209414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5545223) q[0];
sx q[0];
rz(-1.5656492) q[0];
sx q[0];
rz(1.5026305) q[0];
rz(1.5837689) q[1];
sx q[1];
rz(-1.0777892) q[1];
sx q[1];
rz(-0.52660175) q[1];
rz(-0.64939349) q[2];
sx q[2];
rz(-0.55772256) q[2];
sx q[2];
rz(-1.100308) q[2];
rz(1.6526374) q[3];
sx q[3];
rz(-1.428953) q[3];
sx q[3];
rz(-0.088464213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
