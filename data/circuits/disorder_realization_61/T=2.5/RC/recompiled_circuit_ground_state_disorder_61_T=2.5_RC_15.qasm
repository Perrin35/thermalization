OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9899848) q[0];
sx q[0];
rz(4.0816981) q[0];
sx q[0];
rz(8.8844086) q[0];
rz(-2.6481533) q[1];
sx q[1];
rz(2.4210338) q[1];
sx q[1];
rz(8.8109206) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47068025) q[0];
sx q[0];
rz(-1.3865836) q[0];
sx q[0];
rz(-0.6349011) q[0];
x q[1];
rz(0.0059284728) q[2];
sx q[2];
rz(-1.8625755) q[2];
sx q[2];
rz(2.8242982) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2207823) q[1];
sx q[1];
rz(-1.8519004) q[1];
sx q[1];
rz(-0.75741655) q[1];
rz(-3.0560232) q[3];
sx q[3];
rz(-0.24805476) q[3];
sx q[3];
rz(-2.5906467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2962239) q[2];
sx q[2];
rz(-1.8989398) q[2];
sx q[2];
rz(2.5073012) q[2];
rz(-0.93572179) q[3];
sx q[3];
rz(-0.24565419) q[3];
sx q[3];
rz(0.74265695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82035404) q[0];
sx q[0];
rz(-1.5411493) q[0];
sx q[0];
rz(0.4441922) q[0];
rz(2.3815637) q[1];
sx q[1];
rz(-1.1385695) q[1];
sx q[1];
rz(0.98145032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0652577) q[0];
sx q[0];
rz(-1.0629356) q[0];
sx q[0];
rz(1.9064376) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18155988) q[2];
sx q[2];
rz(-1.4019792) q[2];
sx q[2];
rz(3.0491875) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1552298) q[1];
sx q[1];
rz(-1.6587509) q[1];
sx q[1];
rz(-0.4390688) q[1];
rz(-pi) q[2];
rz(0.4325474) q[3];
sx q[3];
rz(-1.8125696) q[3];
sx q[3];
rz(-1.4885308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0280219) q[2];
sx q[2];
rz(-2.5271723) q[2];
sx q[2];
rz(-1.9363972) q[2];
rz(2.6049854) q[3];
sx q[3];
rz(-1.2380995) q[3];
sx q[3];
rz(-1.7369778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2236915) q[0];
sx q[0];
rz(-1.8603928) q[0];
sx q[0];
rz(0.83876383) q[0];
rz(-0.63610786) q[1];
sx q[1];
rz(-1.6517703) q[1];
sx q[1];
rz(0.020523358) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3394525) q[0];
sx q[0];
rz(-1.0217371) q[0];
sx q[0];
rz(0.86717506) q[0];
rz(-pi) q[1];
rz(-0.62600796) q[2];
sx q[2];
rz(-1.1824208) q[2];
sx q[2];
rz(3.1386536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5387303) q[1];
sx q[1];
rz(-2.7285721) q[1];
sx q[1];
rz(-3.1228423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9841524) q[3];
sx q[3];
rz(-1.8529112) q[3];
sx q[3];
rz(-1.071014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8013578) q[2];
sx q[2];
rz(-1.5284208) q[2];
sx q[2];
rz(2.197263) q[2];
rz(1.0229735) q[3];
sx q[3];
rz(-1.2228271) q[3];
sx q[3];
rz(-2.003722) q[3];
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
sx q[0];
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
rz(1.9842072) q[0];
rz(1.6429139) q[1];
sx q[1];
rz(-1.6694371) q[1];
sx q[1];
rz(-1.1882943) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4469915) q[0];
sx q[0];
rz(-2.0216536) q[0];
sx q[0];
rz(0.53038181) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5929668) q[2];
sx q[2];
rz(-1.5323735) q[2];
sx q[2];
rz(0.98603934) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.374158) q[1];
sx q[1];
rz(-2.1302472) q[1];
sx q[1];
rz(2.4591695) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1340302) q[3];
sx q[3];
rz(-2.7705857) q[3];
sx q[3];
rz(1.8414258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1184065) q[2];
sx q[2];
rz(-1.964317) q[2];
sx q[2];
rz(0.6905306) q[2];
rz(1.1150507) q[3];
sx q[3];
rz(-2.3670022) q[3];
sx q[3];
rz(1.5007277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046415064) q[0];
sx q[0];
rz(-1.4354118) q[0];
sx q[0];
rz(1.0272367) q[0];
rz(1.8401624) q[1];
sx q[1];
rz(-2.4899028) q[1];
sx q[1];
rz(-2.8177736) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3394748) q[0];
sx q[0];
rz(-2.3314752) q[0];
sx q[0];
rz(2.1359865) q[0];
x q[1];
rz(2.1236046) q[2];
sx q[2];
rz(-1.5587285) q[2];
sx q[2];
rz(2.7928305) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4253772) q[1];
sx q[1];
rz(-2.1105123) q[1];
sx q[1];
rz(1.9464689) q[1];
rz(-pi) q[2];
rz(1.2965167) q[3];
sx q[3];
rz(-0.52805942) q[3];
sx q[3];
rz(-1.5825281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41574898) q[2];
sx q[2];
rz(-0.21012935) q[2];
sx q[2];
rz(0.48247639) q[2];
rz(-1.9735362) q[3];
sx q[3];
rz(-1.9246293) q[3];
sx q[3];
rz(-0.67679685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2060858) q[0];
sx q[0];
rz(-2.2914903) q[0];
sx q[0];
rz(2.2555943) q[0];
rz(2.50792) q[1];
sx q[1];
rz(-0.88992563) q[1];
sx q[1];
rz(-0.69127965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1544428) q[0];
sx q[0];
rz(-1.6073174) q[0];
sx q[0];
rz(2.4958688) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5070314) q[2];
sx q[2];
rz(-2.261613) q[2];
sx q[2];
rz(1.4379759) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4659652) q[1];
sx q[1];
rz(-1.4838929) q[1];
sx q[1];
rz(-2.3244832) q[1];
x q[2];
rz(-1.8227303) q[3];
sx q[3];
rz(-0.79753424) q[3];
sx q[3];
rz(-3.0424527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0014235) q[2];
sx q[2];
rz(-2.3052577) q[2];
sx q[2];
rz(0.26958618) q[2];
rz(-2.1932898) q[3];
sx q[3];
rz(-1.5323261) q[3];
sx q[3];
rz(-1.3986826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7798994) q[0];
sx q[0];
rz(-2.7044856) q[0];
sx q[0];
rz(-1.0070739) q[0];
rz(0.40564793) q[1];
sx q[1];
rz(-0.59527731) q[1];
sx q[1];
rz(0.55353037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2562228) q[0];
sx q[0];
rz(-1.6281343) q[0];
sx q[0];
rz(-1.8256515) q[0];
rz(-2.361176) q[2];
sx q[2];
rz(-2.4274106) q[2];
sx q[2];
rz(1.2250021) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1055401) q[1];
sx q[1];
rz(-0.80784384) q[1];
sx q[1];
rz(0.64942645) q[1];
rz(-pi) q[2];
rz(2.7615879) q[3];
sx q[3];
rz(-1.5940574) q[3];
sx q[3];
rz(-1.8398566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8290528) q[2];
sx q[2];
rz(-2.0487787) q[2];
sx q[2];
rz(1.9276169) q[2];
rz(2.35516) q[3];
sx q[3];
rz(-1.7460456) q[3];
sx q[3];
rz(-0.023155183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19602747) q[0];
sx q[0];
rz(-2.8493311) q[0];
sx q[0];
rz(-1.8096402) q[0];
rz(2.5281483) q[1];
sx q[1];
rz(-2.1211801) q[1];
sx q[1];
rz(-0.44949284) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17380781) q[0];
sx q[0];
rz(-2.0575843) q[0];
sx q[0];
rz(-0.094102458) q[0];
rz(-pi) q[1];
rz(-1.23486) q[2];
sx q[2];
rz(-1.6523696) q[2];
sx q[2];
rz(0.60816607) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39169381) q[1];
sx q[1];
rz(-2.3428681) q[1];
sx q[1];
rz(2.9500089) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.883955) q[3];
sx q[3];
rz(-1.5143359) q[3];
sx q[3];
rz(-0.27654058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8248262) q[2];
sx q[2];
rz(-0.68244857) q[2];
sx q[2];
rz(-2.5332434) q[2];
rz(1.1634722) q[3];
sx q[3];
rz(-1.3694265) q[3];
sx q[3];
rz(-2.5782862) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0332396) q[0];
sx q[0];
rz(-0.9110564) q[0];
sx q[0];
rz(2.5392927) q[0];
rz(0.96744084) q[1];
sx q[1];
rz(-2.2106705) q[1];
sx q[1];
rz(-2.0379351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33537827) q[0];
sx q[0];
rz(-2.4455482) q[0];
sx q[0];
rz(0.93675128) q[0];
x q[1];
rz(-3.121762) q[2];
sx q[2];
rz(-2.5049372) q[2];
sx q[2];
rz(-3.0002468) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.98075726) q[1];
sx q[1];
rz(-2.6904562) q[1];
sx q[1];
rz(-1.5680997) q[1];
rz(1.6742485) q[3];
sx q[3];
rz(-2.561063) q[3];
sx q[3];
rz(0.79642297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2063107) q[2];
sx q[2];
rz(-1.3056359) q[2];
sx q[2];
rz(0.3375816) q[2];
rz(-2.857699) q[3];
sx q[3];
rz(-2.4604535) q[3];
sx q[3];
rz(-0.18686992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(1.5934481) q[0];
sx q[0];
rz(-1.633506) q[0];
sx q[0];
rz(-3.0737851) q[0];
rz(-1.1495122) q[1];
sx q[1];
rz(-1.5510473) q[1];
sx q[1];
rz(-0.4745208) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10013469) q[0];
sx q[0];
rz(-2.9024419) q[0];
sx q[0];
rz(-1.6275703) q[0];
x q[1];
rz(-2.6340747) q[2];
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
rz(-pi/2) q[0];
rz(-0.671512) q[1];
sx q[1];
rz(-1.5776411) q[1];
sx q[1];
rz(1.4840625) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3274997) q[3];
sx q[3];
rz(-1.0225227) q[3];
sx q[3];
rz(-0.89422885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6600251) q[2];
sx q[2];
rz(-2.2109172) q[2];
sx q[2];
rz(-1.0732667) q[2];
rz(2.977071) q[3];
sx q[3];
rz(-1.8142895) q[3];
sx q[3];
rz(2.8209414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58707033) q[0];
sx q[0];
rz(-1.5759435) q[0];
sx q[0];
rz(-1.6389621) q[0];
rz(1.5578237) q[1];
sx q[1];
rz(-2.0638034) q[1];
sx q[1];
rz(2.6149909) q[1];
rz(-2.4921992) q[2];
sx q[2];
rz(-2.5838701) q[2];
sx q[2];
rz(2.0412847) q[2];
rz(2.9992793) q[3];
sx q[3];
rz(-1.489779) q[3];
sx q[3];
rz(-1.6708556) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
