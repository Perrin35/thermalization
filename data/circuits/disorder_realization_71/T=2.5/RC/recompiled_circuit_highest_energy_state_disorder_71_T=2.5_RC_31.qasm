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
rz(-2.5833997) q[0];
sx q[0];
rz(-0.71891958) q[0];
sx q[0];
rz(2.4655226) q[0];
rz(-2.8683635) q[1];
sx q[1];
rz(-0.9571119) q[1];
sx q[1];
rz(-0.91125429) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34150037) q[0];
sx q[0];
rz(-1.5908634) q[0];
sx q[0];
rz(1.5930727) q[0];
rz(-pi) q[1];
x q[1];
rz(0.010013117) q[2];
sx q[2];
rz(-0.62558936) q[2];
sx q[2];
rz(2.6715476) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76409066) q[1];
sx q[1];
rz(-1.5385748) q[1];
sx q[1];
rz(-0.5733787) q[1];
rz(-1.4449213) q[3];
sx q[3];
rz(-2.3231703) q[3];
sx q[3];
rz(-2.6361806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91997826) q[2];
sx q[2];
rz(-0.45404926) q[2];
sx q[2];
rz(-0.81217074) q[2];
rz(-0.075856097) q[3];
sx q[3];
rz(-2.4935738) q[3];
sx q[3];
rz(-0.59504741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49038637) q[0];
sx q[0];
rz(-2.6631329) q[0];
sx q[0];
rz(-1.2745717) q[0];
rz(-1.6365341) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(-1.9777745) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75117373) q[0];
sx q[0];
rz(-1.112022) q[0];
sx q[0];
rz(2.075688) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69804783) q[2];
sx q[2];
rz(-0.62605941) q[2];
sx q[2];
rz(-1.2241838) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3006176) q[1];
sx q[1];
rz(-1.3625486) q[1];
sx q[1];
rz(0.62842259) q[1];
x q[2];
rz(1.9979565) q[3];
sx q[3];
rz(-1.2795951) q[3];
sx q[3];
rz(1.2808275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85094467) q[2];
sx q[2];
rz(-3.1182351) q[2];
sx q[2];
rz(1.5811496) q[2];
rz(-0.23049878) q[3];
sx q[3];
rz(-1.0483619) q[3];
sx q[3];
rz(0.44052625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5480963) q[0];
sx q[0];
rz(-0.31262147) q[0];
sx q[0];
rz(-0.12259677) q[0];
rz(0.7971881) q[1];
sx q[1];
rz(-1.0120579) q[1];
sx q[1];
rz(-0.0880934) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0038550475) q[0];
sx q[0];
rz(-0.81043506) q[0];
sx q[0];
rz(1.2146945) q[0];
rz(-pi) q[1];
rz(1.9162019) q[2];
sx q[2];
rz(-2.900335) q[2];
sx q[2];
rz(0.119482) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8338203) q[1];
sx q[1];
rz(-1.7227748) q[1];
sx q[1];
rz(0.13843468) q[1];
x q[2];
rz(-2.5949305) q[3];
sx q[3];
rz(-2.529749) q[3];
sx q[3];
rz(2.8091968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88283551) q[2];
sx q[2];
rz(-1.5587403) q[2];
sx q[2];
rz(1.8176414) q[2];
rz(-0.49805182) q[3];
sx q[3];
rz(-2.1784454) q[3];
sx q[3];
rz(2.3948885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.8525456) q[0];
sx q[0];
rz(-2.9883224) q[0];
sx q[0];
rz(-2.9943941) q[0];
rz(2.4920801) q[1];
sx q[1];
rz(-1.5107061) q[1];
sx q[1];
rz(-1.6726327) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29028673) q[0];
sx q[0];
rz(-1.3477579) q[0];
sx q[0];
rz(-0.0073435535) q[0];
rz(1.9542404) q[2];
sx q[2];
rz(-2.6818775) q[2];
sx q[2];
rz(1.5292668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.34079724) q[1];
sx q[1];
rz(-2.1873475) q[1];
sx q[1];
rz(0.430213) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5178403) q[3];
sx q[3];
rz(-0.68188462) q[3];
sx q[3];
rz(2.9870913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9327675) q[2];
sx q[2];
rz(-0.494445) q[2];
sx q[2];
rz(-0.24173582) q[2];
rz(0.73511165) q[3];
sx q[3];
rz(-1.7416411) q[3];
sx q[3];
rz(-2.476695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.46297896) q[0];
sx q[0];
rz(-1.047387) q[0];
sx q[0];
rz(-3.0188634) q[0];
rz(0.8853451) q[1];
sx q[1];
rz(-2.131999) q[1];
sx q[1];
rz(-2.5811894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2764328) q[0];
sx q[0];
rz(-1.7483646) q[0];
sx q[0];
rz(2.0899169) q[0];
x q[1];
rz(-2.5053315) q[2];
sx q[2];
rz(-0.36379746) q[2];
sx q[2];
rz(0.76278245) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0593157) q[1];
sx q[1];
rz(-2.6059285) q[1];
sx q[1];
rz(1.0470864) q[1];
rz(-pi) q[2];
rz(-0.82920977) q[3];
sx q[3];
rz(-2.5253339) q[3];
sx q[3];
rz(-2.5065638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5237328) q[2];
sx q[2];
rz(-0.79404074) q[2];
sx q[2];
rz(-2.9368371) q[2];
rz(2.2023885) q[3];
sx q[3];
rz(-1.3699646) q[3];
sx q[3];
rz(-3.052616) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8646249) q[0];
sx q[0];
rz(-0.11700103) q[0];
sx q[0];
rz(0.65698874) q[0];
rz(-2.9981546) q[1];
sx q[1];
rz(-1.5851494) q[1];
sx q[1];
rz(0.73854804) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1114233) q[0];
sx q[0];
rz(-2.2403754) q[0];
sx q[0];
rz(-1.4233146) q[0];
x q[1];
rz(-1.1656649) q[2];
sx q[2];
rz(-0.52743334) q[2];
sx q[2];
rz(-0.95587522) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2966531) q[1];
sx q[1];
rz(-2.2187738) q[1];
sx q[1];
rz(2.7428127) q[1];
x q[2];
rz(2.5127327) q[3];
sx q[3];
rz(-0.81506461) q[3];
sx q[3];
rz(-3.0185543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2584381) q[2];
sx q[2];
rz(-0.65885764) q[2];
sx q[2];
rz(3.1385885) q[2];
rz(-0.7891807) q[3];
sx q[3];
rz(-2.7572258) q[3];
sx q[3];
rz(-1.9743617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060519144) q[0];
sx q[0];
rz(-0.34927148) q[0];
sx q[0];
rz(0.14807598) q[0];
rz(-2.608346) q[1];
sx q[1];
rz(-2.8956469) q[1];
sx q[1];
rz(-1.6291133) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142198) q[0];
sx q[0];
rz(-2.5974869) q[0];
sx q[0];
rz(-2.6065488) q[0];
x q[1];
rz(-2.2165636) q[2];
sx q[2];
rz(-2.2730977) q[2];
sx q[2];
rz(1.4214732) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5809734) q[1];
sx q[1];
rz(-1.4500014) q[1];
sx q[1];
rz(-1.2289775) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2929271) q[3];
sx q[3];
rz(-1.7300528) q[3];
sx q[3];
rz(-1.4383663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6988354) q[2];
sx q[2];
rz(-1.7008702) q[2];
sx q[2];
rz(0.093078144) q[2];
rz(0.062151521) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(-1.2232346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1208948) q[0];
sx q[0];
rz(-0.59497213) q[0];
sx q[0];
rz(0.66463941) q[0];
rz(1.7918034) q[1];
sx q[1];
rz(-2.2638958) q[1];
sx q[1];
rz(0.68914366) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6751373) q[0];
sx q[0];
rz(-0.27091089) q[0];
sx q[0];
rz(-2.6938426) q[0];
rz(-pi) q[1];
rz(0.99456878) q[2];
sx q[2];
rz(-2.3625018) q[2];
sx q[2];
rz(1.3887203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41708159) q[1];
sx q[1];
rz(-1.8838091) q[1];
sx q[1];
rz(1.723812) q[1];
rz(1.6313305) q[3];
sx q[3];
rz(-1.406296) q[3];
sx q[3];
rz(-1.9609083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.050921116) q[2];
sx q[2];
rz(-0.65616578) q[2];
sx q[2];
rz(-1.3760759) q[2];
rz(1.225166) q[3];
sx q[3];
rz(-0.71198946) q[3];
sx q[3];
rz(-0.041697748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-3.0390778) q[0];
sx q[0];
rz(-3.097105) q[0];
sx q[0];
rz(-2.5503889) q[0];
rz(-0.2419596) q[1];
sx q[1];
rz(-1.0958902) q[1];
sx q[1];
rz(-2.718149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0288182) q[0];
sx q[0];
rz(-1.2679546) q[0];
sx q[0];
rz(1.8870728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4884004) q[2];
sx q[2];
rz(-1.6624221) q[2];
sx q[2];
rz(2.5727401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6953227) q[1];
sx q[1];
rz(-0.62275089) q[1];
sx q[1];
rz(-1.0458617) q[1];
rz(-pi) q[2];
rz(-1.1657646) q[3];
sx q[3];
rz(-0.99294239) q[3];
sx q[3];
rz(-2.809643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5847136) q[2];
sx q[2];
rz(-2.5355279) q[2];
sx q[2];
rz(-1.9752183) q[2];
rz(-2.0139458) q[3];
sx q[3];
rz(-2.1731264) q[3];
sx q[3];
rz(0.52496547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.128085) q[0];
sx q[0];
rz(-0.43376827) q[0];
sx q[0];
rz(-0.67310131) q[0];
rz(0.85064864) q[1];
sx q[1];
rz(-1.3530082) q[1];
sx q[1];
rz(0.32352111) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0304383) q[0];
sx q[0];
rz(-1.9435099) q[0];
sx q[0];
rz(-0.26813676) q[0];
rz(-2.0490626) q[2];
sx q[2];
rz(-1.8169122) q[2];
sx q[2];
rz(0.1637295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9266745) q[1];
sx q[1];
rz(-1.9028544) q[1];
sx q[1];
rz(2.3658793) q[1];
rz(-pi) q[2];
rz(-2.5577522) q[3];
sx q[3];
rz(-1.5642318) q[3];
sx q[3];
rz(-3.1232426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2181776) q[2];
sx q[2];
rz(-0.18109334) q[2];
sx q[2];
rz(-3.0695445) q[2];
rz(-2.1273023) q[3];
sx q[3];
rz(-0.91599661) q[3];
sx q[3];
rz(2.5924276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88076787) q[0];
sx q[0];
rz(-1.7171971) q[0];
sx q[0];
rz(1.1203753) q[0];
rz(-0.33521677) q[1];
sx q[1];
rz(-2.4991279) q[1];
sx q[1];
rz(2.0229708) q[1];
rz(-1.3041244) q[2];
sx q[2];
rz(-0.49378569) q[2];
sx q[2];
rz(2.7616382) q[2];
rz(0.14306457) q[3];
sx q[3];
rz(-1.5456556) q[3];
sx q[3];
rz(1.7329334) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
