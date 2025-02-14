OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.6505721) q[0];
sx q[0];
rz(0.84492961) q[0];
sx q[0];
rz(6.9965811) q[0];
rz(2.9333935) q[1];
sx q[1];
rz(2.9543076) q[1];
sx q[1];
rz(2.0050144) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0381361) q[0];
sx q[0];
rz(-2.2503997) q[0];
sx q[0];
rz(-2.8299299) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6460674) q[2];
sx q[2];
rz(-2.4589361) q[2];
sx q[2];
rz(-0.45479506) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40219992) q[1];
sx q[1];
rz(-0.2570411) q[1];
sx q[1];
rz(-0.57769026) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4376115) q[3];
sx q[3];
rz(-1.2696506) q[3];
sx q[3];
rz(-2.4011321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1897757) q[2];
sx q[2];
rz(-1.2729898) q[2];
sx q[2];
rz(-0.58958685) q[2];
rz(0.014852614) q[3];
sx q[3];
rz(-1.868052) q[3];
sx q[3];
rz(-2.5428037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1349161) q[0];
sx q[0];
rz(-1.0449469) q[0];
sx q[0];
rz(2.6892804) q[0];
rz(2.2882838) q[1];
sx q[1];
rz(-0.49214688) q[1];
sx q[1];
rz(0.36202994) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2483082) q[0];
sx q[0];
rz(-2.4442857) q[0];
sx q[0];
rz(2.6261283) q[0];
rz(-pi) q[1];
rz(-0.49136038) q[2];
sx q[2];
rz(-1.4213014) q[2];
sx q[2];
rz(0.66505611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4072864) q[1];
sx q[1];
rz(-0.9883259) q[1];
sx q[1];
rz(2.4206619) q[1];
rz(-0.70254247) q[3];
sx q[3];
rz(-1.4558534) q[3];
sx q[3];
rz(3.0963932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33874685) q[2];
sx q[2];
rz(-0.57463988) q[2];
sx q[2];
rz(2.23526) q[2];
rz(-1.686442) q[3];
sx q[3];
rz(-0.95821277) q[3];
sx q[3];
rz(-1.5866535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872658) q[0];
sx q[0];
rz(-1.9540906) q[0];
sx q[0];
rz(-0.97386709) q[0];
rz(2.5663238) q[1];
sx q[1];
rz(-1.560805) q[1];
sx q[1];
rz(-1.5781933) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78689903) q[0];
sx q[0];
rz(-1.5074566) q[0];
sx q[0];
rz(0.022149274) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.695486) q[2];
sx q[2];
rz(-0.89274721) q[2];
sx q[2];
rz(2.7715832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0583349) q[1];
sx q[1];
rz(-2.2607723) q[1];
sx q[1];
rz(-0.95466701) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3493912) q[3];
sx q[3];
rz(-2.3292037) q[3];
sx q[3];
rz(-1.0426857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32597184) q[2];
sx q[2];
rz(-1.9591363) q[2];
sx q[2];
rz(1.5284485) q[2];
rz(-3.1260887) q[3];
sx q[3];
rz(-1.1752693) q[3];
sx q[3];
rz(-1.6541245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4267321) q[0];
sx q[0];
rz(-0.69789129) q[0];
sx q[0];
rz(-0.062407169) q[0];
rz(1.1563673) q[1];
sx q[1];
rz(-2.195475) q[1];
sx q[1];
rz(0.83414331) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5204374) q[0];
sx q[0];
rz(-0.59323673) q[0];
sx q[0];
rz(1.0246681) q[0];
x q[1];
rz(-1.1961924) q[2];
sx q[2];
rz(-0.13557493) q[2];
sx q[2];
rz(-1.6983528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3662009) q[1];
sx q[1];
rz(-1.6244666) q[1];
sx q[1];
rz(2.1691775) q[1];
x q[2];
rz(2.9328654) q[3];
sx q[3];
rz(-2.0794915) q[3];
sx q[3];
rz(1.3172447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1674898) q[2];
sx q[2];
rz(-2.9656599) q[2];
sx q[2];
rz(1.2220194) q[2];
rz(-2.7410298) q[3];
sx q[3];
rz(-1.0223072) q[3];
sx q[3];
rz(0.21489531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.062531384) q[0];
sx q[0];
rz(-0.65199861) q[0];
sx q[0];
rz(-2.8849211) q[0];
rz(-1.3447064) q[1];
sx q[1];
rz(-1.2213629) q[1];
sx q[1];
rz(1.7324956) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.870011) q[0];
sx q[0];
rz(-0.74423941) q[0];
sx q[0];
rz(-1.1209473) q[0];
rz(-1.8518894) q[2];
sx q[2];
rz(-1.7651193) q[2];
sx q[2];
rz(1.3240726) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0540648) q[1];
sx q[1];
rz(-1.1502741) q[1];
sx q[1];
rz(1.6144714) q[1];
rz(-pi) q[2];
rz(0.6677169) q[3];
sx q[3];
rz(-1.8103292) q[3];
sx q[3];
rz(-0.67090494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1008489) q[2];
sx q[2];
rz(-0.70793968) q[2];
sx q[2];
rz(2.9948998) q[2];
rz(1.7723068) q[3];
sx q[3];
rz(-0.36543235) q[3];
sx q[3];
rz(-2.9036314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9323102) q[0];
sx q[0];
rz(-2.3754061) q[0];
sx q[0];
rz(-2.3796418) q[0];
rz(-0.85488287) q[1];
sx q[1];
rz(-1.2763174) q[1];
sx q[1];
rz(-0.68861047) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3935034) q[0];
sx q[0];
rz(-2.6660109) q[0];
sx q[0];
rz(1.7548496) q[0];
rz(-pi) q[1];
rz(1.5916078) q[2];
sx q[2];
rz(-1.4574798) q[2];
sx q[2];
rz(-1.7422388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.411158) q[1];
sx q[1];
rz(-1.467492) q[1];
sx q[1];
rz(1.6367326) q[1];
x q[2];
rz(0.9290885) q[3];
sx q[3];
rz(-2.6179465) q[3];
sx q[3];
rz(-1.3646082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7410572) q[2];
sx q[2];
rz(-2.1964938) q[2];
sx q[2];
rz(-2.1607024) q[2];
rz(-0.26997057) q[3];
sx q[3];
rz(-1.910784) q[3];
sx q[3];
rz(-0.009416906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16801676) q[0];
sx q[0];
rz(-1.4485757) q[0];
sx q[0];
rz(0.5434522) q[0];
rz(-1.5501529) q[1];
sx q[1];
rz(-0.88528577) q[1];
sx q[1];
rz(-0.98522225) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0854093) q[0];
sx q[0];
rz(-2.0517573) q[0];
sx q[0];
rz(-1.3113326) q[0];
x q[1];
rz(-1.0205054) q[2];
sx q[2];
rz(-0.63945635) q[2];
sx q[2];
rz(0.24580641) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7581525) q[1];
sx q[1];
rz(-1.1894656) q[1];
sx q[1];
rz(-0.9419946) q[1];
rz(2.6292218) q[3];
sx q[3];
rz(-1.2768043) q[3];
sx q[3];
rz(1.9749354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.33124179) q[2];
sx q[2];
rz(-0.50464973) q[2];
sx q[2];
rz(-1.3387559) q[2];
rz(1.6797558) q[3];
sx q[3];
rz(-1.5906518) q[3];
sx q[3];
rz(-2.4707879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93002334) q[0];
sx q[0];
rz(-1.1299364) q[0];
sx q[0];
rz(-1.2343963) q[0];
rz(-1.7878388) q[1];
sx q[1];
rz(-2.6852971) q[1];
sx q[1];
rz(3.0455468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1559469) q[0];
sx q[0];
rz(-1.1178769) q[0];
sx q[0];
rz(0.80891305) q[0];
x q[1];
rz(0.12577082) q[2];
sx q[2];
rz(-0.54361225) q[2];
sx q[2];
rz(2.8961492) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2224642) q[1];
sx q[1];
rz(-1.9507196) q[1];
sx q[1];
rz(1.6187731) q[1];
rz(1.6173564) q[3];
sx q[3];
rz(-0.73107204) q[3];
sx q[3];
rz(0.94207596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3357521) q[2];
sx q[2];
rz(-2.7554607) q[2];
sx q[2];
rz(-2.0167572) q[2];
rz(2.2624894) q[3];
sx q[3];
rz(-1.9156888) q[3];
sx q[3];
rz(-0.44900289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6833471) q[0];
sx q[0];
rz(-1.6793716) q[0];
sx q[0];
rz(0.71518389) q[0];
rz(-0.31582754) q[1];
sx q[1];
rz(-2.4555457) q[1];
sx q[1];
rz(2.979523) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2020362) q[0];
sx q[0];
rz(-1.9959772) q[0];
sx q[0];
rz(-2.3935938) q[0];
x q[1];
rz(-2.8388073) q[2];
sx q[2];
rz(-2.5368779) q[2];
sx q[2];
rz(-0.59689891) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4906076) q[1];
sx q[1];
rz(-0.7443634) q[1];
sx q[1];
rz(0.26000826) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8211179) q[3];
sx q[3];
rz(-1.8928013) q[3];
sx q[3];
rz(0.97089773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7519303) q[2];
sx q[2];
rz(-1.1952362) q[2];
sx q[2];
rz(-1.8771089) q[2];
rz(1.8090931) q[3];
sx q[3];
rz(-1.2800565) q[3];
sx q[3];
rz(2.8247824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(1.7952809) q[0];
sx q[0];
rz(-2.3667211) q[0];
sx q[0];
rz(-2.0922022) q[0];
rz(0.28114444) q[1];
sx q[1];
rz(-1.0422336) q[1];
sx q[1];
rz(1.3220471) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41073179) q[0];
sx q[0];
rz(-1.3202892) q[0];
sx q[0];
rz(2.5772472) q[0];
rz(1.0820002) q[2];
sx q[2];
rz(-0.99129516) q[2];
sx q[2];
rz(2.3129675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7483298) q[1];
sx q[1];
rz(-2.7460881) q[1];
sx q[1];
rz(1.1480182) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7554402) q[3];
sx q[3];
rz(-1.6623467) q[3];
sx q[3];
rz(-1.1923252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37810024) q[2];
sx q[2];
rz(-1.3214448) q[2];
sx q[2];
rz(3.0901001) q[2];
rz(-1.3242877) q[3];
sx q[3];
rz(-1.4824425) q[3];
sx q[3];
rz(1.2073368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1235724) q[0];
sx q[0];
rz(-1.3326895) q[0];
sx q[0];
rz(1.4154758) q[0];
rz(-2.3208658) q[1];
sx q[1];
rz(-1.4438933) q[1];
sx q[1];
rz(-1.874598) q[1];
rz(-0.61931284) q[2];
sx q[2];
rz(-1.7650585) q[2];
sx q[2];
rz(-1.4256918) q[2];
rz(0.5140082) q[3];
sx q[3];
rz(-1.2759943) q[3];
sx q[3];
rz(-1.3392824) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
