OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.50897941) q[0];
sx q[0];
rz(-0.84492961) q[0];
sx q[0];
rz(0.71339575) q[0];
rz(2.9333935) q[1];
sx q[1];
rz(2.9543076) q[1];
sx q[1];
rz(2.0050144) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644585) q[0];
sx q[0];
rz(-2.4043879) q[0];
sx q[0];
rz(1.2080824) q[0];
rz(-0.061066433) q[2];
sx q[2];
rz(-0.89043987) q[2];
sx q[2];
rz(-0.55168569) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7311095) q[1];
sx q[1];
rz(-1.7100733) q[1];
sx q[1];
rz(2.9248505) q[1];
x q[2];
rz(0.40402255) q[3];
sx q[3];
rz(-0.32846132) q[3];
sx q[3];
rz(2.8253909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95181695) q[2];
sx q[2];
rz(-1.8686029) q[2];
sx q[2];
rz(-0.58958685) q[2];
rz(3.12674) q[3];
sx q[3];
rz(-1.2735406) q[3];
sx q[3];
rz(0.59878892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1349161) q[0];
sx q[0];
rz(-1.0449469) q[0];
sx q[0];
rz(2.6892804) q[0];
rz(0.85330883) q[1];
sx q[1];
rz(-0.49214688) q[1];
sx q[1];
rz(-0.36202994) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2483082) q[0];
sx q[0];
rz(-0.69730699) q[0];
sx q[0];
rz(-0.5154644) q[0];
x q[1];
rz(-0.30899991) q[2];
sx q[2];
rz(-2.6297792) q[2];
sx q[2];
rz(1.1771894) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7343062) q[1];
sx q[1];
rz(-0.9883259) q[1];
sx q[1];
rz(-2.4206619) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4206616) q[3];
sx q[3];
rz(-2.2677652) q[3];
sx q[3];
rz(1.622395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.33874685) q[2];
sx q[2];
rz(-2.5669528) q[2];
sx q[2];
rz(0.90633264) q[2];
rz(1.4551506) q[3];
sx q[3];
rz(-0.95821277) q[3];
sx q[3];
rz(-1.5866535) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543268) q[0];
sx q[0];
rz(-1.187502) q[0];
sx q[0];
rz(-0.97386709) q[0];
rz(0.57526881) q[1];
sx q[1];
rz(-1.560805) q[1];
sx q[1];
rz(-1.5633993) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78689903) q[0];
sx q[0];
rz(-1.5074566) q[0];
sx q[0];
rz(-0.022149274) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0625774) q[2];
sx q[2];
rz(-0.79177982) q[2];
sx q[2];
rz(2.1201652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0583349) q[1];
sx q[1];
rz(-2.2607723) q[1];
sx q[1];
rz(0.95466701) q[1];
rz(-pi) q[2];
rz(0.22776484) q[3];
sx q[3];
rz(-0.78383369) q[3];
sx q[3];
rz(1.3589825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8156208) q[2];
sx q[2];
rz(-1.1824563) q[2];
sx q[2];
rz(1.6131442) q[2];
rz(-0.01550393) q[3];
sx q[3];
rz(-1.9663234) q[3];
sx q[3];
rz(1.4874682) q[3];
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
rz(1.7148606) q[0];
sx q[0];
rz(-0.69789129) q[0];
sx q[0];
rz(3.0791855) q[0];
rz(-1.9852253) q[1];
sx q[1];
rz(-2.195475) q[1];
sx q[1];
rz(-2.3074493) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011403712) q[0];
sx q[0];
rz(-2.0688648) q[0];
sx q[0];
rz(2.8047436) q[0];
x q[1];
rz(1.4445199) q[2];
sx q[2];
rz(-1.6202721) q[2];
sx q[2];
rz(-0.49903185) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4247113) q[1];
sx q[1];
rz(-0.60048999) q[1];
sx q[1];
rz(-1.6658777) q[1];
x q[2];
rz(2.9328654) q[3];
sx q[3];
rz(-2.0794915) q[3];
sx q[3];
rz(1.3172447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1674898) q[2];
sx q[2];
rz(-2.9656599) q[2];
sx q[2];
rz(1.2220194) q[2];
rz(0.40056285) q[3];
sx q[3];
rz(-2.1192854) q[3];
sx q[3];
rz(2.9266973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.9202298) q[1];
sx q[1];
rz(-1.7324956) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2715816) q[0];
sx q[0];
rz(-0.74423941) q[0];
sx q[0];
rz(-2.0206454) q[0];
rz(-1.8518894) q[2];
sx q[2];
rz(-1.7651193) q[2];
sx q[2];
rz(1.3240726) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46542922) q[1];
sx q[1];
rz(-1.6106642) q[1];
sx q[1];
rz(2.7207147) q[1];
rz(-pi) q[2];
rz(2.4738757) q[3];
sx q[3];
rz(-1.3312634) q[3];
sx q[3];
rz(-0.67090494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.04074374) q[2];
sx q[2];
rz(-2.433653) q[2];
sx q[2];
rz(-2.9948998) q[2];
rz(-1.3692859) q[3];
sx q[3];
rz(-0.36543235) q[3];
sx q[3];
rz(-2.9036314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20928243) q[0];
sx q[0];
rz(-2.3754061) q[0];
sx q[0];
rz(-0.76195088) q[0];
rz(0.85488287) q[1];
sx q[1];
rz(-1.8652752) q[1];
sx q[1];
rz(2.4529822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013278339) q[0];
sx q[0];
rz(-1.6546895) q[0];
sx q[0];
rz(1.1021139) q[0];
rz(0.11334085) q[2];
sx q[2];
rz(-1.5501184) q[2];
sx q[2];
rz(0.17379601) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7304346) q[1];
sx q[1];
rz(-1.6741006) q[1];
sx q[1];
rz(1.6367326) q[1];
x q[2];
rz(-1.1375539) q[3];
sx q[3];
rz(-1.2668303) q[3];
sx q[3];
rz(-0.36813018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40053549) q[2];
sx q[2];
rz(-0.94509882) q[2];
sx q[2];
rz(0.98089027) q[2];
rz(2.8716221) q[3];
sx q[3];
rz(-1.910784) q[3];
sx q[3];
rz(3.1321757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9735759) q[0];
sx q[0];
rz(-1.4485757) q[0];
sx q[0];
rz(2.5981405) q[0];
rz(-1.5914397) q[1];
sx q[1];
rz(-2.2563069) q[1];
sx q[1];
rz(-0.98522225) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0561834) q[0];
sx q[0];
rz(-2.0517573) q[0];
sx q[0];
rz(1.3113326) q[0];
rz(-pi) q[1];
rz(-2.7706861) q[2];
sx q[2];
rz(-2.1044253) q[2];
sx q[2];
rz(2.7346389) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66032672) q[1];
sx q[1];
rz(-0.72176343) q[1];
sx q[1];
rz(0.97246171) q[1];
rz(-1.9051294) q[3];
sx q[3];
rz(-1.0823993) q[3];
sx q[3];
rz(2.5758829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33124179) q[2];
sx q[2];
rz(-0.50464973) q[2];
sx q[2];
rz(-1.8028367) q[2];
rz(1.4618368) q[3];
sx q[3];
rz(-1.5906518) q[3];
sx q[3];
rz(2.4707879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.93002334) q[0];
sx q[0];
rz(-1.1299364) q[0];
sx q[0];
rz(-1.2343963) q[0];
rz(1.3537539) q[1];
sx q[1];
rz(-2.6852971) q[1];
sx q[1];
rz(-0.096045883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1609566) q[0];
sx q[0];
rz(-0.90117878) q[0];
sx q[0];
rz(-2.5494844) q[0];
x q[1];
rz(-2.6014843) q[2];
sx q[2];
rz(-1.5058695) q[2];
sx q[2];
rz(-1.924032) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.36613606) q[1];
sx q[1];
rz(-1.526243) q[1];
sx q[1];
rz(-2.7612727) q[1];
rz(3.0998746) q[3];
sx q[3];
rz(-0.8406958) q[3];
sx q[3];
rz(0.87957013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8058406) q[2];
sx q[2];
rz(-0.38613191) q[2];
sx q[2];
rz(-1.1248355) q[2];
rz(-0.8791033) q[3];
sx q[3];
rz(-1.9156888) q[3];
sx q[3];
rz(2.6925898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.45824555) q[0];
sx q[0];
rz(-1.4622211) q[0];
sx q[0];
rz(0.71518389) q[0];
rz(0.31582754) q[1];
sx q[1];
rz(-2.4555457) q[1];
sx q[1];
rz(-2.979523) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26570854) q[0];
sx q[0];
rz(-2.2390463) q[0];
sx q[0];
rz(-1.0174684) q[0];
rz(-pi) q[1];
rz(0.30278532) q[2];
sx q[2];
rz(-2.5368779) q[2];
sx q[2];
rz(2.5446937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.650985) q[1];
sx q[1];
rz(-0.7443634) q[1];
sx q[1];
rz(-2.8815844) q[1];
x q[2];
rz(-0.32047477) q[3];
sx q[3];
rz(-1.2487914) q[3];
sx q[3];
rz(-2.1706949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7519303) q[2];
sx q[2];
rz(-1.9463564) q[2];
sx q[2];
rz(1.8771089) q[2];
rz(1.8090931) q[3];
sx q[3];
rz(-1.8615362) q[3];
sx q[3];
rz(-2.8247824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463117) q[0];
sx q[0];
rz(-0.77487159) q[0];
sx q[0];
rz(-2.0922022) q[0];
rz(-0.28114444) q[1];
sx q[1];
rz(-1.0422336) q[1];
sx q[1];
rz(-1.3220471) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8258749) q[0];
sx q[0];
rz(-2.115504) q[0];
sx q[0];
rz(-1.2767351) q[0];
rz(-pi) q[1];
rz(-0.62236592) q[2];
sx q[2];
rz(-2.4020148) q[2];
sx q[2];
rz(-0.057829521) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57105455) q[1];
sx q[1];
rz(-1.4120543) q[1];
sx q[1];
rz(-1.9345933) q[1];
x q[2];
rz(-0.23909823) q[3];
sx q[3];
rz(-0.39632583) q[3];
sx q[3];
rz(-2.9843085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7634924) q[2];
sx q[2];
rz(-1.3214448) q[2];
sx q[2];
rz(3.0901001) q[2];
rz(1.817305) q[3];
sx q[3];
rz(-1.6591502) q[3];
sx q[3];
rz(1.9342559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0180203) q[0];
sx q[0];
rz(-1.8089031) q[0];
sx q[0];
rz(-1.7261169) q[0];
rz(-0.82072683) q[1];
sx q[1];
rz(-1.6976994) q[1];
sx q[1];
rz(1.2669947) q[1];
rz(2.8148094) q[2];
sx q[2];
rz(-2.4963623) q[2];
sx q[2];
rz(0.4095578) q[2];
rz(-1.2352712) q[3];
sx q[3];
rz(-2.060609) q[3];
sx q[3];
rz(0.068908066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
