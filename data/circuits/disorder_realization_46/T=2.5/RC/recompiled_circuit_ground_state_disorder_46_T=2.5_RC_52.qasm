OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0134861) q[0];
sx q[0];
rz(-0.82733893) q[0];
sx q[0];
rz(-1.3878571) q[0];
rz(0.85343051) q[1];
sx q[1];
rz(5.8813385) q[1];
sx q[1];
rz(11.452236) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6809876) q[0];
sx q[0];
rz(-1.3189684) q[0];
sx q[0];
rz(0.84427174) q[0];
x q[1];
rz(1.3783437) q[2];
sx q[2];
rz(-2.2666605) q[2];
sx q[2];
rz(0.025283289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.6889894) q[1];
sx q[1];
rz(-1.5905979) q[1];
sx q[1];
rz(-0.0048586998) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57350343) q[3];
sx q[3];
rz(-0.57799229) q[3];
sx q[3];
rz(-1.0740394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9806597) q[2];
sx q[2];
rz(-2.7200343) q[2];
sx q[2];
rz(1.6932311) q[2];
rz(-0.24250044) q[3];
sx q[3];
rz(-0.91553965) q[3];
sx q[3];
rz(1.8814794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.563797) q[0];
sx q[0];
rz(-1.8255434) q[0];
sx q[0];
rz(-1.0003566) q[0];
rz(-2.4099804) q[1];
sx q[1];
rz(-1.4294521) q[1];
sx q[1];
rz(-1.9614722) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7655752) q[0];
sx q[0];
rz(-1.4808324) q[0];
sx q[0];
rz(-3.0375809) q[0];
x q[1];
rz(1.6890516) q[2];
sx q[2];
rz(-0.27980556) q[2];
sx q[2];
rz(-2.5733833) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6669162) q[1];
sx q[1];
rz(-1.70494) q[1];
sx q[1];
rz(0.0018906126) q[1];
x q[2];
rz(1.2878706) q[3];
sx q[3];
rz(-1.8370421) q[3];
sx q[3];
rz(0.11582813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0483027) q[2];
sx q[2];
rz(-2.2403109) q[2];
sx q[2];
rz(2.062659) q[2];
rz(3.0537187) q[3];
sx q[3];
rz(-1.7885957) q[3];
sx q[3];
rz(-0.24438508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3946149) q[0];
sx q[0];
rz(-1.9786388) q[0];
sx q[0];
rz(-1.8633457) q[0];
rz(-2.6784015) q[1];
sx q[1];
rz(-1.3563503) q[1];
sx q[1];
rz(-0.82122222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0671135) q[0];
sx q[0];
rz(-2.4601106) q[0];
sx q[0];
rz(1.8675748) q[0];
rz(-pi) q[1];
rz(-0.37126705) q[2];
sx q[2];
rz(-1.5131931) q[2];
sx q[2];
rz(1.0641554) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6602409) q[1];
sx q[1];
rz(-1.7180182) q[1];
sx q[1];
rz(-2.6805356) q[1];
rz(-pi) q[2];
rz(-0.65659237) q[3];
sx q[3];
rz(-2.2664605) q[3];
sx q[3];
rz(0.10564239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2440196) q[2];
sx q[2];
rz(-1.579318) q[2];
sx q[2];
rz(2.139034) q[2];
rz(1.2282486) q[3];
sx q[3];
rz(-1.7398261) q[3];
sx q[3];
rz(1.4209411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9035852) q[0];
sx q[0];
rz(-2.0915732) q[0];
sx q[0];
rz(2.8557657) q[0];
rz(-2.8803275) q[1];
sx q[1];
rz(-1.3328054) q[1];
sx q[1];
rz(-0.030698311) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0470025) q[0];
sx q[0];
rz(-1.6436716) q[0];
sx q[0];
rz(-2.8776309) q[0];
x q[1];
rz(0.07918723) q[2];
sx q[2];
rz(-2.2461366) q[2];
sx q[2];
rz(1.5164708) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5324227) q[1];
sx q[1];
rz(-1.318319) q[1];
sx q[1];
rz(2.0308475) q[1];
x q[2];
rz(-1.8844344) q[3];
sx q[3];
rz(-1.6480069) q[3];
sx q[3];
rz(1.9208637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8135445) q[2];
sx q[2];
rz(-1.8009461) q[2];
sx q[2];
rz(3.1385341) q[2];
rz(0.84117738) q[3];
sx q[3];
rz(-0.58924651) q[3];
sx q[3];
rz(0.35471788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4641814) q[0];
sx q[0];
rz(-0.84742904) q[0];
sx q[0];
rz(1.6743976) q[0];
rz(1.6293619) q[1];
sx q[1];
rz(-2.1758175) q[1];
sx q[1];
rz(-0.95219749) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6124355) q[0];
sx q[0];
rz(-1.5780492) q[0];
sx q[0];
rz(0.070242191) q[0];
rz(-pi) q[1];
rz(1.1811851) q[2];
sx q[2];
rz(-2.2305373) q[2];
sx q[2];
rz(0.3241186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4213688) q[1];
sx q[1];
rz(-1.6744057) q[1];
sx q[1];
rz(-3.1006579) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2837379) q[3];
sx q[3];
rz(-0.80227533) q[3];
sx q[3];
rz(1.0209393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3516922) q[2];
sx q[2];
rz(-1.60195) q[2];
sx q[2];
rz(-2.6969625) q[2];
rz(2.6897258) q[3];
sx q[3];
rz(-0.63932747) q[3];
sx q[3];
rz(3.0795081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.708798) q[0];
sx q[0];
rz(-0.72467703) q[0];
sx q[0];
rz(-0.92025796) q[0];
rz(0.96011773) q[1];
sx q[1];
rz(-1.7476387) q[1];
sx q[1];
rz(0.49930176) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1007784) q[0];
sx q[0];
rz(-1.7216604) q[0];
sx q[0];
rz(-0.34025451) q[0];
x q[1];
rz(-2.9174706) q[2];
sx q[2];
rz(-0.67356985) q[2];
sx q[2];
rz(2.3526255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.074452049) q[1];
sx q[1];
rz(-0.66508355) q[1];
sx q[1];
rz(2.6159899) q[1];
rz(-pi) q[2];
rz(2.4003827) q[3];
sx q[3];
rz(-1.831372) q[3];
sx q[3];
rz(1.8678566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.7509191) q[2];
sx q[2];
rz(-1.6030739) q[2];
sx q[2];
rz(1.0323367) q[2];
rz(0.8647626) q[3];
sx q[3];
rz(-1.7059749) q[3];
sx q[3];
rz(-0.56149703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.432935) q[0];
sx q[0];
rz(-1.189804) q[0];
sx q[0];
rz(-1.7565961) q[0];
rz(-2.2191091) q[1];
sx q[1];
rz(-1.205516) q[1];
sx q[1];
rz(2.9965957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6330616) q[0];
sx q[0];
rz(-0.2446188) q[0];
sx q[0];
rz(-0.82392719) q[0];
rz(3.112193) q[2];
sx q[2];
rz(-0.85762944) q[2];
sx q[2];
rz(0.61427639) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8783826) q[1];
sx q[1];
rz(-3.0621798) q[1];
sx q[1];
rz(1.9694049) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0078405) q[3];
sx q[3];
rz(-0.87742311) q[3];
sx q[3];
rz(1.0348606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21236913) q[2];
sx q[2];
rz(-2.9085458) q[2];
sx q[2];
rz(1.93335) q[2];
rz(-2.3182747) q[3];
sx q[3];
rz(-1.9602785) q[3];
sx q[3];
rz(1.3594782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.061148297) q[0];
sx q[0];
rz(-0.7875945) q[0];
sx q[0];
rz(1.0795235) q[0];
rz(2.1615255) q[1];
sx q[1];
rz(-1.020224) q[1];
sx q[1];
rz(3.0316839) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0281668) q[0];
sx q[0];
rz(-0.47385212) q[0];
sx q[0];
rz(-0.68606152) q[0];
rz(-pi) q[1];
rz(0.031326913) q[2];
sx q[2];
rz(-0.27849012) q[2];
sx q[2];
rz(0.14360076) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7798063) q[1];
sx q[1];
rz(-2.61269) q[1];
sx q[1];
rz(-2.5216549) q[1];
rz(-pi) q[2];
rz(-0.22447666) q[3];
sx q[3];
rz(-2.0894451) q[3];
sx q[3];
rz(-0.26640893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.93426934) q[2];
sx q[2];
rz(-2.0685652) q[2];
sx q[2];
rz(0.70460021) q[2];
rz(-2.8568824) q[3];
sx q[3];
rz(-0.73211089) q[3];
sx q[3];
rz(-2.708191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0076184) q[0];
sx q[0];
rz(-1.824279) q[0];
sx q[0];
rz(0.91868573) q[0];
rz(-2.6559415) q[1];
sx q[1];
rz(-0.44356569) q[1];
sx q[1];
rz(-1.1007016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2393409) q[0];
sx q[0];
rz(-1.4318236) q[0];
sx q[0];
rz(0.13987565) q[0];
rz(-pi) q[1];
rz(1.2687911) q[2];
sx q[2];
rz(-1.7922316) q[2];
sx q[2];
rz(1.0332274) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0079643) q[1];
sx q[1];
rz(-1.6755662) q[1];
sx q[1];
rz(-1.973335) q[1];
rz(-0.63559182) q[3];
sx q[3];
rz(-1.6219536) q[3];
sx q[3];
rz(-0.0035465586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.47621581) q[2];
sx q[2];
rz(-0.85415888) q[2];
sx q[2];
rz(-0.7594792) q[2];
rz(2.1935943) q[3];
sx q[3];
rz(-1.4576603) q[3];
sx q[3];
rz(1.6215526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043902472) q[0];
sx q[0];
rz(-0.49711415) q[0];
sx q[0];
rz(0.33101606) q[0];
rz(0.76200062) q[1];
sx q[1];
rz(-2.8490729) q[1];
sx q[1];
rz(2.5133572) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2866007) q[0];
sx q[0];
rz(-2.1341265) q[0];
sx q[0];
rz(-1.2391633) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8339001) q[2];
sx q[2];
rz(-1.7950222) q[2];
sx q[2];
rz(2.4196408) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.68220968) q[1];
sx q[1];
rz(-1.2222689) q[1];
sx q[1];
rz(-0.63732707) q[1];
rz(-1.9167494) q[3];
sx q[3];
rz(-1.9613492) q[3];
sx q[3];
rz(1.2929163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1198173) q[2];
sx q[2];
rz(-2.6675318) q[2];
sx q[2];
rz(2.9277756) q[2];
rz(-0.98624936) q[3];
sx q[3];
rz(-1.2912368) q[3];
sx q[3];
rz(-3.0333062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66961359) q[0];
sx q[0];
rz(-1.6983953) q[0];
sx q[0];
rz(-1.4629913) q[0];
rz(-1.9646473) q[1];
sx q[1];
rz(-1.6564449) q[1];
sx q[1];
rz(0.0080531837) q[1];
rz(-0.17433737) q[2];
sx q[2];
rz(-1.871289) q[2];
sx q[2];
rz(1.540779) q[2];
rz(-2.4900548) q[3];
sx q[3];
rz(-1.6350411) q[3];
sx q[3];
rz(-1.7512334) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
