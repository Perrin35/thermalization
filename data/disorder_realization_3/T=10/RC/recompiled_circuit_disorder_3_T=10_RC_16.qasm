OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(2.7273942) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1026693) q[0];
sx q[0];
rz(-1.2352714) q[0];
sx q[0];
rz(-1.194792) q[0];
rz(-pi) q[1];
rz(-1.0010927) q[2];
sx q[2];
rz(-1.4909407) q[2];
sx q[2];
rz(-0.18730883) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84250433) q[1];
sx q[1];
rz(-0.43049225) q[1];
sx q[1];
rz(-1.4166142) q[1];
rz(-1.6730509) q[3];
sx q[3];
rz(-1.3073321) q[3];
sx q[3];
rz(1.9139569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(-3.1100173) q[2];
rz(-1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3027705) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-0.1698499) q[0];
rz(0.70392144) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(0.53952113) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44885143) q[0];
sx q[0];
rz(-0.451938) q[0];
sx q[0];
rz(-1.6777722) q[0];
rz(-2.6846634) q[2];
sx q[2];
rz(-0.77605844) q[2];
sx q[2];
rz(1.6844144) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5121465) q[1];
sx q[1];
rz(-2.0680032) q[1];
sx q[1];
rz(2.9260103) q[1];
rz(-pi) q[2];
rz(2.5856421) q[3];
sx q[3];
rz(-0.9405989) q[3];
sx q[3];
rz(2.6549908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66723055) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-0.24307069) q[2];
rz(-0.66611755) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3617525) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(-2.6932122) q[0];
rz(-1.7547296) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(0.2562491) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6305144) q[0];
sx q[0];
rz(-0.96195463) q[0];
sx q[0];
rz(1.8947253) q[0];
x q[1];
rz(2.3149812) q[2];
sx q[2];
rz(-0.6140784) q[2];
sx q[2];
rz(0.6073063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5445404) q[1];
sx q[1];
rz(-1.8982732) q[1];
sx q[1];
rz(-2.271133) q[1];
x q[2];
rz(2.9647397) q[3];
sx q[3];
rz(-1.6079418) q[3];
sx q[3];
rz(2.0166486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8024575) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(2.5668872) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9280076) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(-2.8821049) q[0];
rz(-1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(0.73192275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5314732) q[0];
sx q[0];
rz(-1.1024794) q[0];
sx q[0];
rz(2.5584695) q[0];
rz(-pi) q[1];
rz(-1.6875793) q[2];
sx q[2];
rz(-0.8568961) q[2];
sx q[2];
rz(1.8494948) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4983474) q[1];
sx q[1];
rz(-0.32532641) q[1];
sx q[1];
rz(2.494032) q[1];
rz(-pi) q[2];
rz(-0.22051208) q[3];
sx q[3];
rz(-1.6933428) q[3];
sx q[3];
rz(-1.3361271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.2735294) q[2];
sx q[2];
rz(0.15110061) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(0.37724075) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.8792101) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(0.79777065) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7933465) q[0];
sx q[0];
rz(-0.12013809) q[0];
sx q[0];
rz(-1.5500463) q[0];
rz(-pi) q[1];
rz(0.30349489) q[2];
sx q[2];
rz(-1.7804838) q[2];
sx q[2];
rz(-1.4022624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3623912) q[1];
sx q[1];
rz(-2.5070842) q[1];
sx q[1];
rz(1.4125376) q[1];
rz(-pi) q[2];
rz(0.4425211) q[3];
sx q[3];
rz(-1.2708775) q[3];
sx q[3];
rz(-0.72344852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34565869) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(-0.30203715) q[2];
rz(-1.9942412) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(-2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(1.1557895) q[0];
rz(-2.0571158) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(0.070080431) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87175831) q[0];
sx q[0];
rz(-1.7239128) q[0];
sx q[0];
rz(1.7437177) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78380615) q[2];
sx q[2];
rz(-1.1597826) q[2];
sx q[2];
rz(2.59936) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.49680432) q[1];
sx q[1];
rz(-1.7323238) q[1];
sx q[1];
rz(-2.5700388) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.035404215) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8391116) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(0.62136674) q[2];
rz(1.7012043) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(-2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(0.0079356114) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.313109) q[0];
sx q[0];
rz(-2.5223753) q[0];
sx q[0];
rz(-2.6909268) q[0];
rz(-pi) q[1];
rz(-0.87848778) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(-1.0483339) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.24212813) q[1];
sx q[1];
rz(-2.2367396) q[1];
sx q[1];
rz(-0.66023402) q[1];
rz(-pi) q[2];
rz(-1.1778529) q[3];
sx q[3];
rz(-1.279139) q[3];
sx q[3];
rz(-1.1100811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(0.41637862) q[2];
rz(1.773206) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798582) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(0.36488786) q[0];
rz(-2.2015613) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(-1.5023124) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68740326) q[0];
sx q[0];
rz(-2.8912376) q[0];
sx q[0];
rz(-3.1012015) q[0];
rz(-0.47183581) q[2];
sx q[2];
rz(-0.89315692) q[2];
sx q[2];
rz(-0.0080646947) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8870526) q[1];
sx q[1];
rz(-2.2215543) q[1];
sx q[1];
rz(0.73144967) q[1];
rz(0.092076093) q[3];
sx q[3];
rz(-1.3658938) q[3];
sx q[3];
rz(-2.3931707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(2.0765182) q[2];
rz(2.8403357) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(1.7565953) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(2.7246144) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96796658) q[0];
sx q[0];
rz(-0.91714232) q[0];
sx q[0];
rz(3.0941512) q[0];
x q[1];
rz(-2.4531104) q[2];
sx q[2];
rz(-2.3496029) q[2];
sx q[2];
rz(-2.288523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7855362) q[1];
sx q[1];
rz(-2.0691263) q[1];
sx q[1];
rz(0.15798012) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7008408) q[3];
sx q[3];
rz(-1.4308883) q[3];
sx q[3];
rz(0.035294447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10432648) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(-2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41480961) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(-1.6171932) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(1.8189925) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37227092) q[0];
sx q[0];
rz(-1.1936545) q[0];
sx q[0];
rz(-3.0110714) q[0];
rz(-pi) q[1];
rz(2.0199213) q[2];
sx q[2];
rz(-1.3438864) q[2];
sx q[2];
rz(0.34044468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9433371) q[1];
sx q[1];
rz(-2.2554734) q[1];
sx q[1];
rz(-0.68590045) q[1];
rz(-pi) q[2];
rz(0.15354746) q[3];
sx q[3];
rz(-2.2205177) q[3];
sx q[3];
rz(-0.49680199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3502729) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(1.0160758) q[2];
rz(1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(-0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9983457) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(-1.3636419) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(3.1238363) q[2];
sx q[2];
rz(-0.50928558) q[2];
sx q[2];
rz(-1.7094564) q[2];
rz(-3.0512814) q[3];
sx q[3];
rz(-1.3406546) q[3];
sx q[3];
rz(-1.8484074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
