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
rz(-1.0567226) q[0];
sx q[0];
rz(-1.280008) q[0];
sx q[0];
rz(-0.7260538) q[0];
rz(-2.5441406) q[1];
sx q[1];
rz(-0.51550454) q[1];
sx q[1];
rz(-1.0322303) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9361794) q[0];
sx q[0];
rz(-1.4698615) q[0];
sx q[0];
rz(2.4839782) q[0];
rz(-1.502336) q[2];
sx q[2];
rz(-2.6227544) q[2];
sx q[2];
rz(-1.8585132) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.5121442) q[1];
sx q[1];
rz(-1.2595065) q[1];
sx q[1];
rz(1.4427079) q[1];
x q[2];
rz(1.0606403) q[3];
sx q[3];
rz(-0.83612305) q[3];
sx q[3];
rz(-1.6250798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7008179) q[2];
sx q[2];
rz(-2.5371607) q[2];
sx q[2];
rz(-0.087892858) q[2];
rz(-1.6739738) q[3];
sx q[3];
rz(-1.847495) q[3];
sx q[3];
rz(-2.1427593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1217594) q[0];
sx q[0];
rz(-1.8353945) q[0];
sx q[0];
rz(-1.649296) q[0];
rz(-0.78798931) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(0.96510395) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5934629) q[0];
sx q[0];
rz(-0.735518) q[0];
sx q[0];
rz(3.0043362) q[0];
rz(-pi) q[1];
rz(0.73678686) q[2];
sx q[2];
rz(-1.6766785) q[2];
sx q[2];
rz(2.4491765) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9763042) q[1];
sx q[1];
rz(-1.7563662) q[1];
sx q[1];
rz(-1.1052422) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6613668) q[3];
sx q[3];
rz(-1.1027059) q[3];
sx q[3];
rz(-0.15453574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1284156) q[2];
sx q[2];
rz(-2.8030277) q[2];
sx q[2];
rz(1.3236807) q[2];
rz(0.88661083) q[3];
sx q[3];
rz(-1.1611791) q[3];
sx q[3];
rz(-2.9429341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24432467) q[0];
sx q[0];
rz(-0.97049814) q[0];
sx q[0];
rz(0.6955198) q[0];
rz(2.3059402) q[1];
sx q[1];
rz(-1.7213768) q[1];
sx q[1];
rz(2.4998891) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013896522) q[0];
sx q[0];
rz(-0.75152107) q[0];
sx q[0];
rz(-1.7836545) q[0];
x q[1];
rz(2.9474845) q[2];
sx q[2];
rz(-2.420937) q[2];
sx q[2];
rz(-1.0247165) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1759626) q[1];
sx q[1];
rz(-1.2568736) q[1];
sx q[1];
rz(2.0759275) q[1];
x q[2];
rz(-2.1459177) q[3];
sx q[3];
rz(-0.96700562) q[3];
sx q[3];
rz(-3.1310981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16331638) q[2];
sx q[2];
rz(-1.1288319) q[2];
sx q[2];
rz(-1.3345435) q[2];
rz(-1.4000019) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(2.002772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1432994) q[0];
sx q[0];
rz(-0.87109733) q[0];
sx q[0];
rz(2.368108) q[0];
rz(-3.0553014) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(0.48826826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937259) q[0];
sx q[0];
rz(-2.6742088) q[0];
sx q[0];
rz(-2.9899238) q[0];
rz(-2.1147229) q[2];
sx q[2];
rz(-0.94878093) q[2];
sx q[2];
rz(-2.029947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99566702) q[1];
sx q[1];
rz(-0.85595501) q[1];
sx q[1];
rz(-0.30024528) q[1];
rz(1.9776439) q[3];
sx q[3];
rz(-1.196612) q[3];
sx q[3];
rz(-2.0896623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46828541) q[2];
sx q[2];
rz(-1.1136592) q[2];
sx q[2];
rz(0.16806531) q[2];
rz(0.42260653) q[3];
sx q[3];
rz(-0.97413078) q[3];
sx q[3];
rz(2.9882123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76822686) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(-1.9510829) q[0];
rz(2.6137784) q[1];
sx q[1];
rz(-1.0820791) q[1];
sx q[1];
rz(-0.0091008069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2784257) q[0];
sx q[0];
rz(-0.7896119) q[0];
sx q[0];
rz(-2.4900683) q[0];
rz(-pi) q[1];
rz(2.3528655) q[2];
sx q[2];
rz(-0.085914748) q[2];
sx q[2];
rz(2.9979561) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.02024164) q[1];
sx q[1];
rz(-0.96881908) q[1];
sx q[1];
rz(-0.51687981) q[1];
rz(2.6713773) q[3];
sx q[3];
rz(-2.145916) q[3];
sx q[3];
rz(1.0660389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3147543) q[2];
sx q[2];
rz(-2.645731) q[2];
sx q[2];
rz(0.66519386) q[2];
rz(1.0264171) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(-2.3608666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5196395) q[0];
sx q[0];
rz(-2.6060947) q[0];
sx q[0];
rz(-2.7435379) q[0];
rz(-1.6710501) q[1];
sx q[1];
rz(-0.87659756) q[1];
sx q[1];
rz(-0.7801396) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7791634) q[0];
sx q[0];
rz(-1.0991968) q[0];
sx q[0];
rz(-0.40059025) q[0];
rz(0.69177912) q[2];
sx q[2];
rz(-2.8304812) q[2];
sx q[2];
rz(1.2300613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91785039) q[1];
sx q[1];
rz(-2.3931574) q[1];
sx q[1];
rz(-2.7176932) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1885181) q[3];
sx q[3];
rz(-0.72186379) q[3];
sx q[3];
rz(0.75306276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.857343) q[2];
sx q[2];
rz(-1.1025068) q[2];
sx q[2];
rz(-0.16172376) q[2];
rz(3.0297847) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(-2.4017754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9991456) q[0];
sx q[0];
rz(-1.7113926) q[0];
sx q[0];
rz(-2.4430742) q[0];
rz(1.0922208) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(-0.05365595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2630477) q[0];
sx q[0];
rz(-1.8302655) q[0];
sx q[0];
rz(0.83570133) q[0];
rz(-1.1321819) q[2];
sx q[2];
rz(-1.2234067) q[2];
sx q[2];
rz(1.6288822) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21660575) q[1];
sx q[1];
rz(-2.0408071) q[1];
sx q[1];
rz(-1.5108677) q[1];
rz(1.3732304) q[3];
sx q[3];
rz(-2.4854492) q[3];
sx q[3];
rz(3.0246317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0003537) q[2];
sx q[2];
rz(-2.443479) q[2];
sx q[2];
rz(-0.21256438) q[2];
rz(-2.8138748) q[3];
sx q[3];
rz(-0.98723427) q[3];
sx q[3];
rz(-1.3535961) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.245529) q[0];
sx q[0];
rz(-2.0241757) q[0];
sx q[0];
rz(-2.8118706) q[0];
rz(-0.7715191) q[1];
sx q[1];
rz(-0.78116575) q[1];
sx q[1];
rz(2.3158997) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2742608) q[0];
sx q[0];
rz(-1.6578339) q[0];
sx q[0];
rz(0.66477338) q[0];
rz(-pi) q[1];
rz(1.4554899) q[2];
sx q[2];
rz(-2.5421511) q[2];
sx q[2];
rz(1.2931461) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.030071478) q[1];
sx q[1];
rz(-2.3276405) q[1];
sx q[1];
rz(1.7205283) q[1];
x q[2];
rz(-0.25715526) q[3];
sx q[3];
rz(-1.7622593) q[3];
sx q[3];
rz(1.2068656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9416435) q[2];
sx q[2];
rz(-1.7295094) q[2];
sx q[2];
rz(1.6597718) q[2];
rz(-0.09305772) q[3];
sx q[3];
rz(-1.2991644) q[3];
sx q[3];
rz(0.53946462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.900443) q[0];
sx q[0];
rz(-0.41541442) q[0];
sx q[0];
rz(2.7442617) q[0];
rz(-2.1516402) q[1];
sx q[1];
rz(-1.7918469) q[1];
sx q[1];
rz(-1.0362524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66175211) q[0];
sx q[0];
rz(-1.1660728) q[0];
sx q[0];
rz(-1.5581801) q[0];
x q[1];
rz(-0.36340275) q[2];
sx q[2];
rz(-0.37523233) q[2];
sx q[2];
rz(0.64363942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4084928) q[1];
sx q[1];
rz(-1.6599791) q[1];
sx q[1];
rz(2.9217138) q[1];
rz(-0.37452368) q[3];
sx q[3];
rz(-2.0595042) q[3];
sx q[3];
rz(2.599409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.022126023) q[2];
sx q[2];
rz(-0.72212044) q[2];
sx q[2];
rz(-1.3163346) q[2];
rz(0.25257603) q[3];
sx q[3];
rz(-1.3467237) q[3];
sx q[3];
rz(-2.778497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13755688) q[0];
sx q[0];
rz(-0.79020774) q[0];
sx q[0];
rz(-0.23605119) q[0];
rz(-0.099418489) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(-1.258446) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6981068) q[0];
sx q[0];
rz(-2.5763566) q[0];
sx q[0];
rz(2.9960669) q[0];
rz(-pi) q[1];
rz(0.66345556) q[2];
sx q[2];
rz(-1.1642312) q[2];
sx q[2];
rz(2.6321509) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8736781) q[1];
sx q[1];
rz(-0.89338747) q[1];
sx q[1];
rz(-1.7180841) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.997273) q[3];
sx q[3];
rz(-1.9460554) q[3];
sx q[3];
rz(0.29510185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1666169) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(-2.100259) q[2];
rz(0.58639041) q[3];
sx q[3];
rz(-0.47724637) q[3];
sx q[3];
rz(0.81048107) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7753684) q[0];
sx q[0];
rz(-2.0851705) q[0];
sx q[0];
rz(2.0023517) q[0];
rz(-2.1495023) q[1];
sx q[1];
rz(-1.6012259) q[1];
sx q[1];
rz(1.59902) q[1];
rz(0.31789657) q[2];
sx q[2];
rz(-1.3513202) q[2];
sx q[2];
rz(2.8693595) q[2];
rz(0.59320368) q[3];
sx q[3];
rz(-1.2013669) q[3];
sx q[3];
rz(-0.47040924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
