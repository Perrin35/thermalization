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
rz(-1.2008774) q[0];
sx q[0];
rz(5.3826217) q[0];
sx q[0];
rz(9.1917888) q[0];
rz(-1.5268582) q[1];
sx q[1];
rz(4.2136636) q[1];
sx q[1];
rz(16.827488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4716063) q[0];
sx q[0];
rz(-1.3363839) q[0];
sx q[0];
rz(-3.1283911) q[0];
x q[1];
rz(0.56143151) q[2];
sx q[2];
rz(-1.5806287) q[2];
sx q[2];
rz(1.9468284) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0488551) q[1];
sx q[1];
rz(-1.4953145) q[1];
sx q[1];
rz(-1.9239747) q[1];
rz(-1.2337716) q[3];
sx q[3];
rz(-0.98857461) q[3];
sx q[3];
rz(1.189251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69922525) q[2];
sx q[2];
rz(-2.0825601) q[2];
sx q[2];
rz(3.0287058) q[2];
rz(-0.26116192) q[3];
sx q[3];
rz(-1.3449202) q[3];
sx q[3];
rz(-1.2507218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.79171044) q[0];
sx q[0];
rz(-3.0630906) q[0];
sx q[0];
rz(-0.054340266) q[0];
rz(-0.21866523) q[1];
sx q[1];
rz(-1.4726787) q[1];
sx q[1];
rz(2.7770619) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35848356) q[0];
sx q[0];
rz(-2.3158584) q[0];
sx q[0];
rz(-0.70348212) q[0];
x q[1];
rz(0.40016115) q[2];
sx q[2];
rz(-2.1605943) q[2];
sx q[2];
rz(1.1584182) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4733401) q[1];
sx q[1];
rz(-2.1917731) q[1];
sx q[1];
rz(-1.4771858) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.909799) q[3];
sx q[3];
rz(-1.7425797) q[3];
sx q[3];
rz(-1.671227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24078044) q[2];
sx q[2];
rz(-1.4419) q[2];
sx q[2];
rz(2.0330632) q[2];
rz(-2.945914) q[3];
sx q[3];
rz(-1.207573) q[3];
sx q[3];
rz(-1.4669363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.3942669) q[0];
sx q[0];
rz(-1.0402004) q[0];
sx q[0];
rz(-2.105383) q[0];
rz(-0.58468435) q[1];
sx q[1];
rz(-1.6744637) q[1];
sx q[1];
rz(-2.5856957) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7202458) q[0];
sx q[0];
rz(-1.1555536) q[0];
sx q[0];
rz(2.3930727) q[0];
x q[1];
rz(-2.0228407) q[2];
sx q[2];
rz(-0.085918203) q[2];
sx q[2];
rz(-2.265386) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8396254) q[1];
sx q[1];
rz(-0.71528597) q[1];
sx q[1];
rz(1.3705181) q[1];
x q[2];
rz(2.6494637) q[3];
sx q[3];
rz(-2.6408697) q[3];
sx q[3];
rz(2.2496417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7613775) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(1.9492487) q[2];
rz(-0.10382593) q[3];
sx q[3];
rz(-1.1622279) q[3];
sx q[3];
rz(1.9942358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035606774) q[0];
sx q[0];
rz(-2.5150531) q[0];
sx q[0];
rz(1.71738) q[0];
rz(-2.4183938) q[1];
sx q[1];
rz(-2.9056748) q[1];
sx q[1];
rz(-2.9211488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.047565) q[0];
sx q[0];
rz(-1.936543) q[0];
sx q[0];
rz(2.0630552) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1197243) q[2];
sx q[2];
rz(-2.5325518) q[2];
sx q[2];
rz(2.232617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.087763) q[1];
sx q[1];
rz(-1.2150032) q[1];
sx q[1];
rz(-0.1202113) q[1];
rz(0.8414874) q[3];
sx q[3];
rz(-2.609775) q[3];
sx q[3];
rz(1.8790661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38349884) q[2];
sx q[2];
rz(-1.8269962) q[2];
sx q[2];
rz(0.51551762) q[2];
rz(-0.085144194) q[3];
sx q[3];
rz(-2.4862423) q[3];
sx q[3];
rz(0.83318025) q[3];
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
rz(-0.67007095) q[0];
sx q[0];
rz(-2.9586198) q[0];
sx q[0];
rz(2.4772189) q[0];
rz(-0.5270671) q[1];
sx q[1];
rz(-1.138569) q[1];
sx q[1];
rz(1.1767496) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7648286) q[0];
sx q[0];
rz(-1.485678) q[0];
sx q[0];
rz(0.028193817) q[0];
rz(-pi) q[1];
rz(-2.68133) q[2];
sx q[2];
rz(-2.2027594) q[2];
sx q[2];
rz(3.021701) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.55498153) q[1];
sx q[1];
rz(-2.9025295) q[1];
sx q[1];
rz(3.0596517) q[1];
x q[2];
rz(1.2958636) q[3];
sx q[3];
rz(-0.58919135) q[3];
sx q[3];
rz(-1.9269892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.027017) q[2];
sx q[2];
rz(-2.9288374) q[2];
sx q[2];
rz(-0.93811718) q[2];
rz(2.0959057) q[3];
sx q[3];
rz(-0.18200471) q[3];
sx q[3];
rz(0.81030455) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.806458) q[0];
sx q[0];
rz(-1.198575) q[0];
sx q[0];
rz(1.8916116) q[0];
rz(0.20420034) q[1];
sx q[1];
rz(-2.4649492) q[1];
sx q[1];
rz(2.2883889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51914224) q[0];
sx q[0];
rz(-1.8565489) q[0];
sx q[0];
rz(0.35730548) q[0];
rz(-pi) q[1];
rz(2.3647652) q[2];
sx q[2];
rz(-0.78571253) q[2];
sx q[2];
rz(1.7091027) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.051999747) q[1];
sx q[1];
rz(-1.9885364) q[1];
sx q[1];
rz(-0.79302782) q[1];
rz(-pi) q[2];
rz(-1.2291992) q[3];
sx q[3];
rz(-1.6909684) q[3];
sx q[3];
rz(-0.12334331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0975254) q[2];
sx q[2];
rz(-1.3419469) q[2];
sx q[2];
rz(-1.1996783) q[2];
rz(1.0058282) q[3];
sx q[3];
rz(-0.89322105) q[3];
sx q[3];
rz(0.83542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47578874) q[0];
sx q[0];
rz(-0.26837334) q[0];
sx q[0];
rz(-0.98261181) q[0];
rz(-1.8334552) q[1];
sx q[1];
rz(-1.9255226) q[1];
sx q[1];
rz(-1.7024202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7076232) q[0];
sx q[0];
rz(-0.48309193) q[0];
sx q[0];
rz(-0.77204319) q[0];
rz(-pi) q[1];
rz(1.4365044) q[2];
sx q[2];
rz(-2.6156313) q[2];
sx q[2];
rz(-0.44912042) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.01165302) q[1];
sx q[1];
rz(-0.90637744) q[1];
sx q[1];
rz(0.42099492) q[1];
x q[2];
rz(-2.6643326) q[3];
sx q[3];
rz(-1.4443732) q[3];
sx q[3];
rz(1.8074769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9787489) q[2];
sx q[2];
rz(-0.89357251) q[2];
sx q[2];
rz(2.5992375) q[2];
rz(-2.1584623) q[3];
sx q[3];
rz(-1.977908) q[3];
sx q[3];
rz(-0.94016176) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56383175) q[0];
sx q[0];
rz(-1.4262119) q[0];
sx q[0];
rz(-2.3610709) q[0];
rz(1.1753987) q[1];
sx q[1];
rz(-2.5927717) q[1];
sx q[1];
rz(-2.1072047) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808409) q[0];
sx q[0];
rz(-0.409161) q[0];
sx q[0];
rz(-2.7098814) q[0];
rz(-pi) q[1];
rz(-3.0597866) q[2];
sx q[2];
rz(-1.4282287) q[2];
sx q[2];
rz(-2.7909401) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0757622) q[1];
sx q[1];
rz(-2.7236863) q[1];
sx q[1];
rz(2.8560545) q[1];
rz(2.5938321) q[3];
sx q[3];
rz(-1.8005074) q[3];
sx q[3];
rz(-1.3198157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9913651) q[2];
sx q[2];
rz(-1.4069724) q[2];
sx q[2];
rz(3.0312209) q[2];
rz(1.8433833) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(-2.8492294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89123911) q[0];
sx q[0];
rz(-2.1284916) q[0];
sx q[0];
rz(-1.891834) q[0];
rz(-3.0505782) q[1];
sx q[1];
rz(-2.4200771) q[1];
sx q[1];
rz(-2.4299842) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34214333) q[0];
sx q[0];
rz(-2.3490688) q[0];
sx q[0];
rz(1.5613129) q[0];
rz(-2.3046012) q[2];
sx q[2];
rz(-2.2185433) q[2];
sx q[2];
rz(0.68488065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8597128) q[1];
sx q[1];
rz(-1.9130453) q[1];
sx q[1];
rz(0.48515353) q[1];
rz(0.75562128) q[3];
sx q[3];
rz(-1.3174302) q[3];
sx q[3];
rz(2.8112335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8083501) q[2];
sx q[2];
rz(-1.4834206) q[2];
sx q[2];
rz(-1.9688155) q[2];
rz(-1.6138389) q[3];
sx q[3];
rz(-1.0781735) q[3];
sx q[3];
rz(-3.0465904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.461819) q[0];
sx q[0];
rz(-3.0175896) q[0];
sx q[0];
rz(1.5754196) q[0];
rz(0.35789403) q[1];
sx q[1];
rz(-2.0707668) q[1];
sx q[1];
rz(-2.2391052) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47115147) q[0];
sx q[0];
rz(-0.67865279) q[0];
sx q[0];
rz(1.1524234) q[0];
rz(-0.029377653) q[2];
sx q[2];
rz(-2.4973719) q[2];
sx q[2];
rz(-2.6219896) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1572307) q[1];
sx q[1];
rz(-2.2868726) q[1];
sx q[1];
rz(-2.6916719) q[1];
rz(-pi) q[2];
rz(3.0542576) q[3];
sx q[3];
rz(-1.0221586) q[3];
sx q[3];
rz(2.5256796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0365399) q[2];
sx q[2];
rz(-2.5532711) q[2];
sx q[2];
rz(-1.6416637) q[2];
rz(2.6203652) q[3];
sx q[3];
rz(-2.9977377) q[3];
sx q[3];
rz(0.95411333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70984107) q[0];
sx q[0];
rz(-2.3713645) q[0];
sx q[0];
rz(1.4159528) q[0];
rz(-0.00090986666) q[1];
sx q[1];
rz(-1.6699427) q[1];
sx q[1];
rz(-1.4572399) q[1];
rz(3.0431912) q[2];
sx q[2];
rz(-0.99953166) q[2];
sx q[2];
rz(2.4649323) q[2];
rz(2.2542027) q[3];
sx q[3];
rz(-2.1947464) q[3];
sx q[3];
rz(-1.1098679) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
