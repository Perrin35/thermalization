OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9031653) q[0];
sx q[0];
rz(2.981346) q[0];
sx q[0];
rz(7.5709406) q[0];
rz(-0.84914485) q[1];
sx q[1];
rz(-1.1916817) q[1];
sx q[1];
rz(-2.3828659) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1232226) q[0];
sx q[0];
rz(-1.8022853) q[0];
sx q[0];
rz(-2.8442848) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8992462) q[2];
sx q[2];
rz(-2.0449315) q[2];
sx q[2];
rz(-0.15210064) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10630305) q[1];
sx q[1];
rz(-0.83005164) q[1];
sx q[1];
rz(-1.071005) q[1];
rz(-2.4843744) q[3];
sx q[3];
rz(-2.2003678) q[3];
sx q[3];
rz(2.8384125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3005001) q[2];
sx q[2];
rz(-2.3427343) q[2];
sx q[2];
rz(0.50660261) q[2];
rz(1.6182342) q[3];
sx q[3];
rz(-2.5169499) q[3];
sx q[3];
rz(-3.0860331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8697206) q[0];
sx q[0];
rz(-2.6756918) q[0];
sx q[0];
rz(-0.063657612) q[0];
rz(1.9438538) q[1];
sx q[1];
rz(-1.5143737) q[1];
sx q[1];
rz(-2.1228085) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96589027) q[0];
sx q[0];
rz(-1.7164443) q[0];
sx q[0];
rz(-2.9431828) q[0];
rz(-pi) q[1];
rz(-1.8070142) q[2];
sx q[2];
rz(-1.2389606) q[2];
sx q[2];
rz(2.163909) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6117759) q[1];
sx q[1];
rz(-0.88068286) q[1];
sx q[1];
rz(-1.9326831) q[1];
x q[2];
rz(-1.7131931) q[3];
sx q[3];
rz(-2.3835858) q[3];
sx q[3];
rz(2.6482794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7055052) q[2];
sx q[2];
rz(-0.59385308) q[2];
sx q[2];
rz(-2.2994821) q[2];
rz(-2.0841133) q[3];
sx q[3];
rz(-2.2130241) q[3];
sx q[3];
rz(-1.1697945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975526) q[0];
sx q[0];
rz(-2.121448) q[0];
sx q[0];
rz(1.1919588) q[0];
rz(-0.92487088) q[1];
sx q[1];
rz(-2.4332739) q[1];
sx q[1];
rz(-1.301773) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0353006) q[0];
sx q[0];
rz(-2.3791693) q[0];
sx q[0];
rz(-0.85233217) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19916735) q[2];
sx q[2];
rz(-1.5323497) q[2];
sx q[2];
rz(-0.13800114) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5980581) q[1];
sx q[1];
rz(-0.50331668) q[1];
sx q[1];
rz(2.7642168) q[1];
rz(-pi) q[2];
rz(1.916094) q[3];
sx q[3];
rz(-1.32816) q[3];
sx q[3];
rz(2.7079328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3785582) q[2];
sx q[2];
rz(-3.0578461) q[2];
sx q[2];
rz(-0.087510022) q[2];
rz(-1.3148426) q[3];
sx q[3];
rz(-2.0851236) q[3];
sx q[3];
rz(-2.1480613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4104376) q[0];
sx q[0];
rz(-1.1719828) q[0];
sx q[0];
rz(-0.80899578) q[0];
rz(2.1058829) q[1];
sx q[1];
rz(-2.6782942) q[1];
sx q[1];
rz(2.1083924) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9161578) q[0];
sx q[0];
rz(-1.5317129) q[0];
sx q[0];
rz(-2.2351511) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1076009) q[2];
sx q[2];
rz(-1.189917) q[2];
sx q[2];
rz(-1.2950667) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7455427) q[1];
sx q[1];
rz(-1.1078568) q[1];
sx q[1];
rz(1.3697423) q[1];
rz(-pi) q[2];
rz(-0.54368005) q[3];
sx q[3];
rz(-1.1515695) q[3];
sx q[3];
rz(2.4329684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0674627) q[2];
sx q[2];
rz(-0.16572696) q[2];
sx q[2];
rz(-1.6039675) q[2];
rz(-1.7240723) q[3];
sx q[3];
rz(-2.0310903) q[3];
sx q[3];
rz(2.0541151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23987016) q[0];
sx q[0];
rz(-2.9481695) q[0];
sx q[0];
rz(-1.9523917) q[0];
rz(2.4173648) q[1];
sx q[1];
rz(-1.8502356) q[1];
sx q[1];
rz(1.474818) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3132202) q[0];
sx q[0];
rz(-0.68183866) q[0];
sx q[0];
rz(2.4260055) q[0];
x q[1];
rz(2.3442845) q[2];
sx q[2];
rz(-1.458079) q[2];
sx q[2];
rz(-1.1657451) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2853299) q[1];
sx q[1];
rz(-1.26795) q[1];
sx q[1];
rz(-1.21638) q[1];
rz(-pi) q[2];
rz(-0.83051564) q[3];
sx q[3];
rz(-0.6336824) q[3];
sx q[3];
rz(-0.049402852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7597947) q[2];
sx q[2];
rz(-0.90961027) q[2];
sx q[2];
rz(-2.6749715) q[2];
rz(-1.4383379) q[3];
sx q[3];
rz(-1.2983026) q[3];
sx q[3];
rz(-2.2013825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71317116) q[0];
sx q[0];
rz(-2.3155825) q[0];
sx q[0];
rz(3.1312422) q[0];
rz(-1.5230644) q[1];
sx q[1];
rz(-1.3925939) q[1];
sx q[1];
rz(-1.3759605) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3165163) q[0];
sx q[0];
rz(-2.7211271) q[0];
sx q[0];
rz(-2.1964873) q[0];
rz(-0.59734663) q[2];
sx q[2];
rz(-0.22946363) q[2];
sx q[2];
rz(0.84537431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.350647) q[1];
sx q[1];
rz(-2.6058491) q[1];
sx q[1];
rz(-0.36046268) q[1];
rz(-pi) q[2];
rz(0.74919219) q[3];
sx q[3];
rz(-1.0674879) q[3];
sx q[3];
rz(2.1319413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4202262) q[2];
sx q[2];
rz(-1.0488291) q[2];
sx q[2];
rz(-2.4143207) q[2];
rz(-0.51491245) q[3];
sx q[3];
rz(-1.3313096) q[3];
sx q[3];
rz(-1.7236727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(0.10776831) q[0];
sx q[0];
rz(-1.6078147) q[0];
sx q[0];
rz(-2.7346101) q[0];
rz(-2.175323) q[1];
sx q[1];
rz(-2.2844908) q[1];
sx q[1];
rz(0.13872096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4086558) q[0];
sx q[0];
rz(-1.6870903) q[0];
sx q[0];
rz(0.035712042) q[0];
x q[1];
rz(-0.27757711) q[2];
sx q[2];
rz(-1.7834366) q[2];
sx q[2];
rz(2.9523498) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2695864) q[1];
sx q[1];
rz(-1.2254189) q[1];
sx q[1];
rz(-0.38265444) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2267031) q[3];
sx q[3];
rz(-1.7100681) q[3];
sx q[3];
rz(0.95668281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8653284) q[2];
sx q[2];
rz(-1.3778957) q[2];
sx q[2];
rz(-2.830937) q[2];
rz(1.3079414) q[3];
sx q[3];
rz(-1.5111978) q[3];
sx q[3];
rz(-0.37180296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0887611) q[0];
sx q[0];
rz(-0.3449057) q[0];
sx q[0];
rz(-0.088223591) q[0];
rz(3.131033) q[1];
sx q[1];
rz(-1.6190395) q[1];
sx q[1];
rz(-0.47058502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7475766) q[0];
sx q[0];
rz(-2.2758099) q[0];
sx q[0];
rz(0.88253077) q[0];
rz(-pi) q[1];
rz(-2.9520985) q[2];
sx q[2];
rz(-0.81606401) q[2];
sx q[2];
rz(0.59743728) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6424017) q[1];
sx q[1];
rz(-2.8676053) q[1];
sx q[1];
rz(1.2528595) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0017774) q[3];
sx q[3];
rz(-1.2212409) q[3];
sx q[3];
rz(2.4619554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18870246) q[2];
sx q[2];
rz(-2.1850093) q[2];
sx q[2];
rz(1.0799705) q[2];
rz(-2.2522669) q[3];
sx q[3];
rz(-1.5998799) q[3];
sx q[3];
rz(-1.5573474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66460669) q[0];
sx q[0];
rz(-2.7598858) q[0];
sx q[0];
rz(-0.81740022) q[0];
rz(2.3951702) q[1];
sx q[1];
rz(-0.67459977) q[1];
sx q[1];
rz(-0.93961632) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8259807) q[0];
sx q[0];
rz(-1.570317) q[0];
sx q[0];
rz(1.5879244) q[0];
x q[1];
rz(-0.65931658) q[2];
sx q[2];
rz(-0.43956471) q[2];
sx q[2];
rz(-2.2826113) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92781237) q[1];
sx q[1];
rz(-0.27691388) q[1];
sx q[1];
rz(2.0109948) q[1];
rz(-pi) q[2];
rz(-1.1577179) q[3];
sx q[3];
rz(-1.7677757) q[3];
sx q[3];
rz(2.911681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32999906) q[2];
sx q[2];
rz(-2.1758695) q[2];
sx q[2];
rz(2.7962371) q[2];
rz(1.6164448) q[3];
sx q[3];
rz(-1.350178) q[3];
sx q[3];
rz(-0.42720544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78373194) q[0];
sx q[0];
rz(-1.1564199) q[0];
sx q[0];
rz(-0.082948908) q[0];
rz(2.9792765) q[1];
sx q[1];
rz(-1.6971308) q[1];
sx q[1];
rz(-2.4457795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7617924) q[0];
sx q[0];
rz(-1.6125935) q[0];
sx q[0];
rz(-0.27970741) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7000574) q[2];
sx q[2];
rz(-1.5737572) q[2];
sx q[2];
rz(-1.2128613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3491427) q[1];
sx q[1];
rz(-2.6969243) q[1];
sx q[1];
rz(-0.044388958) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8105177) q[3];
sx q[3];
rz(-2.4088915) q[3];
sx q[3];
rz(2.1520881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42056981) q[2];
sx q[2];
rz(-2.9837954) q[2];
sx q[2];
rz(-1.921462) q[2];
rz(-0.57341352) q[3];
sx q[3];
rz(-1.3308595) q[3];
sx q[3];
rz(2.8964608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86380105) q[0];
sx q[0];
rz(-1.5024804) q[0];
sx q[0];
rz(-0.13308751) q[0];
rz(0.0028903891) q[1];
sx q[1];
rz(-0.4160226) q[1];
sx q[1];
rz(-0.70152534) q[1];
rz(-0.1724381) q[2];
sx q[2];
rz(-1.3256831) q[2];
sx q[2];
rz(-0.55064461) q[2];
rz(1.5118128) q[3];
sx q[3];
rz(-0.77131693) q[3];
sx q[3];
rz(-0.49280096) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
