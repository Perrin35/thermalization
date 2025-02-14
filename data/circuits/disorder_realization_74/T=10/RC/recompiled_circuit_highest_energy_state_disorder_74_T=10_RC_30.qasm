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
rz(1.9407152) q[0];
sx q[0];
rz(-2.241029) q[0];
sx q[0];
rz(-2.9086034) q[0];
rz(-1.5268582) q[1];
sx q[1];
rz(-2.0695217) q[1];
sx q[1];
rz(1.1195247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716063) q[0];
sx q[0];
rz(-1.8052088) q[0];
sx q[0];
rz(-3.1283911) q[0];
rz(-2.5801611) q[2];
sx q[2];
rz(-1.560964) q[2];
sx q[2];
rz(1.1947643) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0927375) q[1];
sx q[1];
rz(-1.4953145) q[1];
sx q[1];
rz(1.2176179) q[1];
rz(0.60910881) q[3];
sx q[3];
rz(-1.850633) q[3];
sx q[3];
rz(-0.5718872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.69922525) q[2];
sx q[2];
rz(-1.0590326) q[2];
sx q[2];
rz(-0.11288682) q[2];
rz(-0.26116192) q[3];
sx q[3];
rz(-1.7966725) q[3];
sx q[3];
rz(1.2507218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79171044) q[0];
sx q[0];
rz(-0.078502027) q[0];
sx q[0];
rz(0.054340266) q[0];
rz(2.9229274) q[1];
sx q[1];
rz(-1.668914) q[1];
sx q[1];
rz(-2.7770619) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4512149) q[0];
sx q[0];
rz(-2.0663107) q[0];
sx q[0];
rz(-0.69083237) q[0];
rz(-pi) q[1];
rz(2.7414315) q[2];
sx q[2];
rz(-0.98099835) q[2];
sx q[2];
rz(1.1584182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5082701) q[1];
sx q[1];
rz(-0.62707096) q[1];
sx q[1];
rz(-3.0116664) q[1];
rz(2.0516615) q[3];
sx q[3];
rz(-0.37853795) q[3];
sx q[3];
rz(-2.7906281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9008122) q[2];
sx q[2];
rz(-1.6996926) q[2];
sx q[2];
rz(-1.1085294) q[2];
rz(0.19567868) q[3];
sx q[3];
rz(-1.207573) q[3];
sx q[3];
rz(-1.4669363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3942669) q[0];
sx q[0];
rz(-2.1013923) q[0];
sx q[0];
rz(2.105383) q[0];
rz(0.58468435) q[1];
sx q[1];
rz(-1.6744637) q[1];
sx q[1];
rz(2.5856957) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5591878) q[0];
sx q[0];
rz(-2.3055861) q[0];
sx q[0];
rz(0.57484267) q[0];
rz(-1.6481208) q[2];
sx q[2];
rz(-1.5333042) q[2];
sx q[2];
rz(-0.24399569) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4208918) q[1];
sx q[1];
rz(-1.4399505) q[1];
sx q[1];
rz(-2.2761005) q[1];
x q[2];
rz(-2.6922052) q[3];
sx q[3];
rz(-1.3419749) q[3];
sx q[3];
rz(0.23923161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7613775) q[2];
sx q[2];
rz(-2.9351202) q[2];
sx q[2];
rz(-1.9492487) q[2];
rz(3.0377667) q[3];
sx q[3];
rz(-1.9793648) q[3];
sx q[3];
rz(1.1473568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035606774) q[0];
sx q[0];
rz(-0.62653956) q[0];
sx q[0];
rz(1.4242127) q[0];
rz(-2.4183938) q[1];
sx q[1];
rz(-0.23591787) q[1];
sx q[1];
rz(2.9211488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.047565) q[0];
sx q[0];
rz(-1.2050497) q[0];
sx q[0];
rz(1.0785375) q[0];
x q[1];
rz(0.29517572) q[2];
sx q[2];
rz(-2.1116426) q[2];
sx q[2];
rz(-1.4424272) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.053829642) q[1];
sx q[1];
rz(-1.2150032) q[1];
sx q[1];
rz(-3.0213814) q[1];
x q[2];
rz(-2.3001053) q[3];
sx q[3];
rz(-2.609775) q[3];
sx q[3];
rz(-1.2625265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7580938) q[2];
sx q[2];
rz(-1.8269962) q[2];
sx q[2];
rz(2.626075) q[2];
rz(-0.085144194) q[3];
sx q[3];
rz(-0.65535039) q[3];
sx q[3];
rz(2.3084124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67007095) q[0];
sx q[0];
rz(-2.9586198) q[0];
sx q[0];
rz(-0.66437379) q[0];
rz(0.5270671) q[1];
sx q[1];
rz(-1.138569) q[1];
sx q[1];
rz(-1.1767496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19163469) q[0];
sx q[0];
rz(-1.5427046) q[0];
sx q[0];
rz(-1.4856443) q[0];
rz(-pi) q[1];
rz(0.46026261) q[2];
sx q[2];
rz(-2.2027594) q[2];
sx q[2];
rz(3.021701) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0954353) q[1];
sx q[1];
rz(-1.5514138) q[1];
sx q[1];
rz(2.9033015) q[1];
rz(-pi) q[2];
rz(0.99915766) q[3];
sx q[3];
rz(-1.7222341) q[3];
sx q[3];
rz(-3.0157523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.027017) q[2];
sx q[2];
rz(-0.21275529) q[2];
sx q[2];
rz(0.93811718) q[2];
rz(1.045687) q[3];
sx q[3];
rz(-0.18200471) q[3];
sx q[3];
rz(2.3312881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.806458) q[0];
sx q[0];
rz(-1.9430176) q[0];
sx q[0];
rz(-1.2499811) q[0];
rz(-2.9373923) q[1];
sx q[1];
rz(-2.4649492) q[1];
sx q[1];
rz(2.2883889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1565022) q[0];
sx q[0];
rz(-1.9130052) q[0];
sx q[0];
rz(-1.8746822) q[0];
x q[1];
rz(2.5218042) q[2];
sx q[2];
rz(-2.0896122) q[2];
sx q[2];
rz(-0.74558115) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1386765) q[1];
sx q[1];
rz(-2.2670548) q[1];
sx q[1];
rz(2.584444) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2291992) q[3];
sx q[3];
rz(-1.4506243) q[3];
sx q[3];
rz(-0.12334331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.044067232) q[2];
sx q[2];
rz(-1.7996457) q[2];
sx q[2];
rz(-1.9419144) q[2];
rz(-2.1357644) q[3];
sx q[3];
rz(-0.89322105) q[3];
sx q[3];
rz(0.83542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.47578874) q[0];
sx q[0];
rz(-0.26837334) q[0];
sx q[0];
rz(-0.98261181) q[0];
rz(1.3081374) q[1];
sx q[1];
rz(-1.2160701) q[1];
sx q[1];
rz(1.7024202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2932201) q[0];
sx q[0];
rz(-1.9008027) q[0];
sx q[0];
rz(2.7820935) q[0];
rz(1.4365044) q[2];
sx q[2];
rz(-0.52596131) q[2];
sx q[2];
rz(0.44912042) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3130668) q[1];
sx q[1];
rz(-1.2432352) q[1];
sx q[1];
rz(0.86159535) q[1];
rz(1.7129219) q[3];
sx q[3];
rz(-1.097659) q[3];
sx q[3];
rz(-0.17156916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16284379) q[2];
sx q[2];
rz(-2.2480201) q[2];
sx q[2];
rz(2.5992375) q[2];
rz(2.1584623) q[3];
sx q[3];
rz(-1.1636846) q[3];
sx q[3];
rz(2.2014309) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56383175) q[0];
sx q[0];
rz(-1.7153808) q[0];
sx q[0];
rz(-2.3610709) q[0];
rz(1.966194) q[1];
sx q[1];
rz(-2.5927717) q[1];
sx q[1];
rz(-1.0343879) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18984737) q[0];
sx q[0];
rz(-1.4035514) q[0];
sx q[0];
rz(2.7664004) q[0];
rz(-pi) q[1];
rz(-1.7138359) q[2];
sx q[2];
rz(-1.489822) q[2];
sx q[2];
rz(-1.2317927) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24496291) q[1];
sx q[1];
rz(-1.9707929) q[1];
sx q[1];
rz(-1.6952312) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54776056) q[3];
sx q[3];
rz(-1.3410853) q[3];
sx q[3];
rz(1.3198157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15022755) q[2];
sx q[2];
rz(-1.4069724) q[2];
sx q[2];
rz(3.0312209) q[2];
rz(-1.8433833) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(-0.29236326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89123911) q[0];
sx q[0];
rz(-2.1284916) q[0];
sx q[0];
rz(1.891834) q[0];
rz(-0.091014422) q[1];
sx q[1];
rz(-0.72151557) q[1];
sx q[1];
rz(0.7116085) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7994493) q[0];
sx q[0];
rz(-0.79252386) q[0];
sx q[0];
rz(-1.5802797) q[0];
rz(-pi) q[1];
rz(2.3046012) q[2];
sx q[2];
rz(-2.2185433) q[2];
sx q[2];
rz(-0.68488065) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8636057) q[1];
sx q[1];
rz(-2.5558439) q[1];
sx q[1];
rz(2.4892155) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36104155) q[3];
sx q[3];
rz(-0.78892498) q[3];
sx q[3];
rz(1.5004683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.3332425) q[2];
sx q[2];
rz(-1.658172) q[2];
sx q[2];
rz(-1.1727772) q[2];
rz(1.5277537) q[3];
sx q[3];
rz(-1.0781735) q[3];
sx q[3];
rz(0.095002256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.461819) q[0];
sx q[0];
rz(-0.12400308) q[0];
sx q[0];
rz(1.5754196) q[0];
rz(-0.35789403) q[1];
sx q[1];
rz(-1.0708258) q[1];
sx q[1];
rz(0.90248743) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99011496) q[0];
sx q[0];
rz(-0.95989908) q[0];
sx q[0];
rz(0.31661242) q[0];
rz(-pi) q[1];
rz(-0.029377653) q[2];
sx q[2];
rz(-0.64422078) q[2];
sx q[2];
rz(-0.51960301) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6186699) q[1];
sx q[1];
rz(-0.82397193) q[1];
sx q[1];
rz(-2.0342779) q[1];
rz(-pi) q[2];
rz(-0.087335056) q[3];
sx q[3];
rz(-2.119434) q[3];
sx q[3];
rz(-2.5256796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0365399) q[2];
sx q[2];
rz(-0.5883216) q[2];
sx q[2];
rz(-1.499929) q[2];
rz(0.52122742) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(-2.1874793) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70984107) q[0];
sx q[0];
rz(-0.7702282) q[0];
sx q[0];
rz(-1.7256398) q[0];
rz(0.00090986666) q[1];
sx q[1];
rz(-1.4716499) q[1];
sx q[1];
rz(1.6843527) q[1];
rz(2.1442689) q[2];
sx q[2];
rz(-1.4880585) q[2];
sx q[2];
rz(-2.1941296) q[2];
rz(-2.3933181) q[3];
sx q[3];
rz(-1.0327485) q[3];
sx q[3];
rz(-2.2365981) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
