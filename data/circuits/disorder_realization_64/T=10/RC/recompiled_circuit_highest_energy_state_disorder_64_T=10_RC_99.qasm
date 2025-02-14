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
rz(2.0848701) q[0];
sx q[0];
rz(-1.8615847) q[0];
sx q[0];
rz(-2.4155389) q[0];
rz(0.59745204) q[1];
sx q[1];
rz(-2.6260881) q[1];
sx q[1];
rz(-2.1093624) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3772349) q[0];
sx q[0];
rz(-0.66417686) q[0];
sx q[0];
rz(2.9773877) q[0];
x q[1];
rz(1.6392567) q[2];
sx q[2];
rz(-2.6227544) q[2];
sx q[2];
rz(1.2830795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6294485) q[1];
sx q[1];
rz(-1.8820861) q[1];
sx q[1];
rz(-1.4427079) q[1];
rz(-pi) q[2];
rz(0.49555669) q[3];
sx q[3];
rz(-2.2751788) q[3];
sx q[3];
rz(0.92951194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44077474) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(-3.0536998) q[2];
rz(1.6739738) q[3];
sx q[3];
rz(-1.847495) q[3];
sx q[3];
rz(2.1427593) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1217594) q[0];
sx q[0];
rz(-1.3061981) q[0];
sx q[0];
rz(-1.4922967) q[0];
rz(0.78798931) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(-0.96510395) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3639587) q[0];
sx q[0];
rz(-2.2978362) q[0];
sx q[0];
rz(-1.6939837) q[0];
rz(-pi) q[1];
rz(-2.4048058) q[2];
sx q[2];
rz(-1.6766785) q[2];
sx q[2];
rz(2.4491765) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75747672) q[1];
sx q[1];
rz(-2.6429477) q[1];
sx q[1];
rz(1.1747141) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6718465) q[3];
sx q[3];
rz(-1.6516017) q[3];
sx q[3];
rz(1.684379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.013177055) q[2];
sx q[2];
rz(-2.8030277) q[2];
sx q[2];
rz(1.817912) q[2];
rz(-0.88661083) q[3];
sx q[3];
rz(-1.1611791) q[3];
sx q[3];
rz(2.9429341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897268) q[0];
sx q[0];
rz(-0.97049814) q[0];
sx q[0];
rz(2.4460728) q[0];
rz(0.8356525) q[1];
sx q[1];
rz(-1.7213768) q[1];
sx q[1];
rz(-2.4998891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27370975) q[0];
sx q[0];
rz(-2.3014304) q[0];
sx q[0];
rz(-0.19489906) q[0];
x q[1];
rz(-1.738606) q[2];
sx q[2];
rz(-0.8665167) q[2];
sx q[2];
rz(-1.2806127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.96563) q[1];
sx q[1];
rz(-1.8847191) q[1];
sx q[1];
rz(-1.0656652) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4738181) q[3];
sx q[3];
rz(-2.3333356) q[3];
sx q[3];
rz(2.2796352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16331638) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(-1.3345435) q[2];
rz(1.4000019) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(-2.002772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1432994) q[0];
sx q[0];
rz(-2.2704953) q[0];
sx q[0];
rz(-2.368108) q[0];
rz(-3.0553014) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(-2.6533244) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047866743) q[0];
sx q[0];
rz(-2.6742088) q[0];
sx q[0];
rz(2.9899238) q[0];
x q[1];
rz(2.4441798) q[2];
sx q[2];
rz(-2.0048755) q[2];
sx q[2];
rz(3.0212439) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.99566702) q[1];
sx q[1];
rz(-0.85595501) q[1];
sx q[1];
rz(-2.8413474) q[1];
rz(-0.4040603) q[3];
sx q[3];
rz(-1.1935788) q[3];
sx q[3];
rz(2.4665143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6733072) q[2];
sx q[2];
rz(-1.1136592) q[2];
sx q[2];
rz(0.16806531) q[2];
rz(0.42260653) q[3];
sx q[3];
rz(-2.1674619) q[3];
sx q[3];
rz(-2.9882123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76822686) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(1.9510829) q[0];
rz(2.6137784) q[1];
sx q[1];
rz(-2.0595136) q[1];
sx q[1];
rz(-3.1324918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21484387) q[0];
sx q[0];
rz(-2.0159449) q[0];
sx q[0];
rz(2.4657004) q[0];
rz(-2.3528655) q[2];
sx q[2];
rz(-3.0556779) q[2];
sx q[2];
rz(2.9979561) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2796265) q[1];
sx q[1];
rz(-1.990296) q[1];
sx q[1];
rz(2.2395833) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1808257) q[3];
sx q[3];
rz(-0.72566635) q[3];
sx q[3];
rz(-0.31455597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3147543) q[2];
sx q[2];
rz(-0.49586168) q[2];
sx q[2];
rz(2.4763988) q[2];
rz(-1.0264171) q[3];
sx q[3];
rz(-2.0668991) q[3];
sx q[3];
rz(-2.3608666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6219532) q[0];
sx q[0];
rz(-2.6060947) q[0];
sx q[0];
rz(-2.7435379) q[0];
rz(-1.6710501) q[1];
sx q[1];
rz(-0.87659756) q[1];
sx q[1];
rz(2.3614531) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7791634) q[0];
sx q[0];
rz(-1.0991968) q[0];
sx q[0];
rz(2.7410024) q[0];
x q[1];
rz(-1.3684811) q[2];
sx q[2];
rz(-1.3328119) q[2];
sx q[2];
rz(-1.1954824) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6717588) q[1];
sx q[1];
rz(-2.2398754) q[1];
sx q[1];
rz(-1.2059187) q[1];
x q[2];
rz(2.8242802) q[3];
sx q[3];
rz(-2.2307591) q[3];
sx q[3];
rz(2.8803006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2842497) q[2];
sx q[2];
rz(-2.0390859) q[2];
sx q[2];
rz(2.9798689) q[2];
rz(0.11180793) q[3];
sx q[3];
rz(-2.4155152) q[3];
sx q[3];
rz(-2.4017754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1424471) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(2.4430742) q[0];
rz(-1.0922208) q[1];
sx q[1];
rz(-2.3869546) q[1];
sx q[1];
rz(3.0879367) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221401) q[0];
sx q[0];
rz(-0.86556731) q[0];
sx q[0];
rz(-2.7979275) q[0];
rz(-pi) q[1];
rz(-1.1321819) q[2];
sx q[2];
rz(-1.918186) q[2];
sx q[2];
rz(-1.6288822) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3813579) q[1];
sx q[1];
rz(-1.6242199) q[1];
sx q[1];
rz(-2.670856) q[1];
rz(0.14999203) q[3];
sx q[3];
rz(-2.2120471) q[3];
sx q[3];
rz(-0.36442001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1412389) q[2];
sx q[2];
rz(-0.69811368) q[2];
sx q[2];
rz(-0.21256438) q[2];
rz(0.32771787) q[3];
sx q[3];
rz(-0.98723427) q[3];
sx q[3];
rz(1.7879965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.245529) q[0];
sx q[0];
rz(-2.0241757) q[0];
sx q[0];
rz(0.32972202) q[0];
rz(-0.7715191) q[1];
sx q[1];
rz(-0.78116575) q[1];
sx q[1];
rz(-0.82569295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77149174) q[0];
sx q[0];
rz(-2.2326062) q[0];
sx q[0];
rz(-1.6812117) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0631376) q[2];
sx q[2];
rz(-2.1657073) q[2];
sx q[2];
rz(-1.1537781) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6439425) q[1];
sx q[1];
rz(-1.6794599) q[1];
sx q[1];
rz(-0.76246271) q[1];
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
rz(-2.9416435) q[2];
sx q[2];
rz(-1.4120833) q[2];
sx q[2];
rz(-1.4818209) q[2];
rz(3.0485349) q[3];
sx q[3];
rz(-1.8424282) q[3];
sx q[3];
rz(-0.53946462) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2411497) q[0];
sx q[0];
rz(-0.41541442) q[0];
sx q[0];
rz(0.39733091) q[0];
rz(2.1516402) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(2.1053402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5118714) q[0];
sx q[0];
rz(-2.7366834) q[0];
sx q[0];
rz(3.1121503) q[0];
rz(-pi) q[1];
rz(1.7099047) q[2];
sx q[2];
rz(-1.9204307) q[2];
sx q[2];
rz(2.1099427) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.73309988) q[1];
sx q[1];
rz(-1.6599791) q[1];
sx q[1];
rz(2.9217138) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.051722) q[3];
sx q[3];
rz(-1.241893) q[3];
sx q[3];
rz(-0.84612209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1194666) q[2];
sx q[2];
rz(-0.72212044) q[2];
sx q[2];
rz(-1.8252581) q[2];
rz(2.8890166) q[3];
sx q[3];
rz(-1.3467237) q[3];
sx q[3];
rz(2.778497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13755688) q[0];
sx q[0];
rz(-0.79020774) q[0];
sx q[0];
rz(0.23605119) q[0];
rz(0.099418489) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(1.258446) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5262622) q[0];
sx q[0];
rz(-2.1293422) q[0];
sx q[0];
rz(1.6625151) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0706606) q[2];
sx q[2];
rz(-2.1720083) q[2];
sx q[2];
rz(-0.76155969) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.395605) q[1];
sx q[1];
rz(-1.4561936) q[1];
sx q[1];
rz(-0.68273165) q[1];
rz(-pi) q[2];
rz(0.80983617) q[3];
sx q[3];
rz(-2.5812529) q[3];
sx q[3];
rz(-0.59691012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1666169) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(1.0413337) q[2];
rz(0.58639041) q[3];
sx q[3];
rz(-2.6643463) q[3];
sx q[3];
rz(-0.81048107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36622421) q[0];
sx q[0];
rz(-2.0851705) q[0];
sx q[0];
rz(2.0023517) q[0];
rz(-0.99209039) q[1];
sx q[1];
rz(-1.5403668) q[1];
sx q[1];
rz(-1.5425727) q[1];
rz(0.61983776) q[2];
sx q[2];
rz(-0.38417338) q[2];
sx q[2];
rz(-1.2585121) q[2];
rz(2.5358148) q[3];
sx q[3];
rz(-2.45469) q[3];
sx q[3];
rz(-1.5495095) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
