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
rz(0.59745204) q[1];
sx q[1];
rz(-2.6260881) q[1];
sx q[1];
rz(-2.1093624) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55695483) q[0];
sx q[0];
rz(-0.9171066) q[0];
sx q[0];
rz(-1.6980706) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0529673) q[2];
sx q[2];
rz(-1.5368688) q[2];
sx q[2];
rz(0.22824333) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5121442) q[1];
sx q[1];
rz(-1.8820861) q[1];
sx q[1];
rz(-1.6988847) q[1];
rz(-pi) q[2];
rz(0.80268245) q[3];
sx q[3];
rz(-1.2000053) q[3];
sx q[3];
rz(0.30457531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7008179) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(-3.0536998) q[2];
rz(-1.6739738) q[3];
sx q[3];
rz(-1.2940977) q[3];
sx q[3];
rz(-0.99883336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1217594) q[0];
sx q[0];
rz(-1.3061981) q[0];
sx q[0];
rz(1.649296) q[0];
rz(2.3536033) q[1];
sx q[1];
rz(-1.9580656) q[1];
sx q[1];
rz(-0.96510395) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5481297) q[0];
sx q[0];
rz(-2.4060747) q[0];
sx q[0];
rz(-3.0043362) q[0];
x q[1];
rz(2.4048058) q[2];
sx q[2];
rz(-1.4649142) q[2];
sx q[2];
rz(-0.69241619) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3841159) q[1];
sx q[1];
rz(-0.49864492) q[1];
sx q[1];
rz(-1.9668786) q[1];
rz(2.6718465) q[3];
sx q[3];
rz(-1.489991) q[3];
sx q[3];
rz(1.684379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.013177055) q[2];
sx q[2];
rz(-0.33856496) q[2];
sx q[2];
rz(1.817912) q[2];
rz(2.2549818) q[3];
sx q[3];
rz(-1.1611791) q[3];
sx q[3];
rz(2.9429341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24432467) q[0];
sx q[0];
rz(-2.1710945) q[0];
sx q[0];
rz(-2.4460728) q[0];
rz(-2.3059402) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(-0.64170352) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013896522) q[0];
sx q[0];
rz(-0.75152107) q[0];
sx q[0];
rz(1.3579382) q[0];
rz(-pi) q[1];
rz(-2.4303248) q[2];
sx q[2];
rz(-1.4431653) q[2];
sx q[2];
rz(0.39943275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56395951) q[1];
sx q[1];
rz(-1.0924939) q[1];
sx q[1];
rz(0.35525124) q[1];
x q[2];
rz(0.66777457) q[3];
sx q[3];
rz(-2.3333356) q[3];
sx q[3];
rz(-0.86195743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16331638) q[2];
sx q[2];
rz(-1.1288319) q[2];
sx q[2];
rz(1.3345435) q[2];
rz(1.4000019) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(-2.002772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1432994) q[0];
sx q[0];
rz(-0.87109733) q[0];
sx q[0];
rz(0.77348462) q[0];
rz(-3.0553014) q[1];
sx q[1];
rz(-0.50420612) q[1];
sx q[1];
rz(-0.48826826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937259) q[0];
sx q[0];
rz(-2.6742088) q[0];
sx q[0];
rz(0.15166883) q[0];
x q[1];
rz(-2.1147229) q[2];
sx q[2];
rz(-2.1928117) q[2];
sx q[2];
rz(2.029947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77535532) q[1];
sx q[1];
rz(-1.7960494) q[1];
sx q[1];
rz(2.3083592) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9776439) q[3];
sx q[3];
rz(-1.196612) q[3];
sx q[3];
rz(2.0896623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46828541) q[2];
sx q[2];
rz(-1.1136592) q[2];
sx q[2];
rz(2.9735273) q[2];
rz(-0.42260653) q[3];
sx q[3];
rz(-2.1674619) q[3];
sx q[3];
rz(2.9882123) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3733658) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(-1.9510829) q[0];
rz(2.6137784) q[1];
sx q[1];
rz(-2.0595136) q[1];
sx q[1];
rz(-3.1324918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2784257) q[0];
sx q[0];
rz(-2.3519807) q[0];
sx q[0];
rz(-2.4900683) q[0];
rz(-pi) q[1];
x q[1];
rz(0.060623364) q[2];
sx q[2];
rz(-1.6317118) q[2];
sx q[2];
rz(2.2140394) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2796265) q[1];
sx q[1];
rz(-1.1512966) q[1];
sx q[1];
rz(-0.90200938) q[1];
rz(-pi) q[2];
rz(2.1808257) q[3];
sx q[3];
rz(-2.4159263) q[3];
sx q[3];
rz(0.31455597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3147543) q[2];
sx q[2];
rz(-2.645731) q[2];
sx q[2];
rz(0.66519386) q[2];
rz(1.0264171) q[3];
sx q[3];
rz(-2.0668991) q[3];
sx q[3];
rz(2.3608666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5196395) q[0];
sx q[0];
rz(-0.53549796) q[0];
sx q[0];
rz(0.39805472) q[0];
rz(1.4705426) q[1];
sx q[1];
rz(-0.87659756) q[1];
sx q[1];
rz(-0.7801396) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0288643) q[0];
sx q[0];
rz(-2.5327589) q[0];
sx q[0];
rz(0.91797773) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69177912) q[2];
sx q[2];
rz(-2.8304812) q[2];
sx q[2];
rz(1.2300613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2237423) q[1];
sx q[1];
rz(-0.7484352) q[1];
sx q[1];
rz(-2.7176932) q[1];
rz(-1.1885181) q[3];
sx q[3];
rz(-2.4197289) q[3];
sx q[3];
rz(-0.75306276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.857343) q[2];
sx q[2];
rz(-1.1025068) q[2];
sx q[2];
rz(2.9798689) q[2];
rz(-3.0297847) q[3];
sx q[3];
rz(-2.4155152) q[3];
sx q[3];
rz(0.73981729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-2.9991456) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(-0.69851843) q[0];
rz(1.0922208) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(-0.05365595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7258436) q[0];
sx q[0];
rz(-0.77139716) q[0];
sx q[0];
rz(1.1939069) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0094107) q[2];
sx q[2];
rz(-1.918186) q[2];
sx q[2];
rz(1.6288822) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.084890993) q[1];
sx q[1];
rz(-2.6680601) q[1];
sx q[1];
rz(-3.0242306) q[1];
x q[2];
rz(0.14999203) q[3];
sx q[3];
rz(-0.92954554) q[3];
sx q[3];
rz(0.36442001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1412389) q[2];
sx q[2];
rz(-2.443479) q[2];
sx q[2];
rz(2.9290283) q[2];
rz(2.8138748) q[3];
sx q[3];
rz(-0.98723427) q[3];
sx q[3];
rz(-1.7879965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89606365) q[0];
sx q[0];
rz(-1.117417) q[0];
sx q[0];
rz(0.32972202) q[0];
rz(-0.7715191) q[1];
sx q[1];
rz(-0.78116575) q[1];
sx q[1];
rz(-0.82569295) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2742608) q[0];
sx q[0];
rz(-1.6578339) q[0];
sx q[0];
rz(0.66477338) q[0];
x q[1];
rz(0.078455047) q[2];
sx q[2];
rz(-2.1657073) q[2];
sx q[2];
rz(-1.9878146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18620423) q[1];
sx q[1];
rz(-2.3729747) q[1];
sx q[1];
rz(0.15665084) q[1];
rz(-pi) q[2];
rz(-1.3729893) q[3];
sx q[3];
rz(-1.8231492) q[3];
sx q[3];
rz(-0.41393241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9416435) q[2];
sx q[2];
rz(-1.7295094) q[2];
sx q[2];
rz(-1.6597718) q[2];
rz(-0.09305772) q[3];
sx q[3];
rz(-1.8424282) q[3];
sx q[3];
rz(2.602128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.2411497) q[0];
sx q[0];
rz(-2.7261782) q[0];
sx q[0];
rz(0.39733091) q[0];
rz(0.98995248) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(-2.1053402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5118714) q[0];
sx q[0];
rz(-2.7366834) q[0];
sx q[0];
rz(-3.1121503) q[0];
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
rz(-pi) q[0];
x q[0];
rz(-0.81779382) q[1];
sx q[1];
rz(-1.3518055) q[1];
sx q[1];
rz(1.6621672) q[1];
rz(-2.0898706) q[3];
sx q[3];
rz(-1.241893) q[3];
sx q[3];
rz(-2.2954706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.022126023) q[2];
sx q[2];
rz(-2.4194722) q[2];
sx q[2];
rz(-1.3163346) q[2];
rz(0.25257603) q[3];
sx q[3];
rz(-1.3467237) q[3];
sx q[3];
rz(0.36309567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0040358) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(0.23605119) q[0];
rz(-0.099418489) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(-1.258446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61533046) q[0];
sx q[0];
rz(-1.0122504) q[0];
sx q[0];
rz(-1.4790776) q[0];
rz(-pi) q[1];
rz(-2.5314517) q[2];
sx q[2];
rz(-2.379866) q[2];
sx q[2];
rz(1.6118647) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.26791456) q[1];
sx q[1];
rz(-0.89338747) q[1];
sx q[1];
rz(-1.4235086) q[1];
x q[2];
rz(1.997273) q[3];
sx q[3];
rz(-1.1955373) q[3];
sx q[3];
rz(0.29510185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9749757) q[2];
sx q[2];
rz(-2.6585343) q[2];
sx q[2];
rz(1.0413337) q[2];
rz(0.58639041) q[3];
sx q[3];
rz(-0.47724637) q[3];
sx q[3];
rz(0.81048107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36622421) q[0];
sx q[0];
rz(-1.0564221) q[0];
sx q[0];
rz(-1.139241) q[0];
rz(-2.1495023) q[1];
sx q[1];
rz(-1.6012259) q[1];
sx q[1];
rz(1.59902) q[1];
rz(0.61983776) q[2];
sx q[2];
rz(-0.38417338) q[2];
sx q[2];
rz(-1.2585121) q[2];
rz(-0.59320368) q[3];
sx q[3];
rz(-1.9402258) q[3];
sx q[3];
rz(2.6711834) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
