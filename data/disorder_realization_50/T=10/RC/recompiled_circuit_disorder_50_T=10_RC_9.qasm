OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2621736) q[0];
sx q[0];
rz(-1.7466495) q[0];
sx q[0];
rz(-0.0012794415) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2389195) q[0];
sx q[0];
rz(-0.15107778) q[0];
sx q[0];
rz(-0.33688776) q[0];
x q[1];
rz(0.84471976) q[2];
sx q[2];
rz(-1.8978999) q[2];
sx q[2];
rz(-2.6127882) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98382681) q[1];
sx q[1];
rz(-0.14769927) q[1];
sx q[1];
rz(-1.8450792) q[1];
rz(1.418387) q[3];
sx q[3];
rz(-1.8525575) q[3];
sx q[3];
rz(1.8834653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0960192) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(1.1958896) q[2];
rz(-1.1536417) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(-1.7529863) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7213223) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(-1.8700245) q[0];
rz(-2.0416073) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(1.3756479) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028582024) q[0];
sx q[0];
rz(-0.41886371) q[0];
sx q[0];
rz(2.9257665) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0963418) q[2];
sx q[2];
rz(-1.8088733) q[2];
sx q[2];
rz(-2.0766052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1609636) q[1];
sx q[1];
rz(-1.3914319) q[1];
sx q[1];
rz(2.7705926) q[1];
x q[2];
rz(-2.1249173) q[3];
sx q[3];
rz(-1.2614054) q[3];
sx q[3];
rz(-0.644899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.21330825) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(0.48970547) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(2.9009853) q[0];
rz(-0.34257564) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(-1.2352357) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2429758) q[0];
sx q[0];
rz(-1.9447864) q[0];
sx q[0];
rz(2.5234733) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3182993) q[2];
sx q[2];
rz(-1.9324979) q[2];
sx q[2];
rz(-2.8351438) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.55360824) q[1];
sx q[1];
rz(-1.8318892) q[1];
sx q[1];
rz(-1.8334465) q[1];
x q[2];
rz(2.5123185) q[3];
sx q[3];
rz(-0.73262239) q[3];
sx q[3];
rz(-0.34793684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3774595) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(1.4366478) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075994611) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(-0.80379379) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0474931) q[0];
sx q[0];
rz(-2.1764628) q[0];
sx q[0];
rz(3.0981307) q[0];
rz(-1.7159963) q[2];
sx q[2];
rz(-1.928066) q[2];
sx q[2];
rz(1.4841339) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67112982) q[1];
sx q[1];
rz(-1.5531335) q[1];
sx q[1];
rz(-1.3385685) q[1];
rz(-2.1060137) q[3];
sx q[3];
rz(-1.03973) q[3];
sx q[3];
rz(-0.5862743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(0.072337739) q[2];
rz(2.7667601) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(-1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(-2.6547292) q[0];
rz(0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(-1.240085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4935721) q[0];
sx q[0];
rz(-2.9020712) q[0];
sx q[0];
rz(-0.055672107) q[0];
rz(0.83112049) q[2];
sx q[2];
rz(-0.34905012) q[2];
sx q[2];
rz(-2.8380053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8023194) q[1];
sx q[1];
rz(-0.60802751) q[1];
sx q[1];
rz(0.1165216) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48360444) q[3];
sx q[3];
rz(-2.1198453) q[3];
sx q[3];
rz(-0.99398127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.038625) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(-2.4385578) q[2];
rz(1.4098343) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(-2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(2.1512206) q[0];
rz(3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.4809158) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71125644) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(-1.3769763) q[0];
x q[1];
rz(2.0159857) q[2];
sx q[2];
rz(-1.5700969) q[2];
sx q[2];
rz(0.47847754) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33038352) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(-1.1793544) q[1];
x q[2];
rz(-3.092993) q[3];
sx q[3];
rz(-2.2489293) q[3];
sx q[3];
rz(-2.5130659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(2.5793502) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(0.26091584) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19787191) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(-1.4690171) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(1.3247103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0425456) q[0];
sx q[0];
rz(-1.5934266) q[0];
sx q[0];
rz(3.0142473) q[0];
x q[1];
rz(-2.1340738) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(-1.8958467) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.087254698) q[1];
sx q[1];
rz(-1.195765) q[1];
sx q[1];
rz(-2.7021728) q[1];
x q[2];
rz(0.38447325) q[3];
sx q[3];
rz(-2.1914346) q[3];
sx q[3];
rz(1.4409161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5801195) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(-0.44357792) q[2];
rz(-0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401944) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(-0.46052128) q[0];
rz(3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(-1.8519648) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3279009) q[0];
sx q[0];
rz(-1.1192338) q[0];
sx q[0];
rz(-0.62664647) q[0];
rz(-0.71596594) q[2];
sx q[2];
rz(-0.95574524) q[2];
sx q[2];
rz(-1.3032608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8045115) q[1];
sx q[1];
rz(-1.6816499) q[1];
sx q[1];
rz(-1.3827419) q[1];
rz(-0.66289785) q[3];
sx q[3];
rz(-1.7015966) q[3];
sx q[3];
rz(2.7152284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9926247) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643395) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(-2.8073231) q[0];
rz(-1.224068) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(0.28265488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6219014) q[0];
sx q[0];
rz(-2.0265409) q[0];
sx q[0];
rz(-0.68463188) q[0];
rz(-pi) q[1];
rz(-2.5731509) q[2];
sx q[2];
rz(-2.754515) q[2];
sx q[2];
rz(-2.8708411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54284755) q[1];
sx q[1];
rz(-1.3951021) q[1];
sx q[1];
rz(2.7302242) q[1];
x q[2];
rz(3.024179) q[3];
sx q[3];
rz(-2.252929) q[3];
sx q[3];
rz(0.62473245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21215542) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(2.4003417) q[2];
rz(2.6397928) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0989477) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(-2.6328971) q[0];
rz(-0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(-0.68181109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1334575) q[0];
sx q[0];
rz(-1.4032161) q[0];
sx q[0];
rz(-1.0159147) q[0];
x q[1];
rz(0.85068662) q[2];
sx q[2];
rz(-2.2518573) q[2];
sx q[2];
rz(2.1249287) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21916418) q[1];
sx q[1];
rz(-1.0370967) q[1];
sx q[1];
rz(1.0132755) q[1];
rz(-pi) q[2];
rz(1.0107332) q[3];
sx q[3];
rz(-2.0252953) q[3];
sx q[3];
rz(1.0487995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8455785) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(2.8005023) q[2];
rz(1.0836481) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(0.60733168) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(1.0495937) q[2];
sx q[2];
rz(-0.97432077) q[2];
sx q[2];
rz(-1.5885098) q[2];
rz(-1.3689465) q[3];
sx q[3];
rz(-1.9484083) q[3];
sx q[3];
rz(1.2231135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
