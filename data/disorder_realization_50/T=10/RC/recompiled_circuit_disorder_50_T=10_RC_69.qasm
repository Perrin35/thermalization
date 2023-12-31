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
rz(3.1403132) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(2.3666518) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89844184) q[0];
sx q[0];
rz(-1.4282707) q[0];
sx q[0];
rz(-1.5205163) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7156419) q[2];
sx q[2];
rz(-2.2507239) q[2];
sx q[2];
rz(-1.8217063) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.70667627) q[1];
sx q[1];
rz(-1.7129363) q[1];
sx q[1];
rz(0.040277004) q[1];
rz(-pi) q[2];
rz(-2.6585456) q[3];
sx q[3];
rz(-0.31937283) q[3];
sx q[3];
rz(-1.3787624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(-1.9457031) q[2];
rz(1.1536417) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.4202704) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(-1.2715682) q[0];
rz(1.0999854) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(1.3756479) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1130106) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(2.9257665) q[0];
rz(-pi) q[1];
x q[1];
rz(0.045250821) q[2];
sx q[2];
rz(-1.8088733) q[2];
sx q[2];
rz(2.0766052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1609636) q[1];
sx q[1];
rz(-1.3914319) q[1];
sx q[1];
rz(-2.7705926) q[1];
x q[2];
rz(0.35956412) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(-2.0294702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-0.48970547) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(-2.074923) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(0.24060732) q[0];
rz(-0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.906357) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19757195) q[0];
sx q[0];
rz(-2.4320142) q[0];
sx q[0];
rz(-0.59528415) q[0];
rz(-1.3182993) q[2];
sx q[2];
rz(-1.9324979) q[2];
sx q[2];
rz(-0.30644882) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0551128) q[1];
sx q[1];
rz(-1.8243454) q[1];
sx q[1];
rz(2.8716645) q[1];
x q[2];
rz(2.0577621) q[3];
sx q[3];
rz(-0.999513) q[3];
sx q[3];
rz(-2.7146102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76413313) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(1.7403587) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(2.5333372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.065598) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(-2.3377989) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(-3.0217357) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940995) q[0];
sx q[0];
rz(-0.96512981) q[0];
sx q[0];
rz(-0.043461965) q[0];
rz(-1.7159963) q[2];
sx q[2];
rz(-1.928066) q[2];
sx q[2];
rz(1.4841339) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.67112982) q[1];
sx q[1];
rz(-1.5531335) q[1];
sx q[1];
rz(-1.3385685) q[1];
rz(-pi) q[2];
rz(2.4265392) q[3];
sx q[3];
rz(-0.73521571) q[3];
sx q[3];
rz(-0.27763593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(-0.48686349) q[0];
rz(0.72987366) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(1.240085) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4935721) q[0];
sx q[0];
rz(-0.23952142) q[0];
sx q[0];
rz(0.055672107) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.901022) q[2];
sx q[2];
rz(-1.3153937) q[2];
sx q[2];
rz(-2.0672928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33927321) q[1];
sx q[1];
rz(-2.5335651) q[1];
sx q[1];
rz(-0.1165216) q[1];
rz(-0.48360444) q[3];
sx q[3];
rz(-1.0217474) q[3];
sx q[3];
rz(2.1476114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(1.4098343) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-2.1773961) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(-2.1512206) q[0];
rz(-0.05274996) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(1.6606768) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71125644) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(1.7646164) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5691721) q[2];
sx q[2];
rz(-2.6964028) q[2];
sx q[2];
rz(-1.0908529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2350125) q[1];
sx q[1];
rz(-1.7791516) q[1];
sx q[1];
rz(2.1085897) q[1];
rz(-0.048599676) q[3];
sx q[3];
rz(-2.2489293) q[3];
sx q[3];
rz(2.5130659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1331553) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(2.0810614) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9437207) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.4690171) q[0];
rz(2.127227) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(1.8168824) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6727407) q[0];
sx q[0];
rz(-1.6981089) q[0];
sx q[0];
rz(-1.5479814) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0075188) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(1.8958467) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8285671) q[1];
sx q[1];
rz(-1.1638068) q[1];
sx q[1];
rz(1.160497) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2277101) q[3];
sx q[3];
rz(-1.2607288) q[3];
sx q[3];
rz(-3.0403746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(2.6980147) q[2];
rz(-0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(-2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396486) q[0];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3562718) q[0];
sx q[0];
rz(-2.3873781) q[0];
sx q[0];
rz(0.69099364) q[0];
rz(0.71596594) q[2];
sx q[2];
rz(-0.95574524) q[2];
sx q[2];
rz(-1.8383319) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.928927) q[1];
sx q[1];
rz(-1.7576828) q[1];
sx q[1];
rz(0.11282632) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66289785) q[3];
sx q[3];
rz(-1.7015966) q[3];
sx q[3];
rz(-0.42636426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.9926247) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(1.0827433) q[2];
rz(0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37725317) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(-0.33426958) q[0];
rz(1.224068) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(-0.28265488) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7455505) q[0];
sx q[0];
rz(-0.96691416) q[0];
sx q[0];
rz(-2.1349483) q[0];
x q[1];
rz(-0.33090584) q[2];
sx q[2];
rz(-1.7754284) q[2];
sx q[2];
rz(-2.3757039) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4947195) q[1];
sx q[1];
rz(-0.44534007) q[1];
sx q[1];
rz(0.41782197) q[1];
rz(0.11741365) q[3];
sx q[3];
rz(-2.252929) q[3];
sx q[3];
rz(2.5168602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(-2.6397928) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989477) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(-0.50869554) q[0];
rz(-3.0264061) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(0.68181109) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17446974) q[0];
sx q[0];
rz(-2.5645064) q[0];
sx q[0];
rz(-1.8814927) q[0];
rz(-pi) q[1];
rz(-2.3186458) q[2];
sx q[2];
rz(-1.032885) q[2];
sx q[2];
rz(-0.049494628) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9224285) q[1];
sx q[1];
rz(-1.0370967) q[1];
sx q[1];
rz(-1.0132755) q[1];
x q[2];
rz(-2.1308594) q[3];
sx q[3];
rz(-1.1162973) q[3];
sx q[3];
rz(2.0927932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8455785) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(2.0579445) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-0.60733168) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(-0.63275679) q[2];
sx q[2];
rz(-2.3709595) q[2];
sx q[2];
rz(-0.79216935) q[2];
rz(0.38469436) q[3];
sx q[3];
rz(-1.3833429) q[3];
sx q[3];
rz(-0.27237567) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
