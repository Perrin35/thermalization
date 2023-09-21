OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(-1.3949431) q[0];
sx q[0];
rz(-3.1403132) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(4.2445634) q[1];
sx q[1];
rz(7.0581262) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66520663) q[0];
sx q[0];
rz(-1.6205661) q[0];
sx q[0];
rz(-2.9988891) q[0];
rz(-pi) q[1];
rz(-2.7156419) q[2];
sx q[2];
rz(-0.89086878) q[2];
sx q[2];
rz(-1.8217063) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1577658) q[1];
sx q[1];
rz(-0.14769927) q[1];
sx q[1];
rz(-1.2965135) q[1];
rz(-pi) q[2];
rz(-1.7232056) q[3];
sx q[3];
rz(-1.2890352) q[3];
sx q[3];
rz(1.2581274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0960192) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(1.9457031) q[2];
rz(1.1536417) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(-1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.4202704) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(1.0999854) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(-1.3756479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028582024) q[0];
sx q[0];
rz(-0.41886371) q[0];
sx q[0];
rz(-2.9257665) q[0];
rz(-1.8091082) q[2];
sx q[2];
rz(-1.5268227) q[2];
sx q[2];
rz(2.6251052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.980629) q[1];
sx q[1];
rz(-1.7501608) q[1];
sx q[1];
rz(0.37100002) q[1];
rz(-2.7820285) q[3];
sx q[3];
rz(-2.095795) q[3];
sx q[3];
rz(-1.1121225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21330825) q[2];
sx q[2];
rz(-1.7958612) q[2];
sx q[2];
rz(-0.48970547) q[2];
rz(2.0080163) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60004822) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(2.9009853) q[0];
rz(0.34257564) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(-1.906357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92634) q[0];
sx q[0];
rz(-2.140576) q[0];
sx q[0];
rz(-2.0195872) q[0];
x q[1];
rz(-1.8232934) q[2];
sx q[2];
rz(-1.9324979) q[2];
sx q[2];
rz(0.30644882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5879844) q[1];
sx q[1];
rz(-1.3097035) q[1];
sx q[1];
rz(1.3081461) q[1];
x q[2];
rz(-2.5123185) q[3];
sx q[3];
rz(-0.73262239) q[3];
sx q[3];
rz(0.34793684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76413313) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(-1.4366478) q[2];
rz(-1.7403587) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(-2.5333372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(-0.11985699) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940995) q[0];
sx q[0];
rz(-2.1764628) q[0];
sx q[0];
rz(3.0981307) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7718133) q[2];
sx q[2];
rz(-2.7571207) q[2];
sx q[2];
rz(-2.0535) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67112982) q[1];
sx q[1];
rz(-1.5884591) q[1];
sx q[1];
rz(1.3385685) q[1];
rz(-pi) q[2];
rz(1.0355789) q[3];
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
rz(-3.0692549) q[2];
rz(2.7667601) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(-1.9434631) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(-0.48686349) q[0];
rz(2.411719) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.240085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2729028) q[0];
sx q[0];
rz(-1.5839974) q[0];
sx q[0];
rz(2.9024283) q[0];
x q[1];
rz(1.8334332) q[2];
sx q[2];
rz(-1.3381759) q[2];
sx q[2];
rz(0.55839415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1357437) q[1];
sx q[1];
rz(-1.6372576) q[1];
sx q[1];
rz(0.60484109) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1754856) q[3];
sx q[3];
rz(-1.1629512) q[3];
sx q[3];
rz(2.2972578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.038625) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(1.7317584) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4617417) q[0];
sx q[0];
rz(-1.6438419) q[0];
sx q[0];
rz(1.9528271) q[0];
x q[1];
rz(-1.5724206) q[2];
sx q[2];
rz(-2.6964028) q[2];
sx q[2];
rz(-1.0908529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78696886) q[1];
sx q[1];
rz(-2.0957392) q[1];
sx q[1];
rz(2.900219) q[1];
rz(1.6310286) q[3];
sx q[3];
rz(-2.4619953) q[3];
sx q[3];
rz(2.4356902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(-2.5793502) q[2];
rz(-1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(-0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19787191) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(1.4690171) q[0];
rz(2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(1.3247103) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64667386) q[0];
sx q[0];
rz(-3.0122628) q[0];
sx q[0];
rz(-2.9652251) q[0];
rz(2.1340738) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(1.8958467) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.054338) q[1];
sx q[1];
rz(-1.195765) q[1];
sx q[1];
rz(0.43941984) q[1];
rz(2.7571194) q[3];
sx q[3];
rz(-2.1914346) q[3];
sx q[3];
rz(1.7006765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(-0.44357792) q[2];
rz(-2.1929072) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.401944) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(-0.46052128) q[0];
rz(3.0415688) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(-1.2896279) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063110654) q[0];
sx q[0];
rz(-1.0149628) q[0];
sx q[0];
rz(1.0312992) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82215674) q[2];
sx q[2];
rz(-0.90688721) q[2];
sx q[2];
rz(0.31809959) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8045115) q[1];
sx q[1];
rz(-1.4599428) q[1];
sx q[1];
rz(1.7588508) q[1];
x q[2];
rz(2.9309978) q[3];
sx q[3];
rz(-2.4678293) q[3];
sx q[3];
rz(1.3099561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.148968) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(1.0827433) q[2];
rz(-3.0454214) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7643395) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(-0.33426958) q[0];
rz(1.224068) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(0.28265488) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7455505) q[0];
sx q[0];
rz(-0.96691416) q[0];
sx q[0];
rz(2.1349483) q[0];
x q[1];
rz(-1.3547782) q[2];
sx q[2];
rz(-1.8945433) q[2];
sx q[2];
rz(2.266989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0375367) q[1];
sx q[1];
rz(-1.9754585) q[1];
sx q[1];
rz(-1.7621103) q[1];
rz(-pi) q[2];
rz(1.4275527) q[3];
sx q[3];
rz(-0.69056615) q[3];
sx q[3];
rz(-2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9294372) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0989477) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(-3.0264061) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(-0.68181109) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9671229) q[0];
sx q[0];
rz(-2.5645064) q[0];
sx q[0];
rz(-1.2601) q[0];
rz(-pi) q[1];
x q[1];
rz(2.290906) q[2];
sx q[2];
rz(-2.2518573) q[2];
sx q[2];
rz(-2.1249287) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21916418) q[1];
sx q[1];
rz(-2.104496) q[1];
sx q[1];
rz(-2.1283172) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0107332) q[3];
sx q[3];
rz(-2.0252953) q[3];
sx q[3];
rz(-2.0927932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8455785) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(1.0836481) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(-2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-0.66424673) q[2];
sx q[2];
rz(-1.9953809) q[2];
sx q[2];
rz(-2.8473163) q[2];
rz(-1.7726462) q[3];
sx q[3];
rz(-1.1931843) q[3];
sx q[3];
rz(-1.9184792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];