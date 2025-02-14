OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.628217) q[0];
sx q[0];
rz(-1.4613419) q[0];
sx q[0];
rz(-2.2267377) q[0];
rz(2.272361) q[1];
sx q[1];
rz(-0.72328049) q[1];
sx q[1];
rz(1.3528612) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992743) q[0];
sx q[0];
rz(-0.87142105) q[0];
sx q[0];
rz(-0.69518881) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3981744) q[2];
sx q[2];
rz(-1.8525436) q[2];
sx q[2];
rz(-2.8670058) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9661508) q[1];
sx q[1];
rz(-1.3805318) q[1];
sx q[1];
rz(-0.89118608) q[1];
rz(-0.21463359) q[3];
sx q[3];
rz(-2.1691536) q[3];
sx q[3];
rz(-2.2325688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.73504084) q[2];
sx q[2];
rz(-2.328379) q[2];
sx q[2];
rz(2.1347894) q[2];
rz(-1.6793647) q[3];
sx q[3];
rz(-1.6961325) q[3];
sx q[3];
rz(-1.538895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2891069) q[0];
sx q[0];
rz(-0.91133457) q[0];
sx q[0];
rz(3.015633) q[0];
rz(1.1360629) q[1];
sx q[1];
rz(-1.2964396) q[1];
sx q[1];
rz(-1.0596641) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019106796) q[0];
sx q[0];
rz(-3.1395432) q[0];
sx q[0];
rz(-0.37197427) q[0];
rz(-pi) q[1];
rz(-2.4074209) q[2];
sx q[2];
rz(-2.3245272) q[2];
sx q[2];
rz(-1.2881607) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8318924) q[1];
sx q[1];
rz(-0.87617481) q[1];
sx q[1];
rz(-2.5599856) q[1];
x q[2];
rz(2.7823388) q[3];
sx q[3];
rz(-1.6174966) q[3];
sx q[3];
rz(-0.31622313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5560567) q[2];
sx q[2];
rz(-1.5289331) q[2];
sx q[2];
rz(0.6297099) q[2];
rz(1.2398047) q[3];
sx q[3];
rz(-2.3972378) q[3];
sx q[3];
rz(-1.121421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12104163) q[0];
sx q[0];
rz(-2.910055) q[0];
sx q[0];
rz(-1.2303906) q[0];
rz(-2.4379099) q[1];
sx q[1];
rz(-0.92862248) q[1];
sx q[1];
rz(1.5812662) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0315097) q[0];
sx q[0];
rz(-0.39387273) q[0];
sx q[0];
rz(-1.5316233) q[0];
rz(-0.21321984) q[2];
sx q[2];
rz(-0.74951321) q[2];
sx q[2];
rz(2.9507278) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3250503) q[1];
sx q[1];
rz(-2.5245164) q[1];
sx q[1];
rz(-2.3338334) q[1];
rz(-pi) q[2];
rz(-0.78572598) q[3];
sx q[3];
rz(-1.2961939) q[3];
sx q[3];
rz(0.532632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8065717) q[2];
sx q[2];
rz(-1.0806012) q[2];
sx q[2];
rz(-3.0873155) q[2];
rz(-2.055376) q[3];
sx q[3];
rz(-0.79038668) q[3];
sx q[3];
rz(2.6495892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7739968) q[0];
sx q[0];
rz(-3.1258899) q[0];
sx q[0];
rz(0.98454654) q[0];
rz(1.7355851) q[1];
sx q[1];
rz(-1.5407298) q[1];
sx q[1];
rz(2.9071009) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67574745) q[0];
sx q[0];
rz(-2.391531) q[0];
sx q[0];
rz(-1.1613599) q[0];
rz(-pi) q[1];
rz(-0.86799829) q[2];
sx q[2];
rz(-2.7722087) q[2];
sx q[2];
rz(-0.59472769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76012661) q[1];
sx q[1];
rz(-2.1153808) q[1];
sx q[1];
rz(-1.8167956) q[1];
rz(0.99458875) q[3];
sx q[3];
rz(-1.595701) q[3];
sx q[3];
rz(-1.8751609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24388127) q[2];
sx q[2];
rz(-0.66456866) q[2];
sx q[2];
rz(-0.74771869) q[2];
rz(-2.054935) q[3];
sx q[3];
rz(-2.1472411) q[3];
sx q[3];
rz(-2.7896816) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9413167) q[0];
sx q[0];
rz(-0.87012297) q[0];
sx q[0];
rz(-1.5927607) q[0];
rz(-3.0039655) q[1];
sx q[1];
rz(-0.96447861) q[1];
sx q[1];
rz(0.67759883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0726822) q[0];
sx q[0];
rz(-0.015791647) q[0];
sx q[0];
rz(2.6591405) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65998544) q[2];
sx q[2];
rz(-0.2435682) q[2];
sx q[2];
rz(3.0483831) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92050925) q[1];
sx q[1];
rz(-1.5649321) q[1];
sx q[1];
rz(1.7076034) q[1];
rz(1.8137535) q[3];
sx q[3];
rz(-0.84872171) q[3];
sx q[3];
rz(0.18096443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.64348334) q[2];
sx q[2];
rz(-0.52916873) q[2];
sx q[2];
rz(0.37437487) q[2];
rz(2.3937285) q[3];
sx q[3];
rz(-1.6467843) q[3];
sx q[3];
rz(1.098746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.13136524) q[0];
sx q[0];
rz(-1.3656728) q[0];
sx q[0];
rz(-2.9170872) q[0];
rz(2.7911216) q[1];
sx q[1];
rz(-1.1843362) q[1];
sx q[1];
rz(-1.1856461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.910584) q[0];
sx q[0];
rz(-1.5155927) q[0];
sx q[0];
rz(1.7189222) q[0];
rz(1.6669438) q[2];
sx q[2];
rz(-2.0324869) q[2];
sx q[2];
rz(-2.8682414) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1846879) q[1];
sx q[1];
rz(-2.3010151) q[1];
sx q[1];
rz(-2.5504677) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4734201) q[3];
sx q[3];
rz(-2.1134317) q[3];
sx q[3];
rz(2.150491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6422358) q[2];
sx q[2];
rz(-0.88358742) q[2];
sx q[2];
rz(2.9767766) q[2];
rz(-2.3757101) q[3];
sx q[3];
rz(-0.82847786) q[3];
sx q[3];
rz(0.59144056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0484475) q[0];
sx q[0];
rz(-3.1146545) q[0];
sx q[0];
rz(-0.41937605) q[0];
rz(1.558149) q[1];
sx q[1];
rz(-0.42834586) q[1];
sx q[1];
rz(1.8289808) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87154138) q[0];
sx q[0];
rz(-1.7333374) q[0];
sx q[0];
rz(0.97515653) q[0];
rz(-pi) q[1];
rz(0.75757005) q[2];
sx q[2];
rz(-0.84572863) q[2];
sx q[2];
rz(0.49882327) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98931354) q[1];
sx q[1];
rz(-1.3817467) q[1];
sx q[1];
rz(1.3192654) q[1];
rz(-2.6978771) q[3];
sx q[3];
rz(-0.29510083) q[3];
sx q[3];
rz(-2.091311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7340362) q[2];
sx q[2];
rz(-0.67136705) q[2];
sx q[2];
rz(-2.7835795) q[2];
rz(-0.19896209) q[3];
sx q[3];
rz(-0.81148654) q[3];
sx q[3];
rz(0.097804047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99893779) q[0];
sx q[0];
rz(-1.0336579) q[0];
sx q[0];
rz(0.79310966) q[0];
rz(2.2456887) q[1];
sx q[1];
rz(-1.9205807) q[1];
sx q[1];
rz(2.6633247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2127599) q[0];
sx q[0];
rz(-2.6750155) q[0];
sx q[0];
rz(-0.27884941) q[0];
x q[1];
rz(-1.2094648) q[2];
sx q[2];
rz(-1.1602957) q[2];
sx q[2];
rz(0.77958661) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.079887159) q[1];
sx q[1];
rz(-1.8338039) q[1];
sx q[1];
rz(2.0531462) q[1];
rz(-0.10981126) q[3];
sx q[3];
rz(-1.0581019) q[3];
sx q[3];
rz(1.5231665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6200977) q[2];
sx q[2];
rz(-1.564881) q[2];
sx q[2];
rz(-0.0077956789) q[2];
rz(-1.1826285) q[3];
sx q[3];
rz(-1.3820442) q[3];
sx q[3];
rz(1.2953322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1496534) q[0];
sx q[0];
rz(-2.6872334) q[0];
sx q[0];
rz(-0.40865189) q[0];
rz(1.0952134) q[1];
sx q[1];
rz(-1.4827671) q[1];
sx q[1];
rz(-0.85831395) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82018554) q[0];
sx q[0];
rz(-1.6371861) q[0];
sx q[0];
rz(1.7184577) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0093391) q[2];
sx q[2];
rz(-2.0200854) q[2];
sx q[2];
rz(1.1960891) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5206362) q[1];
sx q[1];
rz(-2.9163216) q[1];
sx q[1];
rz(2.0849289) q[1];
rz(-pi) q[2];
rz(1.7935221) q[3];
sx q[3];
rz(-2.1568884) q[3];
sx q[3];
rz(3.0243788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3506763) q[2];
sx q[2];
rz(-2.3886267) q[2];
sx q[2];
rz(-2.47993) q[2];
rz(2.3412797) q[3];
sx q[3];
rz(-1.1626264) q[3];
sx q[3];
rz(-0.77320981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8372339) q[0];
sx q[0];
rz(-0.14166129) q[0];
sx q[0];
rz(-0.71453553) q[0];
rz(-0.47572687) q[1];
sx q[1];
rz(-2.5560502) q[1];
sx q[1];
rz(-2.661396) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44549756) q[0];
sx q[0];
rz(-1.0871743) q[0];
sx q[0];
rz(0.66894834) q[0];
rz(-0.33145655) q[2];
sx q[2];
rz(-0.94539795) q[2];
sx q[2];
rz(-1.3783704) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.46327766) q[1];
sx q[1];
rz(-0.92683219) q[1];
sx q[1];
rz(-1.2106845) q[1];
rz(-1.1350182) q[3];
sx q[3];
rz(-1.8284855) q[3];
sx q[3];
rz(2.1439752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6314038) q[2];
sx q[2];
rz(-1.9437342) q[2];
sx q[2];
rz(-2.8655444) q[2];
rz(-2.7035642) q[3];
sx q[3];
rz(-1.5166401) q[3];
sx q[3];
rz(-2.5146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.20090564) q[0];
sx q[0];
rz(-1.154366) q[0];
sx q[0];
rz(0.05703297) q[0];
rz(3.1055462) q[1];
sx q[1];
rz(-1.6135975) q[1];
sx q[1];
rz(1.5560908) q[1];
rz(1.0895928) q[2];
sx q[2];
rz(-2.5960428) q[2];
sx q[2];
rz(-3.0910603) q[2];
rz(1.7071758) q[3];
sx q[3];
rz(-1.3188667) q[3];
sx q[3];
rz(1.8018166) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
