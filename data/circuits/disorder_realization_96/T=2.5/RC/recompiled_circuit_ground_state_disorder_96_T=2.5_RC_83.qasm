OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.51337564) q[0];
sx q[0];
rz(-1.6802508) q[0];
sx q[0];
rz(-0.914855) q[0];
rz(2.272361) q[1];
sx q[1];
rz(2.4183122) q[1];
sx q[1];
rz(8.0719168) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9147089) q[0];
sx q[0];
rz(-1.0585203) q[0];
sx q[0];
rz(-0.73988503) q[0];
rz(-pi) q[1];
rz(1.1961785) q[2];
sx q[2];
rz(-2.2785419) q[2];
sx q[2];
rz(-1.5951235) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1754419) q[1];
sx q[1];
rz(-1.7610608) q[1];
sx q[1];
rz(2.2504066) q[1];
rz(-2.1800024) q[3];
sx q[3];
rz(-1.7477027) q[3];
sx q[3];
rz(0.78395236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4065518) q[2];
sx q[2];
rz(-0.81321365) q[2];
sx q[2];
rz(1.0068033) q[2];
rz(1.6793647) q[3];
sx q[3];
rz(-1.4454602) q[3];
sx q[3];
rz(1.6026976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2891069) q[0];
sx q[0];
rz(-2.2302581) q[0];
sx q[0];
rz(-3.015633) q[0];
rz(2.0055298) q[1];
sx q[1];
rz(-1.845153) q[1];
sx q[1];
rz(-1.0596641) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019106796) q[0];
sx q[0];
rz(-0.0020494941) q[0];
sx q[0];
rz(-0.37197427) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7341717) q[2];
sx q[2];
rz(-0.81706542) q[2];
sx q[2];
rz(1.2881607) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4821533) q[1];
sx q[1];
rz(-2.0065313) q[1];
sx q[1];
rz(2.3547291) q[1];
x q[2];
rz(-0.35925389) q[3];
sx q[3];
rz(-1.5240961) q[3];
sx q[3];
rz(-2.8253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58553592) q[2];
sx q[2];
rz(-1.5289331) q[2];
sx q[2];
rz(-2.5118828) q[2];
rz(1.9017879) q[3];
sx q[3];
rz(-2.3972378) q[3];
sx q[3];
rz(1.121421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(3.020551) q[0];
sx q[0];
rz(-0.23153767) q[0];
sx q[0];
rz(1.9112021) q[0];
rz(-0.70368272) q[1];
sx q[1];
rz(-2.2129702) q[1];
sx q[1];
rz(-1.5603265) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7170555) q[0];
sx q[0];
rz(-1.5858264) q[0];
sx q[0];
rz(-1.9643972) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7652487) q[2];
sx q[2];
rz(-0.84215763) q[2];
sx q[2];
rz(3.0448845) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3250503) q[1];
sx q[1];
rz(-0.61707622) q[1];
sx q[1];
rz(2.3338334) q[1];
x q[2];
rz(0.78572598) q[3];
sx q[3];
rz(-1.8453987) q[3];
sx q[3];
rz(-2.6089607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8065717) q[2];
sx q[2];
rz(-1.0806012) q[2];
sx q[2];
rz(0.054277167) q[2];
rz(1.0862167) q[3];
sx q[3];
rz(-2.351206) q[3];
sx q[3];
rz(-2.6495892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739968) q[0];
sx q[0];
rz(-3.1258899) q[0];
sx q[0];
rz(2.1570461) q[0];
rz(-1.4060075) q[1];
sx q[1];
rz(-1.6008629) q[1];
sx q[1];
rz(0.23449177) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2024883) q[0];
sx q[0];
rz(-1.8456158) q[0];
sx q[0];
rz(-0.86357926) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.245204) q[2];
sx q[2];
rz(-1.8498932) q[2];
sx q[2];
rz(0.14268219) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.76012661) q[1];
sx q[1];
rz(-1.0262118) q[1];
sx q[1];
rz(-1.3247971) q[1];
rz(3.111895) q[3];
sx q[3];
rz(-2.1468024) q[3];
sx q[3];
rz(-2.821049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24388127) q[2];
sx q[2];
rz(-2.477024) q[2];
sx q[2];
rz(-0.74771869) q[2];
rz(-1.0866577) q[3];
sx q[3];
rz(-0.99435157) q[3];
sx q[3];
rz(0.35191107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20027593) q[0];
sx q[0];
rz(-2.2714697) q[0];
sx q[0];
rz(-1.548832) q[0];
rz(-3.0039655) q[1];
sx q[1];
rz(-0.96447861) q[1];
sx q[1];
rz(0.67759883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068910412) q[0];
sx q[0];
rz(-3.125801) q[0];
sx q[0];
rz(-0.48245211) q[0];
rz(0.19385152) q[2];
sx q[2];
rz(-1.7192013) q[2];
sx q[2];
rz(-1.0184763) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4487344) q[1];
sx q[1];
rz(-0.13693196) q[1];
sx q[1];
rz(1.5278234) q[1];
x q[2];
rz(-2.8749488) q[3];
sx q[3];
rz(-2.3867749) q[3];
sx q[3];
rz(-0.17780576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64348334) q[2];
sx q[2];
rz(-0.52916873) q[2];
sx q[2];
rz(-0.37437487) q[2];
rz(-0.74786413) q[3];
sx q[3];
rz(-1.6467843) q[3];
sx q[3];
rz(-2.0428467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0102274) q[0];
sx q[0];
rz(-1.7759198) q[0];
sx q[0];
rz(0.22450547) q[0];
rz(-0.35047105) q[1];
sx q[1];
rz(-1.1843362) q[1];
sx q[1];
rz(-1.1856461) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.910584) q[0];
sx q[0];
rz(-1.5155927) q[0];
sx q[0];
rz(-1.4226705) q[0];
rz(1.4746488) q[2];
sx q[2];
rz(-2.0324869) q[2];
sx q[2];
rz(2.8682414) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19294944) q[1];
sx q[1];
rz(-1.1426326) q[1];
sx q[1];
rz(-2.3938426) q[1];
x q[2];
rz(-0.77187431) q[3];
sx q[3];
rz(-2.3079685) q[3];
sx q[3];
rz(1.1590568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6422358) q[2];
sx q[2];
rz(-2.2580052) q[2];
sx q[2];
rz(-0.16481608) q[2];
rz(-2.3757101) q[3];
sx q[3];
rz(-2.3131148) q[3];
sx q[3];
rz(-0.59144056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0931452) q[0];
sx q[0];
rz(-0.026938139) q[0];
sx q[0];
rz(0.41937605) q[0];
rz(-1.5834437) q[1];
sx q[1];
rz(-2.7132468) q[1];
sx q[1];
rz(-1.8289808) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2079577) q[0];
sx q[0];
rz(-0.61481732) q[0];
sx q[0];
rz(-1.2864248) q[0];
rz(-pi) q[1];
rz(0.75757005) q[2];
sx q[2];
rz(-0.84572863) q[2];
sx q[2];
rz(-2.6427694) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98931354) q[1];
sx q[1];
rz(-1.759846) q[1];
sx q[1];
rz(1.3192654) q[1];
rz(2.873654) q[3];
sx q[3];
rz(-1.4456141) q[3];
sx q[3];
rz(0.94733688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7340362) q[2];
sx q[2];
rz(-0.67136705) q[2];
sx q[2];
rz(2.7835795) q[2];
rz(-2.9426306) q[3];
sx q[3];
rz(-2.3301061) q[3];
sx q[3];
rz(0.097804047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.47826794) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90253622) q[0];
sx q[0];
rz(-1.1235813) q[0];
sx q[0];
rz(-1.7085525) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9321279) q[2];
sx q[2];
rz(-1.981297) q[2];
sx q[2];
rz(2.362006) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9517802) q[1];
sx q[1];
rz(-0.54439044) q[1];
sx q[1];
rz(2.0967099) q[1];
rz(-1.0555154) q[3];
sx q[3];
rz(-1.4751504) q[3];
sx q[3];
rz(-3.1351922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.521495) q[2];
sx q[2];
rz(-1.5767117) q[2];
sx q[2];
rz(-0.0077956789) q[2];
rz(-1.1826285) q[3];
sx q[3];
rz(-1.7595485) q[3];
sx q[3];
rz(-1.2953322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.9919392) q[0];
sx q[0];
rz(-0.45435926) q[0];
sx q[0];
rz(-0.40865189) q[0];
rz(-1.0952134) q[1];
sx q[1];
rz(-1.6588255) q[1];
sx q[1];
rz(-0.85831395) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3811144) q[0];
sx q[0];
rz(-1.71813) q[0];
sx q[0];
rz(3.0744746) q[0];
rz(-pi) q[1];
rz(-1.3038396) q[2];
sx q[2];
rz(-0.46707312) q[2];
sx q[2];
rz(2.2427223) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6950001) q[1];
sx q[1];
rz(-1.4607251) q[1];
sx q[1];
rz(-1.7677444) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32118843) q[3];
sx q[3];
rz(-2.5192755) q[3];
sx q[3];
rz(-0.27143196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3506763) q[2];
sx q[2];
rz(-0.75296593) q[2];
sx q[2];
rz(-2.47993) q[2];
rz(0.800313) q[3];
sx q[3];
rz(-1.9789663) q[3];
sx q[3];
rz(-0.77320981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.4270571) q[0];
rz(0.47572687) q[1];
sx q[1];
rz(-2.5560502) q[1];
sx q[1];
rz(-0.48019662) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4775608) q[0];
sx q[0];
rz(-0.98958221) q[0];
sx q[0];
rz(2.1607735) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91861208) q[2];
sx q[2];
rz(-1.8377853) q[2];
sx q[2];
rz(0.3912386) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0968761) q[1];
sx q[1];
rz(-2.4164817) q[1];
sx q[1];
rz(0.43885751) q[1];
rz(-pi) q[2];
rz(2.0065745) q[3];
sx q[3];
rz(-1.3131071) q[3];
sx q[3];
rz(-2.1439752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6314038) q[2];
sx q[2];
rz(-1.1978585) q[2];
sx q[2];
rz(-2.8655444) q[2];
rz(-0.43802842) q[3];
sx q[3];
rz(-1.5166401) q[3];
sx q[3];
rz(-0.62694222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.20090564) q[0];
sx q[0];
rz(-1.154366) q[0];
sx q[0];
rz(0.05703297) q[0];
rz(-0.036046473) q[1];
sx q[1];
rz(-1.6135975) q[1];
sx q[1];
rz(1.5560908) q[1];
rz(-2.0519999) q[2];
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
