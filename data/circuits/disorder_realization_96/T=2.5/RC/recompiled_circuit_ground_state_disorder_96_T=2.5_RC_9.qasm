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
rz(2.2267377) q[0];
rz(2.272361) q[1];
sx q[1];
rz(-0.72328049) q[1];
sx q[1];
rz(-1.7887315) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3769414) q[0];
sx q[0];
rz(-0.94263715) q[0];
sx q[0];
rz(-2.2215722) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74341821) q[2];
sx q[2];
rz(-1.2890491) q[2];
sx q[2];
rz(-2.8670058) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5946002) q[1];
sx q[1];
rz(-2.2359097) q[1];
sx q[1];
rz(0.2427264) q[1];
rz(-pi) q[2];
rz(-2.1800024) q[3];
sx q[3];
rz(-1.7477027) q[3];
sx q[3];
rz(-2.3576403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73504084) q[2];
sx q[2];
rz(-2.328379) q[2];
sx q[2];
rz(1.0068033) q[2];
rz(1.6793647) q[3];
sx q[3];
rz(-1.4454602) q[3];
sx q[3];
rz(-1.538895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8524858) q[0];
sx q[0];
rz(-2.2302581) q[0];
sx q[0];
rz(-3.015633) q[0];
rz(2.0055298) q[1];
sx q[1];
rz(-1.2964396) q[1];
sx q[1];
rz(-2.0819285) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7505109) q[0];
sx q[0];
rz(-1.5727057) q[0];
sx q[0];
rz(1.5715412) q[0];
rz(2.4074209) q[2];
sx q[2];
rz(-2.3245272) q[2];
sx q[2];
rz(1.2881607) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30970022) q[1];
sx q[1];
rz(-0.87617481) q[1];
sx q[1];
rz(-2.5599856) q[1];
x q[2];
rz(-2.7823388) q[3];
sx q[3];
rz(-1.6174966) q[3];
sx q[3];
rz(0.31622313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5560567) q[2];
sx q[2];
rz(-1.5289331) q[2];
sx q[2];
rz(-2.5118828) q[2];
rz(-1.2398047) q[3];
sx q[3];
rz(-2.3972378) q[3];
sx q[3];
rz(1.121421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12104163) q[0];
sx q[0];
rz(-0.23153767) q[0];
sx q[0];
rz(1.2303906) q[0];
rz(2.4379099) q[1];
sx q[1];
rz(-0.92862248) q[1];
sx q[1];
rz(-1.5812662) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15250041) q[0];
sx q[0];
rz(-1.1772424) q[0];
sx q[0];
rz(3.1253184) q[0];
rz(1.376344) q[2];
sx q[2];
rz(-2.299435) q[2];
sx q[2];
rz(-3.0448845) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.81654233) q[1];
sx q[1];
rz(-0.61707622) q[1];
sx q[1];
rz(0.80775921) q[1];
rz(-pi) q[2];
rz(2.7625691) q[3];
sx q[3];
rz(-0.82250094) q[3];
sx q[3];
rz(-1.3027954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33502093) q[2];
sx q[2];
rz(-1.0806012) q[2];
sx q[2];
rz(-3.0873155) q[2];
rz(1.0862167) q[3];
sx q[3];
rz(-0.79038668) q[3];
sx q[3];
rz(2.6495892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739968) q[0];
sx q[0];
rz(-0.01570276) q[0];
sx q[0];
rz(-2.1570461) q[0];
rz(1.4060075) q[1];
sx q[1];
rz(-1.5407298) q[1];
sx q[1];
rz(0.23449177) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67574745) q[0];
sx q[0];
rz(-0.75006163) q[0];
sx q[0];
rz(-1.1613599) q[0];
x q[1];
rz(1.8580397) q[2];
sx q[2];
rz(-1.3352639) q[2];
sx q[2];
rz(-1.6446496) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.381466) q[1];
sx q[1];
rz(-1.0262118) q[1];
sx q[1];
rz(-1.8167956) q[1];
rz(-pi) q[2];
rz(-3.111895) q[3];
sx q[3];
rz(-2.1468024) q[3];
sx q[3];
rz(2.821049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8977114) q[2];
sx q[2];
rz(-0.66456866) q[2];
sx q[2];
rz(0.74771869) q[2];
rz(1.0866577) q[3];
sx q[3];
rz(-0.99435157) q[3];
sx q[3];
rz(2.7896816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9413167) q[0];
sx q[0];
rz(-0.87012297) q[0];
sx q[0];
rz(1.548832) q[0];
rz(-0.13762711) q[1];
sx q[1];
rz(-2.177114) q[1];
sx q[1];
rz(0.67759883) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9842868) q[0];
sx q[0];
rz(-1.56347) q[0];
sx q[0];
rz(3.1276032) q[0];
rz(2.9477411) q[2];
sx q[2];
rz(-1.7192013) q[2];
sx q[2];
rz(1.0184763) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4487344) q[1];
sx q[1];
rz(-0.13693196) q[1];
sx q[1];
rz(-1.5278234) q[1];
rz(-pi) q[2];
rz(-2.4047071) q[3];
sx q[3];
rz(-1.7523271) q[3];
sx q[3];
rz(1.5893861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4981093) q[2];
sx q[2];
rz(-0.52916873) q[2];
sx q[2];
rz(2.7672178) q[2];
rz(-0.74786413) q[3];
sx q[3];
rz(-1.4948083) q[3];
sx q[3];
rz(2.0428467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0102274) q[0];
sx q[0];
rz(-1.7759198) q[0];
sx q[0];
rz(0.22450547) q[0];
rz(2.7911216) q[1];
sx q[1];
rz(-1.1843362) q[1];
sx q[1];
rz(1.9559466) q[1];
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
rz(1.7189222) q[0];
x q[1];
rz(1.6669438) q[2];
sx q[2];
rz(-1.1091057) q[2];
sx q[2];
rz(-0.27335121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.396186) q[1];
sx q[1];
rz(-2.2377659) q[1];
sx q[1];
rz(-2.1275669) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3697183) q[3];
sx q[3];
rz(-0.83362416) q[3];
sx q[3];
rz(1.9825359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4993569) q[2];
sx q[2];
rz(-2.2580052) q[2];
sx q[2];
rz(0.16481608) q[2];
rz(-0.76588255) q[3];
sx q[3];
rz(-0.82847786) q[3];
sx q[3];
rz(-0.59144056) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0484475) q[0];
sx q[0];
rz(-0.026938139) q[0];
sx q[0];
rz(2.7222166) q[0];
rz(1.558149) q[1];
sx q[1];
rz(-2.7132468) q[1];
sx q[1];
rz(-1.8289808) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2079577) q[0];
sx q[0];
rz(-2.5267753) q[0];
sx q[0];
rz(-1.2864248) q[0];
rz(-pi) q[1];
rz(0.68676853) q[2];
sx q[2];
rz(-1.0305974) q[2];
sx q[2];
rz(1.5093925) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53322843) q[1];
sx q[1];
rz(-1.3238412) q[1];
sx q[1];
rz(-0.19503959) q[1];
rz(-2.6978771) q[3];
sx q[3];
rz(-2.8464918) q[3];
sx q[3];
rz(2.091311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4075564) q[2];
sx q[2];
rz(-2.4702256) q[2];
sx q[2];
rz(-2.7835795) q[2];
rz(-2.9426306) q[3];
sx q[3];
rz(-2.3301061) q[3];
sx q[3];
rz(-3.0437886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99893779) q[0];
sx q[0];
rz(-2.1079347) q[0];
sx q[0];
rz(-2.348483) q[0];
rz(-0.89590394) q[1];
sx q[1];
rz(-1.221012) q[1];
sx q[1];
rz(-2.6633247) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90253622) q[0];
sx q[0];
rz(-2.0180114) q[0];
sx q[0];
rz(1.4330401) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7061126) q[2];
sx q[2];
rz(-1.9009095) q[2];
sx q[2];
rz(-0.94089905) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.079887159) q[1];
sx q[1];
rz(-1.8338039) q[1];
sx q[1];
rz(1.0884465) q[1];
rz(-2.0860772) q[3];
sx q[3];
rz(-1.6664423) q[3];
sx q[3];
rz(-3.1351922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6200977) q[2];
sx q[2];
rz(-1.564881) q[2];
sx q[2];
rz(-0.0077956789) q[2];
rz(-1.1826285) q[3];
sx q[3];
rz(-1.7595485) q[3];
sx q[3];
rz(1.8462605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9919392) q[0];
sx q[0];
rz(-0.45435926) q[0];
sx q[0];
rz(-2.7329408) q[0];
rz(-2.0463792) q[1];
sx q[1];
rz(-1.4827671) q[1];
sx q[1];
rz(2.2832787) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33108157) q[0];
sx q[0];
rz(-0.16180049) q[0];
sx q[0];
rz(1.9952379) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1180664) q[2];
sx q[2];
rz(-1.6898587) q[2];
sx q[2];
rz(0.43242142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6950001) q[1];
sx q[1];
rz(-1.6808676) q[1];
sx q[1];
rz(1.3738483) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7935221) q[3];
sx q[3];
rz(-0.98470426) q[3];
sx q[3];
rz(3.0243788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3506763) q[2];
sx q[2];
rz(-0.75296593) q[2];
sx q[2];
rz(-0.6616627) q[2];
rz(-0.800313) q[3];
sx q[3];
rz(-1.9789663) q[3];
sx q[3];
rz(0.77320981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3043587) q[0];
sx q[0];
rz(-2.9999314) q[0];
sx q[0];
rz(-0.71453553) q[0];
rz(2.6658658) q[1];
sx q[1];
rz(-2.5560502) q[1];
sx q[1];
rz(-2.661396) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4775608) q[0];
sx q[0];
rz(-2.1520104) q[0];
sx q[0];
rz(-0.98081918) q[0];
rz(1.1473898) q[2];
sx q[2];
rz(-2.4443279) q[2];
sx q[2];
rz(-0.84691511) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.678315) q[1];
sx q[1];
rz(-0.92683219) q[1];
sx q[1];
rz(-1.2106845) q[1];
rz(-1.0126635) q[3];
sx q[3];
rz(-2.6395661) q[3];
sx q[3];
rz(-3.0691911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51018888) q[2];
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
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.940687) q[0];
sx q[0];
rz(-1.154366) q[0];
sx q[0];
rz(0.05703297) q[0];
rz(-3.1055462) q[1];
sx q[1];
rz(-1.5279952) q[1];
sx q[1];
rz(-1.5855018) q[1];
rz(1.0771607) q[2];
sx q[2];
rz(-1.8133327) q[2];
sx q[2];
rz(1.2015154) q[2];
rz(1.4344169) q[3];
sx q[3];
rz(-1.8227259) q[3];
sx q[3];
rz(-1.3397761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
