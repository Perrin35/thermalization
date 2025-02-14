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
rz(1.243408) q[0];
sx q[0];
rz(-1.5771447) q[0];
sx q[0];
rz(0.91140437) q[0];
rz(2.4721594) q[1];
sx q[1];
rz(-2.8652006) q[1];
sx q[1];
rz(-0.82958329) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54691168) q[0];
sx q[0];
rz(-1.24238) q[0];
sx q[0];
rz(-1.9921147) q[0];
x q[1];
rz(-1.6638512) q[2];
sx q[2];
rz(-1.3246127) q[2];
sx q[2];
rz(0.024178084) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.098886996) q[1];
sx q[1];
rz(-1.2707842) q[1];
sx q[1];
rz(0.76619958) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8220956) q[3];
sx q[3];
rz(-0.45630672) q[3];
sx q[3];
rz(2.1647477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.36898819) q[2];
sx q[2];
rz(-2.172894) q[2];
sx q[2];
rz(-1.9770835) q[2];
rz(0.72689593) q[3];
sx q[3];
rz(-1.3317069) q[3];
sx q[3];
rz(0.070778457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1233391) q[0];
sx q[0];
rz(-2.4788719) q[0];
sx q[0];
rz(-0.29139274) q[0];
rz(2.5413051) q[1];
sx q[1];
rz(-0.95708668) q[1];
sx q[1];
rz(-2.9765863) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24848973) q[0];
sx q[0];
rz(-1.6418253) q[0];
sx q[0];
rz(1.3079337) q[0];
rz(3.0290452) q[2];
sx q[2];
rz(-2.0172545) q[2];
sx q[2];
rz(3.1242862) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2272039) q[1];
sx q[1];
rz(-1.640135) q[1];
sx q[1];
rz(2.4607865) q[1];
x q[2];
rz(-0.21221186) q[3];
sx q[3];
rz(-2.8877875) q[3];
sx q[3];
rz(-2.5830944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0773641) q[2];
sx q[2];
rz(-0.95244971) q[2];
sx q[2];
rz(-0.23951086) q[2];
rz(-0.32660487) q[3];
sx q[3];
rz(-0.17290792) q[3];
sx q[3];
rz(-3.0467564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8851682) q[0];
sx q[0];
rz(-2.7056077) q[0];
sx q[0];
rz(-1.5727795) q[0];
rz(-0.39372152) q[1];
sx q[1];
rz(-0.55501333) q[1];
sx q[1];
rz(-2.6655925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18506348) q[0];
sx q[0];
rz(-2.1250884) q[0];
sx q[0];
rz(-2.8962062) q[0];
rz(-pi) q[1];
rz(1.7577241) q[2];
sx q[2];
rz(-1.6003967) q[2];
sx q[2];
rz(1.6254978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2326374) q[1];
sx q[1];
rz(-2.401315) q[1];
sx q[1];
rz(-0.93772447) q[1];
rz(-pi) q[2];
rz(2.1952403) q[3];
sx q[3];
rz(-0.40842429) q[3];
sx q[3];
rz(-0.32865903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8126882) q[2];
sx q[2];
rz(-0.71949553) q[2];
sx q[2];
rz(0.94181124) q[2];
rz(-1.4224667) q[3];
sx q[3];
rz(-0.37426451) q[3];
sx q[3];
rz(-2.4678738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068168966) q[0];
sx q[0];
rz(-0.24003679) q[0];
sx q[0];
rz(1.6085251) q[0];
rz(0.34670058) q[1];
sx q[1];
rz(-1.0418714) q[1];
sx q[1];
rz(0.18992058) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1153206) q[0];
sx q[0];
rz(-2.2824725) q[0];
sx q[0];
rz(1.9188966) q[0];
rz(2.0702576) q[2];
sx q[2];
rz(-0.53218666) q[2];
sx q[2];
rz(2.4376873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9393438) q[1];
sx q[1];
rz(-2.1758683) q[1];
sx q[1];
rz(-0.069496198) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40855405) q[3];
sx q[3];
rz(-1.1748649) q[3];
sx q[3];
rz(2.8568639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7803663) q[2];
sx q[2];
rz(-0.85005886) q[2];
sx q[2];
rz(0.38193199) q[2];
rz(3.0319038) q[3];
sx q[3];
rz(-1.0332801) q[3];
sx q[3];
rz(-3.0221353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1031951) q[0];
sx q[0];
rz(-1.9586451) q[0];
sx q[0];
rz(-0.84684816) q[0];
rz(-0.13941828) q[1];
sx q[1];
rz(-2.3250695) q[1];
sx q[1];
rz(0.74904186) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8183621) q[0];
sx q[0];
rz(-1.4233755) q[0];
sx q[0];
rz(-2.4458103) q[0];
rz(-2.8770859) q[2];
sx q[2];
rz(-2.6418984) q[2];
sx q[2];
rz(-0.12910143) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3934259) q[1];
sx q[1];
rz(-1.6977377) q[1];
sx q[1];
rz(2.9462161) q[1];
rz(-pi) q[2];
rz(0.11744515) q[3];
sx q[3];
rz(-2.0694039) q[3];
sx q[3];
rz(-1.487446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.9355115) q[2];
sx q[2];
rz(-1.2568018) q[2];
sx q[2];
rz(2.0743267) q[2];
rz(-2.4308128) q[3];
sx q[3];
rz(-1.658541) q[3];
sx q[3];
rz(1.425364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7988605) q[0];
sx q[0];
rz(-1.4787759) q[0];
sx q[0];
rz(-2.1730098) q[0];
rz(-0.69106483) q[1];
sx q[1];
rz(-1.3422809) q[1];
sx q[1];
rz(-1.8341281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3569865) q[0];
sx q[0];
rz(-0.95080599) q[0];
sx q[0];
rz(1.8852573) q[0];
rz(-0.5796104) q[2];
sx q[2];
rz(-2.7977849) q[2];
sx q[2];
rz(2.0909617) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.71364129) q[1];
sx q[1];
rz(-0.390807) q[1];
sx q[1];
rz(-1.7094088) q[1];
x q[2];
rz(-2.1627681) q[3];
sx q[3];
rz(-2.8205964) q[3];
sx q[3];
rz(2.5227566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6664194) q[2];
sx q[2];
rz(-0.23491493) q[2];
sx q[2];
rz(-1.2972181) q[2];
rz(1.746486) q[3];
sx q[3];
rz(-1.3336983) q[3];
sx q[3];
rz(1.6102128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.7995826) q[0];
sx q[0];
rz(-0.81556773) q[0];
sx q[0];
rz(0.94014257) q[0];
rz(2.4932585) q[1];
sx q[1];
rz(-1.9377361) q[1];
sx q[1];
rz(-2.8727093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.927664) q[0];
sx q[0];
rz(-1.0335796) q[0];
sx q[0];
rz(1.8570333) q[0];
rz(-2.778947) q[2];
sx q[2];
rz(-0.73604903) q[2];
sx q[2];
rz(2.6578052) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9098879) q[1];
sx q[1];
rz(-2.1429166) q[1];
sx q[1];
rz(0.27426274) q[1];
rz(-2.4869789) q[3];
sx q[3];
rz(-0.5236434) q[3];
sx q[3];
rz(0.80010993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0073283422) q[2];
sx q[2];
rz(-1.2790044) q[2];
sx q[2];
rz(-2.4556665) q[2];
rz(2.4053597) q[3];
sx q[3];
rz(-2.1015621) q[3];
sx q[3];
rz(0.22549103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69970423) q[0];
sx q[0];
rz(-1.2208953) q[0];
sx q[0];
rz(0.7487444) q[0];
rz(3.0486095) q[1];
sx q[1];
rz(-2.3817606) q[1];
sx q[1];
rz(-2.0704796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7629246) q[0];
sx q[0];
rz(-2.0708186) q[0];
sx q[0];
rz(-1.0944233) q[0];
rz(-pi) q[1];
rz(-0.65437859) q[2];
sx q[2];
rz(-1.9438231) q[2];
sx q[2];
rz(-1.4747064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2805713) q[1];
sx q[1];
rz(-1.8771267) q[1];
sx q[1];
rz(0.88545274) q[1];
x q[2];
rz(1.6900916) q[3];
sx q[3];
rz(-0.41364851) q[3];
sx q[3];
rz(3.1320437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.45212713) q[2];
sx q[2];
rz(-0.63673821) q[2];
sx q[2];
rz(2.9316736) q[2];
rz(3.012015) q[3];
sx q[3];
rz(-0.94730535) q[3];
sx q[3];
rz(-1.5773704) q[3];
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
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0196446) q[0];
sx q[0];
rz(-2.8804998) q[0];
sx q[0];
rz(-3.1280532) q[0];
rz(2.6059222) q[1];
sx q[1];
rz(-0.56709138) q[1];
sx q[1];
rz(-0.013669107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22937742) q[0];
sx q[0];
rz(-1.3446864) q[0];
sx q[0];
rz(-0.78271981) q[0];
x q[1];
rz(2.2315908) q[2];
sx q[2];
rz(-2.5507413) q[2];
sx q[2];
rz(2.1995418) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0422349) q[1];
sx q[1];
rz(-0.7868979) q[1];
sx q[1];
rz(2.6762647) q[1];
x q[2];
rz(0.39710703) q[3];
sx q[3];
rz(-0.79063168) q[3];
sx q[3];
rz(1.9472029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6232264) q[2];
sx q[2];
rz(-1.3806815) q[2];
sx q[2];
rz(-1.7784485) q[2];
rz(-0.81950435) q[3];
sx q[3];
rz(-0.47911152) q[3];
sx q[3];
rz(0.68377408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226456) q[0];
sx q[0];
rz(-2.2830257) q[0];
sx q[0];
rz(1.488142) q[0];
rz(-0.34291521) q[1];
sx q[1];
rz(-1.5838793) q[1];
sx q[1];
rz(2.3905579) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2011178) q[0];
sx q[0];
rz(-1.6342666) q[0];
sx q[0];
rz(1.8497784) q[0];
rz(-pi) q[1];
rz(-1.1079587) q[2];
sx q[2];
rz(-1.7506536) q[2];
sx q[2];
rz(-3.0015869) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.064709238) q[1];
sx q[1];
rz(-2.5297477) q[1];
sx q[1];
rz(-0.20670273) q[1];
x q[2];
rz(-1.4806101) q[3];
sx q[3];
rz(-1.5377561) q[3];
sx q[3];
rz(0.15109135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.670383) q[2];
sx q[2];
rz(-0.94380108) q[2];
sx q[2];
rz(0.38718265) q[2];
rz(-2.6662628) q[3];
sx q[3];
rz(-1.1673704) q[3];
sx q[3];
rz(-0.79132426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.249007) q[0];
sx q[0];
rz(-0.96207608) q[0];
sx q[0];
rz(0.67208653) q[0];
rz(-2.7723906) q[1];
sx q[1];
rz(-0.75405706) q[1];
sx q[1];
rz(-1.3934607) q[1];
rz(-0.21239352) q[2];
sx q[2];
rz(-0.72827374) q[2];
sx q[2];
rz(-0.49401415) q[2];
rz(1.8440856) q[3];
sx q[3];
rz(-1.4617625) q[3];
sx q[3];
rz(1.1980496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
