OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(3.677877) q[0];
sx q[0];
rz(10.372547) q[0];
rz(-1.3287969) q[1];
sx q[1];
rz(-1.8741908) q[1];
sx q[1];
rz(1.0277494) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2078903) q[0];
sx q[0];
rz(-0.37651248) q[0];
sx q[0];
rz(-0.062113751) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1711636) q[2];
sx q[2];
rz(-0.52414775) q[2];
sx q[2];
rz(1.580796) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9103968) q[1];
sx q[1];
rz(-2.3379571) q[1];
sx q[1];
rz(-2.5956144) q[1];
rz(0.58524744) q[3];
sx q[3];
rz(-0.71551502) q[3];
sx q[3];
rz(1.1701208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.7069858) q[2];
sx q[2];
rz(-1.0502846) q[2];
rz(-1.1132647) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(-0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072409078) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(-0.29775277) q[0];
rz(-0.61966664) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(1.108095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28716921) q[0];
sx q[0];
rz(-0.65432917) q[0];
sx q[0];
rz(1.9186583) q[0];
x q[1];
rz(2.7116398) q[2];
sx q[2];
rz(-1.0064831) q[2];
sx q[2];
rz(-1.0831837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.260028) q[1];
sx q[1];
rz(-1.2177055) q[1];
sx q[1];
rz(2.7295223) q[1];
rz(1.3012582) q[3];
sx q[3];
rz(-2.2552935) q[3];
sx q[3];
rz(-1.0052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(0.97529808) q[2];
rz(-2.2235928) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(2.8454034) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.179203) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-0.96484819) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436413) q[0];
sx q[0];
rz(-2.7086341) q[0];
sx q[0];
rz(0.71822449) q[0];
rz(-0.29693551) q[2];
sx q[2];
rz(-1.6116217) q[2];
sx q[2];
rz(-2.9253935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.112903) q[1];
sx q[1];
rz(-1.5064872) q[1];
sx q[1];
rz(1.6921922) q[1];
rz(-pi) q[2];
rz(1.1448907) q[3];
sx q[3];
rz(-2.689038) q[3];
sx q[3];
rz(1.2061314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0470011) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(2.0084521) q[2];
rz(-2.9099693) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8621181) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.7650771) q[0];
rz(-2.6230985) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(-0.24212295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3926485) q[0];
sx q[0];
rz(-2.2376275) q[0];
sx q[0];
rz(1.7783661) q[0];
rz(-2.7119615) q[2];
sx q[2];
rz(-1.5826844) q[2];
sx q[2];
rz(-0.55693835) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0603234) q[1];
sx q[1];
rz(-2.6594866) q[1];
sx q[1];
rz(0.37270765) q[1];
rz(-pi) q[2];
rz(2.8321213) q[3];
sx q[3];
rz(-2.7052393) q[3];
sx q[3];
rz(0.63758028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.018192856) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-0.094853178) q[2];
rz(1.799396) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(-0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(-0.36079303) q[0];
rz(1.3882673) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(-1.1345908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73571262) q[0];
sx q[0];
rz(-2.1942733) q[0];
sx q[0];
rz(0.9491802) q[0];
rz(-pi) q[1];
rz(1.1410494) q[2];
sx q[2];
rz(-1.8998713) q[2];
sx q[2];
rz(2.3171901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0134296) q[1];
sx q[1];
rz(-1.1268106) q[1];
sx q[1];
rz(-1.7432937) q[1];
rz(-pi) q[2];
rz(0.17825019) q[3];
sx q[3];
rz(-2.6986487) q[3];
sx q[3];
rz(2.6435542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0052884) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(-2.5689382) q[2];
rz(-0.92875656) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(1.1675534) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3271493) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(-1.2493398) q[0];
rz(1.918474) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(1.1522326) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845247) q[0];
sx q[0];
rz(-0.10604924) q[0];
sx q[0];
rz(1.690879) q[0];
x q[1];
rz(-0.50243369) q[2];
sx q[2];
rz(-1.2049434) q[2];
sx q[2];
rz(-2.1937214) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7631543) q[1];
sx q[1];
rz(-1.3707146) q[1];
sx q[1];
rz(1.6775908) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96529393) q[3];
sx q[3];
rz(-2.1552342) q[3];
sx q[3];
rz(-2.8402929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6727009) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(-2.1506298) q[2];
rz(-2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794849) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(3.0793072) q[0];
rz(-2.9557872) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(-2.7468162) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6900401) q[0];
sx q[0];
rz(-1.3344904) q[0];
sx q[0];
rz(0.4060181) q[0];
rz(1.2497181) q[2];
sx q[2];
rz(-1.9294538) q[2];
sx q[2];
rz(1.9888339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2784087) q[1];
sx q[1];
rz(-1.3235958) q[1];
sx q[1];
rz(-0.43988887) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23468252) q[3];
sx q[3];
rz(-0.15771401) q[3];
sx q[3];
rz(1.7705256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(2.8685692) q[2];
rz(1.8388883) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(2.8295512) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7438695) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(-1.2851108) q[0];
rz(-1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96269755) q[0];
sx q[0];
rz(-2.6216051) q[0];
sx q[0];
rz(-2.2682701) q[0];
x q[1];
rz(1.6216535) q[2];
sx q[2];
rz(-2.2659781) q[2];
sx q[2];
rz(0.67957544) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33412877) q[1];
sx q[1];
rz(-0.51528105) q[1];
sx q[1];
rz(0.79828243) q[1];
rz(-pi) q[2];
rz(2.0128485) q[3];
sx q[3];
rz(-1.8574517) q[3];
sx q[3];
rz(0.096404508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0354707) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(-1.0236053) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(-1.5203083) q[0];
rz(-2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(-2.3044589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1820113) q[0];
sx q[0];
rz(-0.89996979) q[0];
sx q[0];
rz(1.3120033) q[0];
x q[1];
rz(2.5651921) q[2];
sx q[2];
rz(-1.8216368) q[2];
sx q[2];
rz(0.62945156) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9311034) q[1];
sx q[1];
rz(-1.466202) q[1];
sx q[1];
rz(-1.4180257) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0666396) q[3];
sx q[3];
rz(-1.3131724) q[3];
sx q[3];
rz(-0.50592929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.60124406) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(-2.1458697) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-0.76374617) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(1.0572222) q[0];
rz(3.0341042) q[1];
sx q[1];
rz(-1.8881533) q[1];
sx q[1];
rz(0.9799788) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5532593) q[0];
sx q[0];
rz(-0.13561121) q[0];
sx q[0];
rz(-1.9562264) q[0];
x q[1];
rz(-2.6238742) q[2];
sx q[2];
rz(-2.0132408) q[2];
sx q[2];
rz(0.25416086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2154448) q[1];
sx q[1];
rz(-1.1415592) q[1];
sx q[1];
rz(0.26305671) q[1];
rz(1.2471334) q[3];
sx q[3];
rz(-1.666288) q[3];
sx q[3];
rz(3.132706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3048627) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(-2.1137962) q[2];
rz(-1.3868388) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733611) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(2.3095619) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(-1.8503415) q[2];
sx q[2];
rz(-2.5646979) q[2];
sx q[2];
rz(-0.23451351) q[2];
rz(2.6408623) q[3];
sx q[3];
rz(-1.9095608) q[3];
sx q[3];
rz(0.6469971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
