OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4176183) q[0];
sx q[0];
rz(-1.4899878) q[0];
sx q[0];
rz(2.2111501) q[0];
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(-2.0342483) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410755) q[0];
sx q[0];
rz(-0.76671769) q[0];
sx q[0];
rz(0.47857743) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6644457) q[2];
sx q[2];
rz(-2.2331182) q[2];
sx q[2];
rz(-1.9805816) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0336696) q[1];
sx q[1];
rz(-2.8406997) q[1];
sx q[1];
rz(1.9574907) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7928365) q[3];
sx q[3];
rz(-0.96491279) q[3];
sx q[3];
rz(0.44266686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.063623108) q[2];
sx q[2];
rz(-2.412553) q[2];
sx q[2];
rz(-1.3280274) q[2];
rz(-2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(-0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.48224738) q[0];
sx q[0];
rz(-3.0292065) q[0];
sx q[0];
rz(0.88062084) q[0];
rz(-1.847514) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(-0.81726384) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9416695) q[0];
sx q[0];
rz(-2.7669853) q[0];
sx q[0];
rz(-1.9648212) q[0];
rz(0.21821071) q[2];
sx q[2];
rz(-2.0514224) q[2];
sx q[2];
rz(-2.6618119) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4546928) q[1];
sx q[1];
rz(-2.2001839) q[1];
sx q[1];
rz(1.252974) q[1];
x q[2];
rz(-0.18873429) q[3];
sx q[3];
rz(-1.1047603) q[3];
sx q[3];
rz(2.1621494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2237504) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(-0.36402738) q[2];
rz(-2.1552127) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.36689511) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(-2.1752775) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(-1.6945217) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9054748) q[0];
sx q[0];
rz(-2.9258699) q[0];
sx q[0];
rz(-0.082213684) q[0];
rz(-pi) q[1];
rz(0.082712163) q[2];
sx q[2];
rz(-1.2809922) q[2];
sx q[2];
rz(2.2207584) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9008873) q[1];
sx q[1];
rz(-1.9449688) q[1];
sx q[1];
rz(-2.4918873) q[1];
rz(-pi) q[2];
rz(-2.0855911) q[3];
sx q[3];
rz(-0.47009531) q[3];
sx q[3];
rz(-1.521829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7029999) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(1.1052216) q[2];
rz(2.3953719) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352017) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(1.4659457) q[0];
rz(-0.28494596) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(0.38898653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60880946) q[0];
sx q[0];
rz(-0.36768498) q[0];
sx q[0];
rz(0.68433783) q[0];
x q[1];
rz(2.94611) q[2];
sx q[2];
rz(-2.0494378) q[2];
sx q[2];
rz(-0.10749707) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3868689) q[1];
sx q[1];
rz(-1.1479953) q[1];
sx q[1];
rz(-2.5368284) q[1];
rz(0.22104927) q[3];
sx q[3];
rz(-0.76023686) q[3];
sx q[3];
rz(-1.5628712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7230364) q[2];
sx q[2];
rz(-1.1636461) q[2];
sx q[2];
rz(-0.45670613) q[2];
rz(-1.5151954) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(-0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91963768) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(0.36002457) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(2.6470851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3798824) q[0];
sx q[0];
rz(-0.30711353) q[0];
sx q[0];
rz(0.42219992) q[0];
x q[1];
rz(-1.8525271) q[2];
sx q[2];
rz(-3.0006471) q[2];
sx q[2];
rz(0.81705392) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33489409) q[1];
sx q[1];
rz(-0.94167275) q[1];
sx q[1];
rz(-1.4966399) q[1];
x q[2];
rz(-0.0534119) q[3];
sx q[3];
rz(-2.6969389) q[3];
sx q[3];
rz(-2.2821471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4313724) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(2.8430856) q[2];
rz(0.31202894) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(-2.503094) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1453778) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(-0.19590713) q[0];
rz(-3.1205102) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.235199) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7346749) q[0];
sx q[0];
rz(-1.2376681) q[0];
sx q[0];
rz(-1.8963277) q[0];
rz(-pi) q[1];
rz(-1.9094798) q[2];
sx q[2];
rz(-2.2465696) q[2];
sx q[2];
rz(2.0896926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.27855733) q[1];
sx q[1];
rz(-2.1457991) q[1];
sx q[1];
rz(2.0726191) q[1];
x q[2];
rz(0.23630948) q[3];
sx q[3];
rz(-1.4893388) q[3];
sx q[3];
rz(2.8636275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.49332508) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(2.7098999) q[2];
rz(-1.7290944) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(-0.13599642) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7504808) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(3.0294763) q[0];
rz(2.926459) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(2.0281866) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8047785) q[0];
sx q[0];
rz(-1.3114197) q[0];
sx q[0];
rz(2.74385) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8261975) q[2];
sx q[2];
rz(-1.8176259) q[2];
sx q[2];
rz(2.18404) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0472764) q[1];
sx q[1];
rz(-0.89135209) q[1];
sx q[1];
rz(1.9654771) q[1];
x q[2];
rz(0.21344276) q[3];
sx q[3];
rz(-2.0091972) q[3];
sx q[3];
rz(1.9813117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.001174288) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(0.12602885) q[2];
rz(1.0472939) q[3];
sx q[3];
rz(-1.8208561) q[3];
sx q[3];
rz(-2.4333911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66013181) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(-1.460176) q[0];
rz(1.8966282) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(-1.9326899) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35248127) q[0];
sx q[0];
rz(-1.4915823) q[0];
sx q[0];
rz(-1.3080025) q[0];
rz(2.3102343) q[2];
sx q[2];
rz(-2.0908329) q[2];
sx q[2];
rz(1.780873) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3371256) q[1];
sx q[1];
rz(-2.2837688) q[1];
sx q[1];
rz(1.1285524) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94001694) q[3];
sx q[3];
rz(-1.3925843) q[3];
sx q[3];
rz(-2.7959787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(1.8219927) q[2];
rz(2.54946) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(3.106451) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2329344) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(1.3611025) q[0];
rz(1.2110442) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(2.7499054) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2101333) q[0];
sx q[0];
rz(-1.1362695) q[0];
sx q[0];
rz(1.5772485) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7708771) q[2];
sx q[2];
rz(-2.0152425) q[2];
sx q[2];
rz(-1.0619628) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0605676) q[1];
sx q[1];
rz(-2.276366) q[1];
sx q[1];
rz(-0.4391567) q[1];
rz(-pi) q[2];
rz(1.4976981) q[3];
sx q[3];
rz(-2.0953296) q[3];
sx q[3];
rz(-2.227921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.52035511) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(-1.3367782) q[2];
rz(-0.76198602) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(15/(14*pi)) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(-0.57089943) q[0];
rz(-1.7123429) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(0.16194078) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25625944) q[0];
sx q[0];
rz(-0.67978871) q[0];
sx q[0];
rz(1.1023561) q[0];
rz(-1.5230721) q[2];
sx q[2];
rz(-0.18880162) q[2];
sx q[2];
rz(1.4023086) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3871799) q[1];
sx q[1];
rz(-1.3993235) q[1];
sx q[1];
rz(0.13053723) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5374244) q[3];
sx q[3];
rz(-2.0936692) q[3];
sx q[3];
rz(-3.0748622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95742115) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(-1.7133678) q[2];
rz(1.9421633) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(-1.7709581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4409055) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.3700925) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(-0.11907555) q[2];
sx q[2];
rz(-1.4243813) q[2];
sx q[2];
rz(1.4204155) q[2];
rz(0.42967038) q[3];
sx q[3];
rz(-1.3702787) q[3];
sx q[3];
rz(-2.7313781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
