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
rz(2.4438357) q[0];
sx q[0];
rz(-1.9480167) q[0];
sx q[0];
rz(0.20456631) q[0];
rz(2.3808631) q[1];
sx q[1];
rz(-1.7525571) q[1];
sx q[1];
rz(-1.3995481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9153109) q[0];
sx q[0];
rz(-1.4797512) q[0];
sx q[0];
rz(1.7928726) q[0];
rz(0.84487652) q[2];
sx q[2];
rz(-2.5683218) q[2];
sx q[2];
rz(2.5115761) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7508937) q[1];
sx q[1];
rz(-1.1502153) q[1];
sx q[1];
rz(1.4750359) q[1];
rz(-pi) q[2];
rz(0.70943009) q[3];
sx q[3];
rz(-1.7755847) q[3];
sx q[3];
rz(0.11573175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43510398) q[2];
sx q[2];
rz(-3.1092643) q[2];
sx q[2];
rz(2.6522563) q[2];
rz(2.1183744) q[3];
sx q[3];
rz(-0.018298572) q[3];
sx q[3];
rz(1.1255012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045687549) q[0];
sx q[0];
rz(-2.4850595) q[0];
sx q[0];
rz(-0.80279654) q[0];
rz(3.070201) q[1];
sx q[1];
rz(-2.8749021) q[1];
sx q[1];
rz(0.057770483) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18081576) q[0];
sx q[0];
rz(-2.1102) q[0];
sx q[0];
rz(-2.8404854) q[0];
x q[1];
rz(2.7947756) q[2];
sx q[2];
rz(-2.8379734) q[2];
sx q[2];
rz(0.70250073) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7747468) q[1];
sx q[1];
rz(-0.54301942) q[1];
sx q[1];
rz(-1.1096891) q[1];
x q[2];
rz(0.85912786) q[3];
sx q[3];
rz(-1.458408) q[3];
sx q[3];
rz(-0.64765644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1959261) q[2];
sx q[2];
rz(-2.0464996) q[2];
sx q[2];
rz(-1.2954953) q[2];
rz(0.96674353) q[3];
sx q[3];
rz(-2.3714378) q[3];
sx q[3];
rz(-0.75743341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029723786) q[0];
sx q[0];
rz(-1.7787378) q[0];
sx q[0];
rz(1.7080074) q[0];
rz(3.0729821) q[1];
sx q[1];
rz(-1.5741293) q[1];
sx q[1];
rz(0.57919085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96087181) q[0];
sx q[0];
rz(-1.5128883) q[0];
sx q[0];
rz(-1.6625893) q[0];
x q[1];
rz(-3.1379329) q[2];
sx q[2];
rz(-2.6658667) q[2];
sx q[2];
rz(1.9329247) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3350796) q[1];
sx q[1];
rz(-1.6066222) q[1];
sx q[1];
rz(-2.9153009) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0491568) q[3];
sx q[3];
rz(-1.1598827) q[3];
sx q[3];
rz(0.26518566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27089831) q[2];
sx q[2];
rz(-1.8879994) q[2];
sx q[2];
rz(-2.9581621) q[2];
rz(-0.80426788) q[3];
sx q[3];
rz(-0.99884123) q[3];
sx q[3];
rz(-2.6075294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428228) q[0];
sx q[0];
rz(-3.0267921) q[0];
sx q[0];
rz(-2.542069) q[0];
rz(-2.7291258) q[1];
sx q[1];
rz(-0.020655276) q[1];
sx q[1];
rz(2.1309158) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090792716) q[0];
sx q[0];
rz(-2.7677892) q[0];
sx q[0];
rz(-0.35275491) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2029271) q[2];
sx q[2];
rz(-1.6118587) q[2];
sx q[2];
rz(1.2031885) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4356723) q[1];
sx q[1];
rz(-1.8413891) q[1];
sx q[1];
rz(0.95939221) q[1];
x q[2];
rz(-2.9948576) q[3];
sx q[3];
rz(-0.76125604) q[3];
sx q[3];
rz(-0.066731922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3707054) q[2];
sx q[2];
rz(-2.7949896) q[2];
sx q[2];
rz(-0.31751219) q[2];
rz(2.5192449) q[3];
sx q[3];
rz(-2.205866) q[3];
sx q[3];
rz(2.5573825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7790826) q[0];
sx q[0];
rz(-0.91479397) q[0];
sx q[0];
rz(-2.1146178) q[0];
rz(2.5912071) q[1];
sx q[1];
rz(-0.06409476) q[1];
sx q[1];
rz(1.2170353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37792045) q[0];
sx q[0];
rz(-0.91912133) q[0];
sx q[0];
rz(-2.4125189) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0116357) q[2];
sx q[2];
rz(-2.5995147) q[2];
sx q[2];
rz(-2.1810069) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7407806) q[1];
sx q[1];
rz(-0.88140872) q[1];
sx q[1];
rz(-2.9184398) q[1];
x q[2];
rz(-0.36320477) q[3];
sx q[3];
rz(-2.2480559) q[3];
sx q[3];
rz(1.1006608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43651849) q[2];
sx q[2];
rz(-1.3631835) q[2];
sx q[2];
rz(-2.3684033) q[2];
rz(-3.0298722) q[3];
sx q[3];
rz(-1.819928) q[3];
sx q[3];
rz(-0.93938655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47950995) q[0];
sx q[0];
rz(-2.918512) q[0];
sx q[0];
rz(-2.7146085) q[0];
rz(0.93049479) q[1];
sx q[1];
rz(-0.016914802) q[1];
sx q[1];
rz(-0.46447909) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9432968) q[0];
sx q[0];
rz(-1.3751327) q[0];
sx q[0];
rz(0.29384675) q[0];
rz(-pi) q[1];
rz(0.33916766) q[2];
sx q[2];
rz(-0.99771032) q[2];
sx q[2];
rz(0.040405191) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.528094) q[1];
sx q[1];
rz(-1.6069176) q[1];
sx q[1];
rz(-2.7725459) q[1];
rz(-1.6397912) q[3];
sx q[3];
rz(-2.7164) q[3];
sx q[3];
rz(-1.0591648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53133196) q[2];
sx q[2];
rz(-1.320763) q[2];
sx q[2];
rz(-0.28826928) q[2];
rz(2.098295) q[3];
sx q[3];
rz(-0.61560029) q[3];
sx q[3];
rz(-2.402795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4814602) q[0];
sx q[0];
rz(-1.2876502) q[0];
sx q[0];
rz(-0.75755358) q[0];
rz(3.0689012) q[1];
sx q[1];
rz(-3.1158267) q[1];
sx q[1];
rz(3.0933948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5023247) q[0];
sx q[0];
rz(-1.5412962) q[0];
sx q[0];
rz(2.155199) q[0];
rz(0.23400314) q[2];
sx q[2];
rz(-2.1562139) q[2];
sx q[2];
rz(-2.8319179) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.99302247) q[1];
sx q[1];
rz(-2.3766368) q[1];
sx q[1];
rz(0.92924849) q[1];
rz(-pi) q[2];
rz(-3.0890325) q[3];
sx q[3];
rz(-2.0400999) q[3];
sx q[3];
rz(1.8183501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66182071) q[2];
sx q[2];
rz(-1.4996108) q[2];
sx q[2];
rz(3.0721967) q[2];
rz(1.5875459) q[3];
sx q[3];
rz(-2.3543251) q[3];
sx q[3];
rz(2.9230996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.7540392) q[0];
sx q[0];
rz(-1.0559005) q[0];
sx q[0];
rz(-1.7283424) q[0];
rz(0.82855254) q[1];
sx q[1];
rz(-0.041752432) q[1];
sx q[1];
rz(0.54263306) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4601645) q[0];
sx q[0];
rz(-1.4519454) q[0];
sx q[0];
rz(2.9959841) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37658875) q[2];
sx q[2];
rz(-1.9053725) q[2];
sx q[2];
rz(-1.1194816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97510135) q[1];
sx q[1];
rz(-1.5169739) q[1];
sx q[1];
rz(3.0256998) q[1];
x q[2];
rz(1.9493136) q[3];
sx q[3];
rz(-1.1633486) q[3];
sx q[3];
rz(0.70127869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1594306) q[2];
sx q[2];
rz(-0.40951481) q[2];
sx q[2];
rz(-2.8816667) q[2];
rz(-2.121117) q[3];
sx q[3];
rz(-2.8838938) q[3];
sx q[3];
rz(0.79403383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4644153) q[0];
sx q[0];
rz(-2.9554415) q[0];
sx q[0];
rz(-1.6625241) q[0];
rz(-1.5326477) q[1];
sx q[1];
rz(-2.1023991) q[1];
sx q[1];
rz(-2.398568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047516454) q[0];
sx q[0];
rz(-2.0819944) q[0];
sx q[0];
rz(-0.1316977) q[0];
rz(-pi) q[1];
rz(-2.7588927) q[2];
sx q[2];
rz(-1.0596794) q[2];
sx q[2];
rz(1.2148884) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.31752095) q[1];
sx q[1];
rz(-1.5050384) q[1];
sx q[1];
rz(-3.0865248) q[1];
rz(-pi) q[2];
rz(3.0366692) q[3];
sx q[3];
rz(-1.6015953) q[3];
sx q[3];
rz(1.9178176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4674025) q[2];
sx q[2];
rz(-0.82618606) q[2];
sx q[2];
rz(0.80545938) q[2];
rz(1.449466) q[3];
sx q[3];
rz(-1.9129246) q[3];
sx q[3];
rz(2.3384371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8886803) q[0];
sx q[0];
rz(-2.5997933) q[0];
sx q[0];
rz(0.75206494) q[0];
rz(-1.1532785) q[1];
sx q[1];
rz(-2.2572932) q[1];
sx q[1];
rz(0.28336743) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3842938) q[0];
sx q[0];
rz(-1.481825) q[0];
sx q[0];
rz(-2.6020537) q[0];
x q[1];
rz(-2.3226363) q[2];
sx q[2];
rz(-1.6791846) q[2];
sx q[2];
rz(3.0367231) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.99774573) q[1];
sx q[1];
rz(-2.4334868) q[1];
sx q[1];
rz(-2.6216402) q[1];
x q[2];
rz(1.9612736) q[3];
sx q[3];
rz(-0.9538528) q[3];
sx q[3];
rz(0.65333593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2470384) q[2];
sx q[2];
rz(-0.082823195) q[2];
sx q[2];
rz(-1.438633) q[2];
rz(2.8476207) q[3];
sx q[3];
rz(-0.014466244) q[3];
sx q[3];
rz(-2.1132052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033584874) q[0];
sx q[0];
rz(-1.4172194) q[0];
sx q[0];
rz(-1.5237756) q[0];
rz(-0.53957466) q[1];
sx q[1];
rz(-2.3537666) q[1];
sx q[1];
rz(-2.981577) q[1];
rz(0.14847427) q[2];
sx q[2];
rz(-1.9273026) q[2];
sx q[2];
rz(-2.9782563) q[2];
rz(2.5359375) q[3];
sx q[3];
rz(-2.6867821) q[3];
sx q[3];
rz(2.3985779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
