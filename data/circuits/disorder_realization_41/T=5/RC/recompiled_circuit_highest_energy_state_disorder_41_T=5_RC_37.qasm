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
rz(2.3693585) q[0];
sx q[0];
rz(-1.9263664) q[0];
sx q[0];
rz(-1.9908494) q[0];
rz(2.8590705) q[1];
sx q[1];
rz(-1.1257659) q[1];
sx q[1];
rz(2.0292625) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.756506) q[0];
sx q[0];
rz(-1.6876564) q[0];
sx q[0];
rz(-2.8782842) q[0];
rz(-pi) q[1];
rz(3.0289828) q[2];
sx q[2];
rz(-2.1499584) q[2];
sx q[2];
rz(-2.4899768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67206269) q[1];
sx q[1];
rz(-2.1066252) q[1];
sx q[1];
rz(-2.5754926) q[1];
rz(-pi) q[2];
rz(0.78948094) q[3];
sx q[3];
rz(-1.0320569) q[3];
sx q[3];
rz(1.0305424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34899601) q[2];
sx q[2];
rz(-2.0670321) q[2];
sx q[2];
rz(-1.7462771) q[2];
rz(-0.058517728) q[3];
sx q[3];
rz(-1.4165712) q[3];
sx q[3];
rz(-1.3099366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7175452) q[0];
sx q[0];
rz(-0.40322867) q[0];
sx q[0];
rz(-0.44977093) q[0];
rz(2.6768118) q[1];
sx q[1];
rz(-1.8315146) q[1];
sx q[1];
rz(-2.8890077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3656986) q[0];
sx q[0];
rz(-1.5868385) q[0];
sx q[0];
rz(1.2195862) q[0];
rz(-pi) q[1];
rz(-0.66030432) q[2];
sx q[2];
rz(-1.2129158) q[2];
sx q[2];
rz(-1.344327) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77923939) q[1];
sx q[1];
rz(-1.5679411) q[1];
sx q[1];
rz(1.289202) q[1];
rz(-pi) q[2];
rz(0.81124311) q[3];
sx q[3];
rz(-1.1170804) q[3];
sx q[3];
rz(-2.8123901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8978591) q[2];
sx q[2];
rz(-2.2458138) q[2];
sx q[2];
rz(3.1254357) q[2];
rz(0.92205087) q[3];
sx q[3];
rz(-2.1478896) q[3];
sx q[3];
rz(3.0097358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048412662) q[0];
sx q[0];
rz(-1.484363) q[0];
sx q[0];
rz(-2.4600001) q[0];
rz(-0.59773481) q[1];
sx q[1];
rz(-1.6860516) q[1];
sx q[1];
rz(-2.8173503) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77406835) q[0];
sx q[0];
rz(-1.2307967) q[0];
sx q[0];
rz(0.46214477) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6382255) q[2];
sx q[2];
rz(-2.0766269) q[2];
sx q[2];
rz(-2.1492599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5432918) q[1];
sx q[1];
rz(-0.17379119) q[1];
sx q[1];
rz(-1.325118) q[1];
x q[2];
rz(-2.1799654) q[3];
sx q[3];
rz(-1.7097843) q[3];
sx q[3];
rz(-2.9442996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2250259) q[2];
sx q[2];
rz(-0.39174199) q[2];
sx q[2];
rz(-2.0908835) q[2];
rz(-2.4010036) q[3];
sx q[3];
rz(-1.8480443) q[3];
sx q[3];
rz(0.061323969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.4968313) q[0];
sx q[0];
rz(-0.83503857) q[0];
sx q[0];
rz(-0.9541676) q[0];
rz(1.1369811) q[1];
sx q[1];
rz(-1.515772) q[1];
sx q[1];
rz(-2.383393) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2187463) q[0];
sx q[0];
rz(-0.2740261) q[0];
sx q[0];
rz(1.42204) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9304843) q[2];
sx q[2];
rz(-0.75497113) q[2];
sx q[2];
rz(1.3364707) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8710701) q[1];
sx q[1];
rz(-2.4409373) q[1];
sx q[1];
rz(1.4187705) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5605035) q[3];
sx q[3];
rz(-0.83196466) q[3];
sx q[3];
rz(-2.9266299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0884462) q[2];
sx q[2];
rz(-0.55067486) q[2];
sx q[2];
rz(1.3147563) q[2];
rz(2.8407319) q[3];
sx q[3];
rz(-1.8010537) q[3];
sx q[3];
rz(-0.68527591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.143173) q[0];
sx q[0];
rz(-1.5851861) q[0];
sx q[0];
rz(1.5013303) q[0];
rz(-2.7550664) q[1];
sx q[1];
rz(-2.0730348) q[1];
sx q[1];
rz(1.5466669) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37171897) q[0];
sx q[0];
rz(-1.7434911) q[0];
sx q[0];
rz(-2.4358389) q[0];
rz(-pi) q[1];
rz(0.96225454) q[2];
sx q[2];
rz(-1.8402765) q[2];
sx q[2];
rz(-2.7121921) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9846696) q[1];
sx q[1];
rz(-2.700173) q[1];
sx q[1];
rz(-0.19766165) q[1];
rz(-pi) q[2];
rz(1.5210454) q[3];
sx q[3];
rz(-2.1440268) q[3];
sx q[3];
rz(1.075268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74653643) q[2];
sx q[2];
rz(-0.48607963) q[2];
sx q[2];
rz(2.0716095) q[2];
rz(-2.2561031) q[3];
sx q[3];
rz(-2.077379) q[3];
sx q[3];
rz(1.989367) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2162061) q[0];
sx q[0];
rz(-0.33846551) q[0];
sx q[0];
rz(-2.0881407) q[0];
rz(-2.9523051) q[1];
sx q[1];
rz(-0.82866755) q[1];
sx q[1];
rz(-1.4922356) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20331377) q[0];
sx q[0];
rz(-0.39394475) q[0];
sx q[0];
rz(-1.2039494) q[0];
rz(-1.8771421) q[2];
sx q[2];
rz(-1.7736846) q[2];
sx q[2];
rz(0.63837793) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0227532) q[1];
sx q[1];
rz(-0.51278881) q[1];
sx q[1];
rz(0.20365479) q[1];
rz(1.6118649) q[3];
sx q[3];
rz(-2.3784935) q[3];
sx q[3];
rz(-1.3535318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3377127) q[2];
sx q[2];
rz(-2.1264075) q[2];
sx q[2];
rz(-2.0886776) q[2];
rz(0.11689154) q[3];
sx q[3];
rz(-2.0785073) q[3];
sx q[3];
rz(-1.4611999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18148024) q[0];
sx q[0];
rz(-2.5616665) q[0];
sx q[0];
rz(-0.54452288) q[0];
rz(2.8879884) q[1];
sx q[1];
rz(-0.2519775) q[1];
sx q[1];
rz(0.25467083) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9601128) q[0];
sx q[0];
rz(-1.5770265) q[0];
sx q[0];
rz(-1.64272) q[0];
x q[1];
rz(-1.6572701) q[2];
sx q[2];
rz(-1.9625798) q[2];
sx q[2];
rz(2.0287152) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0905625) q[1];
sx q[1];
rz(-1.6335618) q[1];
sx q[1];
rz(0.81845567) q[1];
x q[2];
rz(-2.3319844) q[3];
sx q[3];
rz(-1.1407462) q[3];
sx q[3];
rz(2.3187217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2250681) q[2];
sx q[2];
rz(-0.93776408) q[2];
sx q[2];
rz(-0.41898215) q[2];
rz(-0.032912832) q[3];
sx q[3];
rz(-1.907828) q[3];
sx q[3];
rz(-2.029443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5527749) q[0];
sx q[0];
rz(-0.7426312) q[0];
sx q[0];
rz(-1.8137929) q[0];
rz(2.1090419) q[1];
sx q[1];
rz(-1.3464709) q[1];
sx q[1];
rz(-0.58951497) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96846554) q[0];
sx q[0];
rz(-2.0636981) q[0];
sx q[0];
rz(-2.1987134) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74331626) q[2];
sx q[2];
rz(-2.7545815) q[2];
sx q[2];
rz(2.1402556) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8679598) q[1];
sx q[1];
rz(-2.3936749) q[1];
sx q[1];
rz(-0.91435097) q[1];
rz(-0.91752136) q[3];
sx q[3];
rz(-0.60169501) q[3];
sx q[3];
rz(1.2334005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0699658) q[2];
sx q[2];
rz(-2.4667141) q[2];
sx q[2];
rz(-0.91342941) q[2];
rz(-2.6895798) q[3];
sx q[3];
rz(-1.2830696) q[3];
sx q[3];
rz(1.7970596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87208676) q[0];
sx q[0];
rz(-2.1156613) q[0];
sx q[0];
rz(1.6312697) q[0];
rz(1.3144846) q[1];
sx q[1];
rz(-1.038082) q[1];
sx q[1];
rz(0.41183919) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3519234) q[0];
sx q[0];
rz(-0.69286197) q[0];
sx q[0];
rz(-1.1293651) q[0];
rz(-0.14526315) q[2];
sx q[2];
rz(-2.8194339) q[2];
sx q[2];
rz(-1.3346498) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2862368) q[1];
sx q[1];
rz(-2.0986631) q[1];
sx q[1];
rz(-1.1638454) q[1];
rz(1.8985073) q[3];
sx q[3];
rz(-1.2043038) q[3];
sx q[3];
rz(-1.9043363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0044378) q[2];
sx q[2];
rz(-1.6496481) q[2];
sx q[2];
rz(-0.23492661) q[2];
rz(-2.2233985) q[3];
sx q[3];
rz(-2.4166959) q[3];
sx q[3];
rz(1.0880067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10201564) q[0];
sx q[0];
rz(-0.35440847) q[0];
sx q[0];
rz(-1.8812195) q[0];
rz(1.6873987) q[1];
sx q[1];
rz(-1.6834384) q[1];
sx q[1];
rz(1.4448602) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4871976) q[0];
sx q[0];
rz(-2.5154468) q[0];
sx q[0];
rz(-2.3930413) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2606662) q[2];
sx q[2];
rz(-1.147715) q[2];
sx q[2];
rz(-2.1264875) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9493172) q[1];
sx q[1];
rz(-0.85964627) q[1];
sx q[1];
rz(-1.4059195) q[1];
rz(-0.71717307) q[3];
sx q[3];
rz(-2.6568718) q[3];
sx q[3];
rz(-2.6074099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5789648) q[2];
sx q[2];
rz(-2.8426888) q[2];
sx q[2];
rz(-3.01801) q[2];
rz(0.28892162) q[3];
sx q[3];
rz(-1.8653423) q[3];
sx q[3];
rz(0.45946521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73744437) q[0];
sx q[0];
rz(-1.3442208) q[0];
sx q[0];
rz(-1.5303045) q[0];
rz(0.96724802) q[1];
sx q[1];
rz(-2.1687242) q[1];
sx q[1];
rz(0.8868934) q[1];
rz(-1.5161472) q[2];
sx q[2];
rz(-0.9722295) q[2];
sx q[2];
rz(0.35884919) q[2];
rz(-0.88913118) q[3];
sx q[3];
rz(-2.0778772) q[3];
sx q[3];
rz(-2.7419013) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
