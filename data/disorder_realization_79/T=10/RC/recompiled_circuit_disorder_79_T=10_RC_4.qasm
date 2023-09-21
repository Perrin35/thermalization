OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(7.1516501) q[0];
sx q[0];
rz(9.2317543) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(2.4584682) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1238414) q[0];
sx q[0];
rz(-0.089028247) q[0];
sx q[0];
rz(2.428399) q[0];
x q[1];
rz(2.9002951) q[2];
sx q[2];
rz(-1.9120875) q[2];
sx q[2];
rz(1.2933921) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9477387) q[1];
sx q[1];
rz(-0.96809909) q[1];
sx q[1];
rz(-2.1147453) q[1];
x q[2];
rz(-2.1300415) q[3];
sx q[3];
rz(-2.4992001) q[3];
sx q[3];
rz(1.1005644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(-1.0985628) q[2];
rz(-2.0627608) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-2.0400955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(2.9336477) q[0];
rz(0.57693276) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(1.4651441) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10509051) q[0];
sx q[0];
rz(-2.4173792) q[0];
sx q[0];
rz(-0.0037395517) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0826153) q[2];
sx q[2];
rz(-1.3435875) q[2];
sx q[2];
rz(0.28085923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0456902) q[1];
sx q[1];
rz(-1.1046788) q[1];
sx q[1];
rz(0.63506605) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7705455) q[3];
sx q[3];
rz(-2.3608532) q[3];
sx q[3];
rz(-2.3649529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1229822) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(0.62260735) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(-2.6254568) q[0];
rz(-2.5667045) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(2.3410472) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24414177) q[0];
sx q[0];
rz(-1.4276917) q[0];
sx q[0];
rz(-1.98154) q[0];
rz(-pi) q[1];
rz(2.4750701) q[2];
sx q[2];
rz(-0.98072532) q[2];
sx q[2];
rz(-0.18596622) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87677466) q[1];
sx q[1];
rz(-0.41300981) q[1];
sx q[1];
rz(3.1289711) q[1];
rz(-2.27182) q[3];
sx q[3];
rz(-1.3372034) q[3];
sx q[3];
rz(-1.9529395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4642554) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.1490705) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(-2.676679) q[0];
rz(-2.7930296) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(2.0565313) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61468716) q[0];
sx q[0];
rz(-0.58262107) q[0];
sx q[0];
rz(2.7132062) q[0];
x q[1];
rz(-0.24809804) q[2];
sx q[2];
rz(-1.2144226) q[2];
sx q[2];
rz(0.46596371) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3738721) q[1];
sx q[1];
rz(-1.671119) q[1];
sx q[1];
rz(2.9994681) q[1];
x q[2];
rz(-0.55235858) q[3];
sx q[3];
rz(-0.68867749) q[3];
sx q[3];
rz(2.3104582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(2.4604649) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(2.8779023) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8624449) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(-1.408668) q[0];
rz(-0.43235835) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(-0.98168215) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0504793) q[0];
sx q[0];
rz(-1.0266773) q[0];
sx q[0];
rz(-2.1168461) q[0];
x q[1];
rz(2.1443411) q[2];
sx q[2];
rz(-0.28595668) q[2];
sx q[2];
rz(-0.85948247) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7179467) q[1];
sx q[1];
rz(-1.7182933) q[1];
sx q[1];
rz(-0.46848483) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6793628) q[3];
sx q[3];
rz(-2.541399) q[3];
sx q[3];
rz(1.997228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(2.126746) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(-1.3283407) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4846102) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(0.16528027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76089345) q[0];
sx q[0];
rz(-1.6000415) q[0];
sx q[0];
rz(-0.0951035) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.500962) q[2];
sx q[2];
rz(-2.2572821) q[2];
sx q[2];
rz(-1.3697461) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.936441) q[1];
sx q[1];
rz(-1.1679822) q[1];
sx q[1];
rz(-2.5564297) q[1];
x q[2];
rz(-2.9061756) q[3];
sx q[3];
rz(-2.0093007) q[3];
sx q[3];
rz(2.8188843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(-1.4432663) q[2];
rz(-2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(-2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7845602) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(0.7473942) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-0.73648891) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1897141) q[0];
sx q[0];
rz(-1.5218381) q[0];
sx q[0];
rz(0.017107054) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2320802) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(0.92600694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2939824) q[1];
sx q[1];
rz(-0.81969205) q[1];
sx q[1];
rz(1.3241029) q[1];
rz(-2.6408259) q[3];
sx q[3];
rz(-1.1874677) q[3];
sx q[3];
rz(-2.4367743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(2.6100256) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(-1.8435562) q[0];
rz(2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(-0.92179006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67547922) q[0];
sx q[0];
rz(-0.70444324) q[0];
sx q[0];
rz(2.9324313) q[0];
rz(-1.2577031) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(1.8956172) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8340048) q[1];
sx q[1];
rz(-1.7760135) q[1];
sx q[1];
rz(0.25968857) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3619625) q[3];
sx q[3];
rz(-0.71912557) q[3];
sx q[3];
rz(1.5428839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90199295) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(-1.1716589) q[2];
rz(-1.7840067) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25306025) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(-1.6660447) q[0];
rz(1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(2.4618861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1595759) q[0];
sx q[0];
rz(-1.3681108) q[0];
sx q[0];
rz(2.6152339) q[0];
rz(-0.87551261) q[2];
sx q[2];
rz(-1.7098223) q[2];
sx q[2];
rz(-1.306844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3956086) q[1];
sx q[1];
rz(-1.2320227) q[1];
sx q[1];
rz(-2.4376274) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67919517) q[3];
sx q[3];
rz(-2.0742356) q[3];
sx q[3];
rz(-2.5356843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(2.2909686) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(0.7146548) q[0];
rz(-2.4275298) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(1.9649327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.048621) q[0];
sx q[0];
rz(-1.4074568) q[0];
sx q[0];
rz(-0.049157814) q[0];
rz(-pi) q[1];
rz(0.54729692) q[2];
sx q[2];
rz(-2.1314869) q[2];
sx q[2];
rz(-1.0123569) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0050051) q[1];
sx q[1];
rz(-2.1274381) q[1];
sx q[1];
rz(3.0535166) q[1];
x q[2];
rz(-1.6586187) q[3];
sx q[3];
rz(-2.0801968) q[3];
sx q[3];
rz(-1.0305962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24361336) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(-2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(2.7453616) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(1.8489807) q[2];
sx q[2];
rz(-0.85360151) q[2];
sx q[2];
rz(-3.1384946) q[2];
rz(-0.35076326) q[3];
sx q[3];
rz(-1.5784932) q[3];
sx q[3];
rz(-2.2898883) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
