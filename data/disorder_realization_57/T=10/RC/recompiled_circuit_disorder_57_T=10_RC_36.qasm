OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6141619) q[0];
sx q[0];
rz(-0.80978137) q[0];
sx q[0];
rz(-0.53139395) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(1.3637435) q[1];
sx q[1];
rz(10.627828) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8578313) q[0];
sx q[0];
rz(-2.8537822) q[0];
sx q[0];
rz(-0.83588155) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3538023) q[2];
sx q[2];
rz(-1.8197682) q[2];
sx q[2];
rz(-0.72460246) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0002999) q[1];
sx q[1];
rz(-1.3360026) q[1];
sx q[1];
rz(-1.512212) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8486512) q[3];
sx q[3];
rz(-1.5353234) q[3];
sx q[3];
rz(-0.5637416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8649341) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(0.19908389) q[2];
rz(1.52786) q[3];
sx q[3];
rz(-0.49694967) q[3];
sx q[3];
rz(0.73408192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99877015) q[0];
sx q[0];
rz(-3.0144189) q[0];
sx q[0];
rz(2.3221827) q[0];
rz(-2.8557414) q[1];
sx q[1];
rz(-0.91393036) q[1];
sx q[1];
rz(-1.2664638) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74840141) q[0];
sx q[0];
rz(-1.1921492) q[0];
sx q[0];
rz(0.65403954) q[0];
x q[1];
rz(-2.0882294) q[2];
sx q[2];
rz(-2.6839163) q[2];
sx q[2];
rz(0.44331726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8394168) q[1];
sx q[1];
rz(-1.9315047) q[1];
sx q[1];
rz(1.9990666) q[1];
rz(-pi) q[2];
rz(0.5760848) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(-1.2036122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9849898) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(-1.9796237) q[2];
rz(-3.0544288) q[3];
sx q[3];
rz(-1.6379084) q[3];
sx q[3];
rz(3.0676837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-0.53482985) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(2.8379922) q[0];
rz(1.3820232) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(2.4940925) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8436463) q[0];
sx q[0];
rz(-0.93695153) q[0];
sx q[0];
rz(-1.6266536) q[0];
x q[1];
rz(0.71543773) q[2];
sx q[2];
rz(-1.9999265) q[2];
sx q[2];
rz(-1.5646343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.69933575) q[1];
sx q[1];
rz(-1.934821) q[1];
sx q[1];
rz(-2.060021) q[1];
rz(-pi) q[2];
rz(0.59433523) q[3];
sx q[3];
rz(-2.1777275) q[3];
sx q[3];
rz(2.2567574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8422164) q[2];
sx q[2];
rz(-2.3589578) q[2];
sx q[2];
rz(1.5807318) q[2];
rz(0.98172274) q[3];
sx q[3];
rz(-1.862062) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9075539) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(-2.2696944) q[0];
rz(0.38726989) q[1];
sx q[1];
rz(-0.64278066) q[1];
sx q[1];
rz(-2.8505468) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8431393) q[0];
sx q[0];
rz(-1.7994587) q[0];
sx q[0];
rz(-0.97897162) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.307345) q[2];
sx q[2];
rz(-2.7195889) q[2];
sx q[2];
rz(0.01854245) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2070771) q[1];
sx q[1];
rz(-1.0672489) q[1];
sx q[1];
rz(0.7657004) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99445677) q[3];
sx q[3];
rz(-1.7061702) q[3];
sx q[3];
rz(1.4801814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7724093) q[2];
sx q[2];
rz(-0.79493752) q[2];
sx q[2];
rz(-2.5975361) q[2];
rz(2.6323075) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(2.2756186) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13957025) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(2.496526) q[0];
rz(-0.6257261) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(-0.63017875) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1030578) q[0];
sx q[0];
rz(-1.6928506) q[0];
sx q[0];
rz(-0.5794258) q[0];
rz(1.1516045) q[2];
sx q[2];
rz(-0.79608166) q[2];
sx q[2];
rz(0.58005709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4634134) q[1];
sx q[1];
rz(-2.5555875) q[1];
sx q[1];
rz(-1.1581808) q[1];
x q[2];
rz(2.9003521) q[3];
sx q[3];
rz(-2.5706228) q[3];
sx q[3];
rz(0.2824479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3877635) q[2];
sx q[2];
rz(-1.8176273) q[2];
sx q[2];
rz(-1.3641664) q[2];
rz(-1.8917313) q[3];
sx q[3];
rz(-1.2156237) q[3];
sx q[3];
rz(-0.92145872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.0867778) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(-0.71682799) q[0];
rz(-2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(-2.3278918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5649367) q[0];
sx q[0];
rz(-1.7840958) q[0];
sx q[0];
rz(0.06421406) q[0];
rz(-pi) q[1];
rz(0.72197638) q[2];
sx q[2];
rz(-2.0275896) q[2];
sx q[2];
rz(1.7320088) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3354213) q[1];
sx q[1];
rz(-2.2911982) q[1];
sx q[1];
rz(2.6909188) q[1];
rz(-pi) q[2];
rz(-2.8343809) q[3];
sx q[3];
rz(-1.9091515) q[3];
sx q[3];
rz(0.070778155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8852691) q[2];
sx q[2];
rz(-0.43748125) q[2];
sx q[2];
rz(-1.1876594) q[2];
rz(2.2096283) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(-2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9880144) q[0];
sx q[0];
rz(-2.1126641) q[0];
sx q[0];
rz(-2.4555092) q[0];
rz(1.5230806) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(2.4783321) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30090573) q[0];
sx q[0];
rz(-1.5468883) q[0];
sx q[0];
rz(-0.92047128) q[0];
x q[1];
rz(-1.7887605) q[2];
sx q[2];
rz(-2.1969165) q[2];
sx q[2];
rz(1.8916212) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7546332) q[1];
sx q[1];
rz(-2.6135751) q[1];
sx q[1];
rz(-1.1623043) q[1];
x q[2];
rz(-2.0459941) q[3];
sx q[3];
rz(-1.0605269) q[3];
sx q[3];
rz(2.1197532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66701165) q[2];
sx q[2];
rz(-0.8374927) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(2.1298501) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(0.57141465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052977234) q[0];
sx q[0];
rz(-2.4089854) q[0];
sx q[0];
rz(0.10738871) q[0];
rz(-2.836851) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(-0.4577786) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25154197) q[0];
sx q[0];
rz(-0.72351447) q[0];
sx q[0];
rz(2.0680122) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0395983) q[2];
sx q[2];
rz(-2.2441412) q[2];
sx q[2];
rz(-1.7216059) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0486054) q[1];
sx q[1];
rz(-2.7106206) q[1];
sx q[1];
rz(-2.155817) q[1];
x q[2];
rz(-0.55191603) q[3];
sx q[3];
rz(-1.2593927) q[3];
sx q[3];
rz(-2.6375254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7765939) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(-1.7101074) q[2];
rz(1.4252023) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(-1.2861929) q[3];
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
rz(-0.15437056) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(-0.95348683) q[0];
rz(-1.8158688) q[1];
sx q[1];
rz(-2.6486501) q[1];
sx q[1];
rz(0.58473933) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7493233) q[0];
sx q[0];
rz(-1.682878) q[0];
sx q[0];
rz(-0.49082503) q[0];
rz(1.936475) q[2];
sx q[2];
rz(-1.056991) q[2];
sx q[2];
rz(0.50525451) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1251724) q[1];
sx q[1];
rz(-2.6965045) q[1];
sx q[1];
rz(-2.9801286) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3034079) q[3];
sx q[3];
rz(-1.7087473) q[3];
sx q[3];
rz(-1.7351868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8087625) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(2.8653223) q[2];
rz(2.590498) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(2.7579257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(1.0693414) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(1.9845225) q[1];
sx q[1];
rz(-1.3815222) q[1];
sx q[1];
rz(-1.6190593) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9690543) q[0];
sx q[0];
rz(-1.2612871) q[0];
sx q[0];
rz(0.68185101) q[0];
rz(-pi) q[1];
rz(2.4535975) q[2];
sx q[2];
rz(-1.3032459) q[2];
sx q[2];
rz(-0.5984532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9138227) q[1];
sx q[1];
rz(-1.4064944) q[1];
sx q[1];
rz(-0.63197244) q[1];
rz(-1.9736245) q[3];
sx q[3];
rz(-1.1790923) q[3];
sx q[3];
rz(-0.21564461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.15554252) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(-0.84214169) q[2];
rz(1.5367674) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(-3.1150637) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90606541) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(2.4189667) q[1];
sx q[1];
rz(-2.2650748) q[1];
sx q[1];
rz(1.5320019) q[1];
rz(1.6956383) q[2];
sx q[2];
rz(-2.9308133) q[2];
sx q[2];
rz(1.4598893) q[2];
rz(-1.6668672) q[3];
sx q[3];
rz(-0.69312743) q[3];
sx q[3];
rz(3.0467924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];