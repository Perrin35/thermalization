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
rz(0.61007845) q[0];
sx q[0];
rz(3.371513) q[0];
sx q[0];
rz(11.846677) q[0];
rz(2.5201058) q[1];
sx q[1];
rz(-1.7401594) q[1];
sx q[1];
rz(-3.016234) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014356456) q[0];
sx q[0];
rz(-1.171021) q[0];
sx q[0];
rz(-0.46443224) q[0];
rz(-pi) q[1];
rz(-1.4761476) q[2];
sx q[2];
rz(-1.2834335) q[2];
sx q[2];
rz(-2.8146625) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6290063) q[1];
sx q[1];
rz(-0.0014481469) q[1];
sx q[1];
rz(1.8667029) q[1];
rz(-pi) q[2];
rz(-0.8475581) q[3];
sx q[3];
rz(-1.8430018) q[3];
sx q[3];
rz(-2.6295626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7980935) q[2];
sx q[2];
rz(-0.40830475) q[2];
sx q[2];
rz(-2.2957392) q[2];
rz(-2.3503303) q[3];
sx q[3];
rz(-0.013412272) q[3];
sx q[3];
rz(-3.0748034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.7283519) q[0];
sx q[0];
rz(-0.46717307) q[0];
sx q[0];
rz(-3.0979284) q[0];
rz(1.5751669) q[1];
sx q[1];
rz(-1.3730201) q[1];
sx q[1];
rz(1.498819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18066809) q[0];
sx q[0];
rz(-2.5378413) q[0];
sx q[0];
rz(-1.5636428) q[0];
rz(1.5909219) q[2];
sx q[2];
rz(-0.57148904) q[2];
sx q[2];
rz(0.0062985346) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0017569204) q[1];
sx q[1];
rz(-0.082232177) q[1];
sx q[1];
rz(-2.7539192) q[1];
rz(-pi) q[2];
rz(3.1083634) q[3];
sx q[3];
rz(-1.2265612) q[3];
sx q[3];
rz(1.1297463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0955536) q[2];
sx q[2];
rz(-2.9914896) q[2];
sx q[2];
rz(-0.52774876) q[2];
rz(2.2857417) q[3];
sx q[3];
rz(-3.1400883) q[3];
sx q[3];
rz(1.9201479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7432778) q[0];
sx q[0];
rz(-2.1687431) q[0];
sx q[0];
rz(-1.1463746) q[0];
rz(-1.758681) q[1];
sx q[1];
rz(-2.8490366) q[1];
sx q[1];
rz(0.10398277) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87097528) q[0];
sx q[0];
rz(-1.7954602) q[0];
sx q[0];
rz(-1.490834) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2984593) q[2];
sx q[2];
rz(-1.587279) q[2];
sx q[2];
rz(-2.4695244) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1749035) q[1];
sx q[1];
rz(-0.82920584) q[1];
sx q[1];
rz(1.8494383) q[1];
rz(-pi) q[2];
rz(0.35211925) q[3];
sx q[3];
rz(-0.30856347) q[3];
sx q[3];
rz(-1.4328014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.18342239) q[2];
sx q[2];
rz(-0.0068181097) q[2];
sx q[2];
rz(0.56162322) q[2];
rz(0.084847458) q[3];
sx q[3];
rz(-0.0054797879) q[3];
sx q[3];
rz(3.1408299) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7247923) q[0];
sx q[0];
rz(-0.264072) q[0];
sx q[0];
rz(2.6208139) q[0];
rz(-2.9837823) q[1];
sx q[1];
rz(-2.4746555) q[1];
sx q[1];
rz(3.0657943) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6822081) q[0];
sx q[0];
rz(-2.5494116) q[0];
sx q[0];
rz(-0.20488744) q[0];
rz(-pi) q[1];
rz(1.5694322) q[2];
sx q[2];
rz(-1.5701541) q[2];
sx q[2];
rz(-3.0027892) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0204649) q[1];
sx q[1];
rz(-1.1942099) q[1];
sx q[1];
rz(2.9894606) q[1];
rz(-1.789572) q[3];
sx q[3];
rz(-0.46275381) q[3];
sx q[3];
rz(-2.3298711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6303404) q[2];
sx q[2];
rz(-3.1255836) q[2];
sx q[2];
rz(-1.7208257) q[2];
rz(-3.1347647) q[3];
sx q[3];
rz(-3.112401) q[3];
sx q[3];
rz(1.6543057) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1348949) q[0];
sx q[0];
rz(-0.20773523) q[0];
sx q[0];
rz(2.8790706) q[0];
rz(2.1916356) q[1];
sx q[1];
rz(-0.078153178) q[1];
sx q[1];
rz(-2.9238759) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11853795) q[0];
sx q[0];
rz(-0.41363198) q[0];
sx q[0];
rz(1.8819811) q[0];
x q[1];
rz(3.0543572) q[2];
sx q[2];
rz(-1.4873184) q[2];
sx q[2];
rz(-1.6501476) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7694665) q[1];
sx q[1];
rz(-1.7431152) q[1];
sx q[1];
rz(-0.0021688633) q[1];
x q[2];
rz(1.6895377) q[3];
sx q[3];
rz(-2.6104527) q[3];
sx q[3];
rz(1.0815999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.21506423) q[2];
sx q[2];
rz(-0.0096409163) q[2];
sx q[2];
rz(2.8936774) q[2];
rz(-1.7667814) q[3];
sx q[3];
rz(-0.050516613) q[3];
sx q[3];
rz(-1.4352528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0404102) q[0];
sx q[0];
rz(-1.8721767) q[0];
sx q[0];
rz(2.5293479) q[0];
rz(-0.170389) q[1];
sx q[1];
rz(-0.082516106) q[1];
sx q[1];
rz(-1.4917397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1338324) q[0];
sx q[0];
rz(-2.2585224) q[0];
sx q[0];
rz(2.4439815) q[0];
x q[1];
rz(-1.5749919) q[2];
sx q[2];
rz(-1.5853264) q[2];
sx q[2];
rz(0.57455237) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9656758) q[1];
sx q[1];
rz(-1.2781483) q[1];
sx q[1];
rz(-1.501237) q[1];
rz(-pi) q[2];
rz(0.18317353) q[3];
sx q[3];
rz(-1.8889931) q[3];
sx q[3];
rz(-3.0903926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6461569) q[2];
sx q[2];
rz(-0.010224552) q[2];
sx q[2];
rz(2.9025027) q[2];
rz(-2.332989) q[3];
sx q[3];
rz(-3.1294398) q[3];
sx q[3];
rz(1.808572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072153) q[0];
sx q[0];
rz(-3.0476397) q[0];
sx q[0];
rz(0.77675003) q[0];
rz(-0.010919318) q[1];
sx q[1];
rz(-0.24591406) q[1];
sx q[1];
rz(1.4728665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0162139) q[0];
sx q[0];
rz(-0.74170602) q[0];
sx q[0];
rz(1.2522231) q[0];
x q[1];
rz(-2.8773715) q[2];
sx q[2];
rz(-3.1172385) q[2];
sx q[2];
rz(3.0517677) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57629055) q[1];
sx q[1];
rz(-1.6716752) q[1];
sx q[1];
rz(2.3957344) q[1];
rz(-0.075447791) q[3];
sx q[3];
rz(-1.8561072) q[3];
sx q[3];
rz(-2.4782654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9933068) q[2];
sx q[2];
rz(-0.017227087) q[2];
sx q[2];
rz(-2.7255608) q[2];
rz(-2.3038583) q[3];
sx q[3];
rz(-0.13844027) q[3];
sx q[3];
rz(1.6637038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96529043) q[0];
sx q[0];
rz(-0.039881341) q[0];
sx q[0];
rz(2.1862929) q[0];
rz(-1.924986) q[1];
sx q[1];
rz(-2.80426) q[1];
sx q[1];
rz(0.40562707) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8389616) q[0];
sx q[0];
rz(-1.8895188) q[0];
sx q[0];
rz(-0.74619729) q[0];
x q[1];
rz(2.9826016) q[2];
sx q[2];
rz(-1.6588677) q[2];
sx q[2];
rz(-0.51364567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15524023) q[1];
sx q[1];
rz(-1.6472247) q[1];
sx q[1];
rz(-0.16096154) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38281103) q[3];
sx q[3];
rz(-2.9124526) q[3];
sx q[3];
rz(2.7349796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2847298) q[2];
sx q[2];
rz(-3.1046107) q[2];
sx q[2];
rz(-2.3178318) q[2];
rz(-3.01037) q[3];
sx q[3];
rz(-2.8516912) q[3];
sx q[3];
rz(0.32740617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71309483) q[0];
sx q[0];
rz(-2.9976124) q[0];
sx q[0];
rz(2.4289828) q[0];
rz(0.54829848) q[1];
sx q[1];
rz(-0.24990853) q[1];
sx q[1];
rz(2.9329494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31748235) q[0];
sx q[0];
rz(-1.1718996) q[0];
sx q[0];
rz(-2.0092416) q[0];
rz(-pi) q[1];
rz(-2.0698572) q[2];
sx q[2];
rz(-0.03943561) q[2];
sx q[2];
rz(-2.5372504) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5447158) q[1];
sx q[1];
rz(-1.5660389) q[1];
sx q[1];
rz(2.1265043) q[1];
rz(-pi) q[2];
rz(-2.194368) q[3];
sx q[3];
rz(-1.0290909) q[3];
sx q[3];
rz(0.052963363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6930406) q[2];
sx q[2];
rz(-0.081144944) q[2];
sx q[2];
rz(-0.21716675) q[2];
rz(-0.36755696) q[3];
sx q[3];
rz(-3.1072072) q[3];
sx q[3];
rz(0.99678451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17598584) q[0];
sx q[0];
rz(-3.0526057) q[0];
sx q[0];
rz(0.17372818) q[0];
rz(-1.4088176) q[1];
sx q[1];
rz(-1.6830187) q[1];
sx q[1];
rz(-1.6857612) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.78409) q[0];
sx q[0];
rz(-1.0561868) q[0];
sx q[0];
rz(1.3823439) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7782434) q[2];
sx q[2];
rz(-1.7073627) q[2];
sx q[2];
rz(-2.2547124) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.10911726) q[1];
sx q[1];
rz(-2.2665503) q[1];
sx q[1];
rz(1.4182219) q[1];
x q[2];
rz(1.3221856) q[3];
sx q[3];
rz(-1.4668307) q[3];
sx q[3];
rz(-1.6373529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9659861) q[2];
sx q[2];
rz(-0.49314988) q[2];
sx q[2];
rz(-1.3803587) q[2];
rz(0.68584758) q[3];
sx q[3];
rz(-0.0018456056) q[3];
sx q[3];
rz(0.67870158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7420409) q[0];
sx q[0];
rz(-0.69235943) q[0];
sx q[0];
rz(1.6211744) q[0];
rz(-1.5581268) q[1];
sx q[1];
rz(-1.635066) q[1];
sx q[1];
rz(0.20546694) q[1];
rz(0.016318446) q[2];
sx q[2];
rz(-1.6961799) q[2];
sx q[2];
rz(0.21519306) q[2];
rz(-0.11859433) q[3];
sx q[3];
rz(-1.2731324) q[3];
sx q[3];
rz(-0.29407339) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
