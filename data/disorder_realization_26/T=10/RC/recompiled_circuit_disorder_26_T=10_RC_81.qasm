OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(4.8298782) q[0];
sx q[0];
rz(9.7363135) q[0];
rz(-0.43752924) q[1];
sx q[1];
rz(-1.8234) q[1];
sx q[1];
rz(0.55895609) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5200978) q[0];
sx q[0];
rz(-1.1557475) q[0];
sx q[0];
rz(0.15226224) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5842701) q[2];
sx q[2];
rz(-1.5601336) q[2];
sx q[2];
rz(2.9023841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47145876) q[1];
sx q[1];
rz(-2.11073) q[1];
sx q[1];
rz(0.88175168) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5823703) q[3];
sx q[3];
rz(-2.6359574) q[3];
sx q[3];
rz(1.3666183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7493593) q[2];
sx q[2];
rz(-1.2831251) q[2];
sx q[2];
rz(-0.63670811) q[2];
rz(-0.84896815) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(-2.9076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7648776) q[0];
sx q[0];
rz(-0.24704084) q[0];
sx q[0];
rz(2.9887181) q[0];
rz(0.75694594) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(0.98639948) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5312885) q[0];
sx q[0];
rz(-3.0695519) q[0];
sx q[0];
rz(2.5487367) q[0];
rz(-0.7631626) q[2];
sx q[2];
rz(-1.9689416) q[2];
sx q[2];
rz(-2.5415908) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4437372) q[1];
sx q[1];
rz(-1.9471696) q[1];
sx q[1];
rz(1.8259551) q[1];
x q[2];
rz(-1.8996703) q[3];
sx q[3];
rz(-1.4423014) q[3];
sx q[3];
rz(2.0957029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5793005) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(-0.78312773) q[2];
rz(0.018571818) q[3];
sx q[3];
rz(-1.5037856) q[3];
sx q[3];
rz(0.40772453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0531533) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(0.96167481) q[0];
rz(-0.36034521) q[1];
sx q[1];
rz(-2.0397489) q[1];
sx q[1];
rz(-0.12869421) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51552396) q[0];
sx q[0];
rz(-1.4883853) q[0];
sx q[0];
rz(-3.0936196) q[0];
x q[1];
rz(2.6016597) q[2];
sx q[2];
rz(-0.88368249) q[2];
sx q[2];
rz(-0.5772669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.62540185) q[1];
sx q[1];
rz(-1.1266202) q[1];
sx q[1];
rz(1.2605915) q[1];
rz(-pi) q[2];
rz(0.80273654) q[3];
sx q[3];
rz(-1.4736796) q[3];
sx q[3];
rz(3.0320398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.06015691) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(-1.241768) q[2];
rz(2.5545819) q[3];
sx q[3];
rz(-2.185052) q[3];
sx q[3];
rz(0.97755066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63242763) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(-2.0571016) q[0];
rz(-1.658461) q[1];
sx q[1];
rz(-2.5741534) q[1];
sx q[1];
rz(-3.049057) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0412484) q[0];
sx q[0];
rz(-2.7462602) q[0];
sx q[0];
rz(0.80907099) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1763457) q[2];
sx q[2];
rz(-0.56376981) q[2];
sx q[2];
rz(-3.0639067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1558518) q[1];
sx q[1];
rz(-0.34793138) q[1];
sx q[1];
rz(-2.9702529) q[1];
rz(-pi) q[2];
rz(2.2257005) q[3];
sx q[3];
rz(-1.6989143) q[3];
sx q[3];
rz(1.7803943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3080421) q[2];
sx q[2];
rz(-1.7001067) q[2];
sx q[2];
rz(0.33205024) q[2];
rz(-1.0559233) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(-0.61029303) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32245359) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(-0.21155587) q[0];
rz(1.3062723) q[1];
sx q[1];
rz(-1.897656) q[1];
sx q[1];
rz(0.64770118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0603795) q[0];
sx q[0];
rz(-0.19506422) q[0];
sx q[0];
rz(2.5293406) q[0];
rz(-0.34824246) q[2];
sx q[2];
rz(-1.7160545) q[2];
sx q[2];
rz(0.20656221) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60407818) q[1];
sx q[1];
rz(-1.068183) q[1];
sx q[1];
rz(-2.9167049) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1224498) q[3];
sx q[3];
rz(-1.585841) q[3];
sx q[3];
rz(-0.5152094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5806879) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(0.69331759) q[2];
rz(-2.4723315) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82419056) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(-1.2325226) q[0];
rz(-1.0725853) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(-0.17428621) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0347736) q[0];
sx q[0];
rz(-0.79402059) q[0];
sx q[0];
rz(1.130571) q[0];
rz(-pi) q[1];
rz(2.2767378) q[2];
sx q[2];
rz(-1.0360498) q[2];
sx q[2];
rz(2.6851482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6433405) q[1];
sx q[1];
rz(-2.2689515) q[1];
sx q[1];
rz(0.50486418) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11630451) q[3];
sx q[3];
rz(-1.260584) q[3];
sx q[3];
rz(2.5935964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8217414) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(0.19763395) q[2];
rz(-0.28891426) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(1.4060085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0777247) q[0];
sx q[0];
rz(-2.6517695) q[0];
sx q[0];
rz(2.9329964) q[0];
rz(-2.1754307) q[1];
sx q[1];
rz(-1.1599133) q[1];
sx q[1];
rz(-1.6360412) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9150328) q[0];
sx q[0];
rz(-1.0462927) q[0];
sx q[0];
rz(-0.71136186) q[0];
x q[1];
rz(0.43039544) q[2];
sx q[2];
rz(-1.5376523) q[2];
sx q[2];
rz(0.37552777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58563102) q[1];
sx q[1];
rz(-2.0913843) q[1];
sx q[1];
rz(0.52557892) q[1];
rz(-0.3397293) q[3];
sx q[3];
rz(-2.8390084) q[3];
sx q[3];
rz(1.8031977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2839526) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(-0.097578438) q[2];
rz(-1.7476667) q[3];
sx q[3];
rz(-1.7912309) q[3];
sx q[3];
rz(2.424749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39032787) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(-0.51399291) q[0];
rz(0.12318525) q[1];
sx q[1];
rz(-0.24736483) q[1];
sx q[1];
rz(-2.2095912) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65864627) q[0];
sx q[0];
rz(-2.6771149) q[0];
sx q[0];
rz(1.8770201) q[0];
rz(-0.93512647) q[2];
sx q[2];
rz(-1.3526275) q[2];
sx q[2];
rz(-0.93014923) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8646647) q[1];
sx q[1];
rz(-2.118646) q[1];
sx q[1];
rz(2.1879556) q[1];
rz(-pi) q[2];
rz(2.9429432) q[3];
sx q[3];
rz(-0.86411398) q[3];
sx q[3];
rz(-2.228565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0294068) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(-0.4294447) q[2];
rz(1.2094234) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(2.4485574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747303) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(-1.8027579) q[0];
rz(2.4354637) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(0.95058092) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4336193) q[0];
sx q[0];
rz(-2.3971359) q[0];
sx q[0];
rz(2.3696193) q[0];
rz(3.0663475) q[2];
sx q[2];
rz(-1.9855472) q[2];
sx q[2];
rz(0.097397734) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49950019) q[1];
sx q[1];
rz(-1.7576808) q[1];
sx q[1];
rz(-2.1128113) q[1];
x q[2];
rz(1.1589963) q[3];
sx q[3];
rz(-1.8048865) q[3];
sx q[3];
rz(1.7115418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6616228) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(-0.92710036) q[3];
sx q[3];
rz(-1.3092594) q[3];
sx q[3];
rz(2.5203729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9973688) q[0];
sx q[0];
rz(-3.0529418) q[0];
sx q[0];
rz(2.9123059) q[0];
rz(0.43481049) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(2.4226709) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8573498) q[0];
sx q[0];
rz(-1.4765258) q[0];
sx q[0];
rz(-2.1988792) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.747379) q[2];
sx q[2];
rz(-0.59497661) q[2];
sx q[2];
rz(2.8240311) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1246008) q[1];
sx q[1];
rz(-1.4233839) q[1];
sx q[1];
rz(0.37313811) q[1];
rz(-2.1925681) q[3];
sx q[3];
rz(-1.9927295) q[3];
sx q[3];
rz(-2.4728647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3832613) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(-0.74679217) q[2];
rz(-2.2693999) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(2.1993568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.223021) q[0];
sx q[0];
rz(-1.7445607) q[0];
sx q[0];
rz(1.8013409) q[0];
rz(0.37721286) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(2.626426) q[2];
sx q[2];
rz(-0.58044051) q[2];
sx q[2];
rz(2.6945111) q[2];
rz(-2.4242998) q[3];
sx q[3];
rz(-1.4581231) q[3];
sx q[3];
rz(0.067618528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];