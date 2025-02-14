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
rz(-1.4827363) q[0];
sx q[0];
rz(-2.1602477) q[0];
sx q[0];
rz(2.0438097) q[0];
rz(2.3535347) q[1];
sx q[1];
rz(-1.0507974) q[1];
sx q[1];
rz(0.0069590574) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.033125) q[0];
sx q[0];
rz(-1.5146717) q[0];
sx q[0];
rz(-2.2214487) q[0];
rz(-pi) q[1];
rz(0.2351951) q[2];
sx q[2];
rz(-0.93738745) q[2];
sx q[2];
rz(-0.051284479) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6485868) q[1];
sx q[1];
rz(-1.3351591) q[1];
sx q[1];
rz(-3.1070903) q[1];
x q[2];
rz(0.29421774) q[3];
sx q[3];
rz(-2.7876884) q[3];
sx q[3];
rz(2.317441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9030582) q[2];
sx q[2];
rz(-1.7944585) q[2];
sx q[2];
rz(1.6713589) q[2];
rz(-3.0155731) q[3];
sx q[3];
rz(-0.44243789) q[3];
sx q[3];
rz(2.457705) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0266492) q[0];
sx q[0];
rz(-2.6925955) q[0];
sx q[0];
rz(-2.5057416) q[0];
rz(-2.3682829) q[1];
sx q[1];
rz(-2.4644303) q[1];
sx q[1];
rz(1.6962475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.792004) q[0];
sx q[0];
rz(-2.3161962) q[0];
sx q[0];
rz(0.94657739) q[0];
rz(-pi) q[1];
rz(2.6670806) q[2];
sx q[2];
rz(-0.84174985) q[2];
sx q[2];
rz(3.0345033) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.066582908) q[1];
sx q[1];
rz(-2.420418) q[1];
sx q[1];
rz(-2.9020792) q[1];
rz(-1.4392733) q[3];
sx q[3];
rz(-1.5029385) q[3];
sx q[3];
rz(2.0236286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9952952) q[2];
sx q[2];
rz(-1.7701021) q[2];
sx q[2];
rz(3.0038317) q[2];
rz(-0.45608258) q[3];
sx q[3];
rz(-2.5465953) q[3];
sx q[3];
rz(0.47541398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59738961) q[0];
sx q[0];
rz(-1.5837639) q[0];
sx q[0];
rz(-2.5624045) q[0];
rz(-0.10313615) q[1];
sx q[1];
rz(-1.3214279) q[1];
sx q[1];
rz(-0.78027049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1121169) q[0];
sx q[0];
rz(-1.5519841) q[0];
sx q[0];
rz(-3.0464059) q[0];
rz(-pi) q[1];
rz(-0.19635503) q[2];
sx q[2];
rz(-0.69075023) q[2];
sx q[2];
rz(-3.1121814) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5580328) q[1];
sx q[1];
rz(-1.7906902) q[1];
sx q[1];
rz(-1.1660485) q[1];
rz(-pi) q[2];
rz(-1.8336201) q[3];
sx q[3];
rz(-2.4287831) q[3];
sx q[3];
rz(1.7607016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4812193) q[2];
sx q[2];
rz(-1.6796835) q[2];
sx q[2];
rz(3.050728) q[2];
rz(-1.4800492) q[3];
sx q[3];
rz(-1.3250947) q[3];
sx q[3];
rz(2.3265694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.018983) q[0];
sx q[0];
rz(-2.9538587) q[0];
sx q[0];
rz(-0.51938272) q[0];
rz(-1.1795562) q[1];
sx q[1];
rz(-2.4603381) q[1];
sx q[1];
rz(-3.093241) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3607442) q[0];
sx q[0];
rz(-2.0791302) q[0];
sx q[0];
rz(-0.8698277) q[0];
x q[1];
rz(-0.91421336) q[2];
sx q[2];
rz(-2.5541452) q[2];
sx q[2];
rz(2.812831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0569339) q[1];
sx q[1];
rz(-0.48843436) q[1];
sx q[1];
rz(-0.059367511) q[1];
x q[2];
rz(1.3062258) q[3];
sx q[3];
rz(-2.3481784) q[3];
sx q[3];
rz(-2.1891264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68253303) q[2];
sx q[2];
rz(-0.30632633) q[2];
sx q[2];
rz(1.691386) q[2];
rz(-1.1059443) q[3];
sx q[3];
rz(-1.148843) q[3];
sx q[3];
rz(2.2753184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3806216) q[0];
sx q[0];
rz(-2.3411317) q[0];
sx q[0];
rz(-0.77734429) q[0];
rz(-0.49579534) q[1];
sx q[1];
rz(-0.40758857) q[1];
sx q[1];
rz(-0.19283238) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042943311) q[0];
sx q[0];
rz(-1.0355562) q[0];
sx q[0];
rz(-1.9283867) q[0];
rz(-pi) q[1];
rz(0.2482581) q[2];
sx q[2];
rz(-2.358486) q[2];
sx q[2];
rz(-0.017489028) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1950394) q[1];
sx q[1];
rz(-1.4832388) q[1];
sx q[1];
rz(-1.273543) q[1];
x q[2];
rz(0.16280547) q[3];
sx q[3];
rz(-1.1737524) q[3];
sx q[3];
rz(-0.59418488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5567646) q[2];
sx q[2];
rz(-0.77797055) q[2];
sx q[2];
rz(2.3053816) q[2];
rz(2.1861475) q[3];
sx q[3];
rz(-1.4232057) q[3];
sx q[3];
rz(2.6098765) q[3];
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
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8504976) q[0];
sx q[0];
rz(-1.9170772) q[0];
sx q[0];
rz(-3.1078872) q[0];
rz(-2.2233502) q[1];
sx q[1];
rz(-0.70437175) q[1];
sx q[1];
rz(-2.4423626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7875154) q[0];
sx q[0];
rz(-1.0654385) q[0];
sx q[0];
rz(-0.92886852) q[0];
rz(2.8644833) q[2];
sx q[2];
rz(-0.19453262) q[2];
sx q[2];
rz(2.3152318) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.67656006) q[1];
sx q[1];
rz(-1.1313712) q[1];
sx q[1];
rz(1.7656754) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27856234) q[3];
sx q[3];
rz(-1.8980366) q[3];
sx q[3];
rz(2.183941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0596727) q[2];
sx q[2];
rz(-2.4884188) q[2];
sx q[2];
rz(-0.095452249) q[2];
rz(1.0517612) q[3];
sx q[3];
rz(-1.8973408) q[3];
sx q[3];
rz(0.89422798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28602257) q[0];
sx q[0];
rz(-0.48208553) q[0];
sx q[0];
rz(-1.9739738) q[0];
rz(2.7883912) q[1];
sx q[1];
rz(-0.51274061) q[1];
sx q[1];
rz(2.8599427) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7324556) q[0];
sx q[0];
rz(-1.0971709) q[0];
sx q[0];
rz(3.043463) q[0];
x q[1];
rz(-0.64729624) q[2];
sx q[2];
rz(-1.656083) q[2];
sx q[2];
rz(2.6814987) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0822009) q[1];
sx q[1];
rz(-1.4363465) q[1];
sx q[1];
rz(-1.4217292) q[1];
x q[2];
rz(1.0135039) q[3];
sx q[3];
rz(-1.1508599) q[3];
sx q[3];
rz(0.98972631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1436651) q[2];
sx q[2];
rz(-1.9196271) q[2];
sx q[2];
rz(-1.3227051) q[2];
rz(-1.6392684) q[3];
sx q[3];
rz(-1.084525) q[3];
sx q[3];
rz(3.0433906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99763501) q[0];
sx q[0];
rz(-1.8956381) q[0];
sx q[0];
rz(2.8676046) q[0];
rz(0.59448376) q[1];
sx q[1];
rz(-1.4130054) q[1];
sx q[1];
rz(0.53057539) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6903821) q[0];
sx q[0];
rz(-1.57461) q[0];
sx q[0];
rz(3.1022152) q[0];
x q[1];
rz(0.12315788) q[2];
sx q[2];
rz(-2.5854857) q[2];
sx q[2];
rz(-1.7983939) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.12996) q[1];
sx q[1];
rz(-1.7372565) q[1];
sx q[1];
rz(0.37977207) q[1];
rz(-pi) q[2];
rz(0.16333111) q[3];
sx q[3];
rz(-2.2769351) q[3];
sx q[3];
rz(2.284019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24255594) q[2];
sx q[2];
rz(-1.9344923) q[2];
sx q[2];
rz(-2.8846018) q[2];
rz(1.1993923) q[3];
sx q[3];
rz(-1.896984) q[3];
sx q[3];
rz(-1.6283584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5523819) q[0];
sx q[0];
rz(-2.1871545) q[0];
sx q[0];
rz(0.5300262) q[0];
rz(1.5026622) q[1];
sx q[1];
rz(-1.9866147) q[1];
sx q[1];
rz(-0.93856215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8535271) q[0];
sx q[0];
rz(-2.1479283) q[0];
sx q[0];
rz(1.0017299) q[0];
x q[1];
rz(2.9979894) q[2];
sx q[2];
rz(-1.4521027) q[2];
sx q[2];
rz(2.375756) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0238513) q[1];
sx q[1];
rz(-1.3582025) q[1];
sx q[1];
rz(2.0076058) q[1];
rz(0.13646941) q[3];
sx q[3];
rz(-0.9896864) q[3];
sx q[3];
rz(-0.438353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3878801) q[2];
sx q[2];
rz(-2.2215999) q[2];
sx q[2];
rz(2.06125) q[2];
rz(-0.8148109) q[3];
sx q[3];
rz(-1.1956513) q[3];
sx q[3];
rz(1.4174392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1431047) q[0];
sx q[0];
rz(-1.657722) q[0];
sx q[0];
rz(-1.0888354) q[0];
rz(3.0740652) q[1];
sx q[1];
rz(-1.2779526) q[1];
sx q[1];
rz(0.75928226) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1812173) q[0];
sx q[0];
rz(-1.8255673) q[0];
sx q[0];
rz(-3.0187458) q[0];
rz(-pi) q[1];
rz(0.22589182) q[2];
sx q[2];
rz(-1.6917657) q[2];
sx q[2];
rz(-2.8793983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2226505) q[1];
sx q[1];
rz(-0.83178751) q[1];
sx q[1];
rz(1.8636501) q[1];
rz(-pi) q[2];
rz(-1.8005695) q[3];
sx q[3];
rz(-1.2342576) q[3];
sx q[3];
rz(1.1024818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8074983) q[2];
sx q[2];
rz(-1.0958025) q[2];
sx q[2];
rz(-2.5661772) q[2];
rz(0.56257644) q[3];
sx q[3];
rz(-2.176599) q[3];
sx q[3];
rz(1.7687198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0811049) q[0];
sx q[0];
rz(-1.2860379) q[0];
sx q[0];
rz(2.2338569) q[0];
rz(2.0545215) q[1];
sx q[1];
rz(-2.0094951) q[1];
sx q[1];
rz(3.1241945) q[1];
rz(-0.74093735) q[2];
sx q[2];
rz(-0.20812427) q[2];
sx q[2];
rz(-0.32231449) q[2];
rz(2.5431332) q[3];
sx q[3];
rz(-1.2597407) q[3];
sx q[3];
rz(1.5516439) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
