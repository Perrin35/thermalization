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
rz(0.084963381) q[0];
sx q[0];
rz(-2.8391916) q[0];
sx q[0];
rz(3.1095355) q[0];
rz(-0.0078460296) q[1];
sx q[1];
rz(-0.50385952) q[1];
sx q[1];
rz(2.388968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7330878) q[0];
sx q[0];
rz(-1.2512387) q[0];
sx q[0];
rz(-1.7825141) q[0];
rz(-2.7787894) q[2];
sx q[2];
rz(-1.9224836) q[2];
sx q[2];
rz(1.2376518) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6258273) q[1];
sx q[1];
rz(-1.2316362) q[1];
sx q[1];
rz(-1.304342) q[1];
x q[2];
rz(2.9273622) q[3];
sx q[3];
rz(-1.2204683) q[3];
sx q[3];
rz(0.080834724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1746615) q[2];
sx q[2];
rz(-1.9661247) q[2];
sx q[2];
rz(-0.82007972) q[2];
rz(2.7005633) q[3];
sx q[3];
rz(-1.1038154) q[3];
sx q[3];
rz(2.5681514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012861982) q[0];
sx q[0];
rz(-2.2282889) q[0];
sx q[0];
rz(0.41020694) q[0];
rz(0.38093105) q[1];
sx q[1];
rz(-1.5313287) q[1];
sx q[1];
rz(1.7832696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84285865) q[0];
sx q[0];
rz(-2.4141387) q[0];
sx q[0];
rz(2.1841315) q[0];
rz(2.1141421) q[2];
sx q[2];
rz(-2.0110235) q[2];
sx q[2];
rz(-1.5433951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4315003) q[1];
sx q[1];
rz(-2.9572801) q[1];
sx q[1];
rz(-1.1460365) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7776599) q[3];
sx q[3];
rz(-0.65012041) q[3];
sx q[3];
rz(0.93322414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35839781) q[2];
sx q[2];
rz(-2.7326549) q[2];
sx q[2];
rz(-0.50296339) q[2];
rz(1.2991692) q[3];
sx q[3];
rz(-0.75853577) q[3];
sx q[3];
rz(-2.9097596) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4101039) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(2.1657535) q[0];
rz(-2.2052235) q[1];
sx q[1];
rz(-1.5872637) q[1];
sx q[1];
rz(-0.13793129) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9693404) q[0];
sx q[0];
rz(-1.9346049) q[0];
sx q[0];
rz(1.0602289) q[0];
x q[1];
rz(-0.39769002) q[2];
sx q[2];
rz(-1.4639336) q[2];
sx q[2];
rz(1.7025089) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8616069) q[1];
sx q[1];
rz(-2.2642038) q[1];
sx q[1];
rz(-2.9704291) q[1];
rz(-1.8432328) q[3];
sx q[3];
rz(-1.1198634) q[3];
sx q[3];
rz(0.48905269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58933538) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(-1.4790347) q[2];
rz(-0.79834437) q[3];
sx q[3];
rz(-0.84404293) q[3];
sx q[3];
rz(-3.121283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.289157) q[0];
sx q[0];
rz(-0.11141369) q[0];
sx q[0];
rz(1.0668466) q[0];
rz(1.5486708) q[1];
sx q[1];
rz(-1.2499864) q[1];
sx q[1];
rz(0.41935316) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3187554) q[0];
sx q[0];
rz(-1.4117804) q[0];
sx q[0];
rz(-1.8029638) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9871171) q[2];
sx q[2];
rz(-0.53336582) q[2];
sx q[2];
rz(-1.6090924) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6450582) q[1];
sx q[1];
rz(-2.7784111) q[1];
sx q[1];
rz(1.7523604) q[1];
x q[2];
rz(-2.7953496) q[3];
sx q[3];
rz(-1.1107363) q[3];
sx q[3];
rz(-1.3390883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7833917) q[2];
sx q[2];
rz(-0.50938598) q[2];
sx q[2];
rz(1.5436714) q[2];
rz(1.2265497) q[3];
sx q[3];
rz(-1.9251325) q[3];
sx q[3];
rz(1.290087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.876038) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(0.75620404) q[0];
rz(-1.8887695) q[1];
sx q[1];
rz(-2.2105261) q[1];
sx q[1];
rz(1.7255712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21391103) q[0];
sx q[0];
rz(-2.1577303) q[0];
sx q[0];
rz(-1.3089433) q[0];
rz(-pi) q[1];
rz(-2.8019287) q[2];
sx q[2];
rz(-0.36298266) q[2];
sx q[2];
rz(1.5659005) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.071347) q[1];
sx q[1];
rz(-1.5778549) q[1];
sx q[1];
rz(0.53877212) q[1];
rz(1.3802212) q[3];
sx q[3];
rz(-1.9280199) q[3];
sx q[3];
rz(-0.009454184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34117928) q[2];
sx q[2];
rz(-0.74113733) q[2];
sx q[2];
rz(0.18079147) q[2];
rz(1.6527269) q[3];
sx q[3];
rz(-1.7427665) q[3];
sx q[3];
rz(1.5159877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7086696) q[0];
sx q[0];
rz(-2.0661418) q[0];
sx q[0];
rz(-2.3413626) q[0];
rz(0.284614) q[1];
sx q[1];
rz(-1.0601284) q[1];
sx q[1];
rz(-0.12399331) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1127045) q[0];
sx q[0];
rz(-1.7005973) q[0];
sx q[0];
rz(0.028667269) q[0];
x q[1];
rz(-1.8058067) q[2];
sx q[2];
rz(-1.5127276) q[2];
sx q[2];
rz(-0.11561671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1734487) q[1];
sx q[1];
rz(-0.3160797) q[1];
sx q[1];
rz(2.8760002) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1670951) q[3];
sx q[3];
rz(-1.8829405) q[3];
sx q[3];
rz(2.6016935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93512145) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(3.0653811) q[2];
rz(1.291409) q[3];
sx q[3];
rz(-2.0468057) q[3];
sx q[3];
rz(-2.5230303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718219) q[0];
sx q[0];
rz(-1.5795647) q[0];
sx q[0];
rz(2.3845657) q[0];
rz(-0.79611671) q[1];
sx q[1];
rz(-1.7346104) q[1];
sx q[1];
rz(0.98181358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2404838) q[0];
sx q[0];
rz(-0.79190688) q[0];
sx q[0];
rz(1.1674985) q[0];
rz(-pi) q[1];
rz(-0.070008833) q[2];
sx q[2];
rz(-1.677779) q[2];
sx q[2];
rz(-1.2194404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.422324) q[1];
sx q[1];
rz(-1.1012205) q[1];
sx q[1];
rz(2.6828241) q[1];
rz(1.8514093) q[3];
sx q[3];
rz(-1.7095437) q[3];
sx q[3];
rz(0.25158238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41165274) q[2];
sx q[2];
rz(-0.81092683) q[2];
sx q[2];
rz(0.5438424) q[2];
rz(-0.020921556) q[3];
sx q[3];
rz(-1.436751) q[3];
sx q[3];
rz(-2.3232443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92886096) q[0];
sx q[0];
rz(-0.91101557) q[0];
sx q[0];
rz(2.5182356) q[0];
rz(-0.32132545) q[1];
sx q[1];
rz(-0.61505452) q[1];
sx q[1];
rz(0.26434937) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661115) q[0];
sx q[0];
rz(-0.74680416) q[0];
sx q[0];
rz(-1.2561428) q[0];
rz(-1.9064205) q[2];
sx q[2];
rz(-1.610029) q[2];
sx q[2];
rz(0.36612636) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5541861) q[1];
sx q[1];
rz(-0.79508143) q[1];
sx q[1];
rz(1.1752179) q[1];
rz(-pi) q[2];
rz(1.2530302) q[3];
sx q[3];
rz(-1.4444286) q[3];
sx q[3];
rz(1.4082091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0048206) q[2];
sx q[2];
rz(-0.68237582) q[2];
sx q[2];
rz(2.6668059) q[2];
rz(-2.7759806) q[3];
sx q[3];
rz(-0.48706278) q[3];
sx q[3];
rz(0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30597618) q[0];
sx q[0];
rz(-2.7662179) q[0];
sx q[0];
rz(2.6557652) q[0];
rz(-1.0659418) q[1];
sx q[1];
rz(-1.424788) q[1];
sx q[1];
rz(-0.33445439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7793559) q[0];
sx q[0];
rz(-1.4051172) q[0];
sx q[0];
rz(-0.30596531) q[0];
x q[1];
rz(2.7031519) q[2];
sx q[2];
rz(-1.4435569) q[2];
sx q[2];
rz(1.3194989) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.124467) q[1];
sx q[1];
rz(-1.4216058) q[1];
sx q[1];
rz(-2.6890432) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3452737) q[3];
sx q[3];
rz(-1.2346754) q[3];
sx q[3];
rz(1.4848061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36755422) q[2];
sx q[2];
rz(-2.2944821) q[2];
sx q[2];
rz(-0.96662194) q[2];
rz(-1.6537846) q[3];
sx q[3];
rz(-0.32854587) q[3];
sx q[3];
rz(0.070913471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8129355) q[0];
sx q[0];
rz(-2.2870977) q[0];
sx q[0];
rz(-0.27035126) q[0];
rz(1.6290172) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(-0.62634748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26890818) q[0];
sx q[0];
rz(-1.88694) q[0];
sx q[0];
rz(0.85319467) q[0];
rz(-2.2081991) q[2];
sx q[2];
rz(-1.2167756) q[2];
sx q[2];
rz(0.4590946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2992363) q[1];
sx q[1];
rz(-1.5299712) q[1];
sx q[1];
rz(0.090260669) q[1];
rz(-pi) q[2];
rz(-2.7553431) q[3];
sx q[3];
rz(-1.2674517) q[3];
sx q[3];
rz(1.4331499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7865929) q[2];
sx q[2];
rz(-2.2248416) q[2];
sx q[2];
rz(1.5691441) q[2];
rz(2.3908424) q[3];
sx q[3];
rz(-0.34199491) q[3];
sx q[3];
rz(-2.2641838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5759721) q[0];
sx q[0];
rz(-1.2979869) q[0];
sx q[0];
rz(-0.57814231) q[0];
rz(-1.7987953) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(-3.1411896) q[2];
sx q[2];
rz(-0.12231356) q[2];
sx q[2];
rz(-0.077153645) q[2];
rz(-2.9407637) q[3];
sx q[3];
rz(-2.8786537) q[3];
sx q[3];
rz(-2.636551) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
